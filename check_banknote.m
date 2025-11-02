% --- check_banknote.m (Wersja 3.1 - Ograniczenie Cech) ---
%
% Cel: Ograniczenie liczby cech na obrazie testowym.
% Wymagania: Nowy plik "templateFeatures.mat" (stworzony przez v3.1)

% Allow calling code to set `testImagePath`, `MAX_FEATURES` or `FAST_MODE` before running.
% If not provided, script will run interactively as before.
% Note: script still clears most variables to keep behavior predictable.
% Keep any externally provided configuration variables.
preserveVars = {'testImagePath','MAX_FEATURES','FAST_MODE'};
keep = false;
for k = 1:numel(preserveVars)
    if exist(preserveVars{k}, 'var')
        keep = true; break;
    end
end
if ~keep
    clear; clc; close all;
else
    % remove everything except preserved variables
    vars = who;
    remove = setdiff(vars, preserveVars);
    clear(remove{:}); clc; close all;
end

fprintf('--- System Wizyjnej Kontroli Banknotów (v3.1 Ograniczenie Cech) ---\n\n');

% --- Stałe Konfiguracyjne ---
MATCH_THRESHOLD = 0.8; 
MIN_DEFECT_AREA = 150; 
ECCENTRICITY_THRESHOLD = 0.92;
THRESHOLD_MULTIPLIER = 1.8; 
% Reguły akceptacji (tolerancja drobnych wad)
ACCEPTANCE = struct( ...
    'MINOR_MAX_COUNT', 3, ...            % maks. liczba drobnych wad
    'MINOR_TOTAL_AREA_PCT', 0.006, ...   % łączna powierzchnia drobnych wad (0.6%)
    'MAJOR_AREA_PCT', 0.015, ...         % pojedyncza wada >1.5% = poważna
    'EDGE_MARGIN', 12 ...                % margines przy krawędzi dla rozdarć
);
% Allow overriding MAX_FEATURES and FAST_MODE from the caller to speed up runs
if ~exist('MAX_FEATURES','var') || isempty(MAX_FEATURES)
    MAX_FEATURES = 10000; % Default limit for test image features
end
if ~exist('FAST_MODE','var') || isempty(FAST_MODE)
    FAST_MODE = false; % When true, reduce work (fewer features)
end
if FAST_MODE
    MAX_FEATURES = min(MAX_FEATURES, 3000);
    fprintf('FAST_MODE: ON (MAX_FEATURES=%d)\n', MAX_FEATURES);
end

% --- Główna Logika ---
try
    data = load('templateFeatures.mat');
    templateData = data.templateData;
    fprintf('Załadowano bazę %d wzorców banknotów.\n', length(templateData));
    
    % Zawsze pytaj o ścieżkę (ignoruj wcześniej ustawioną zmienną)
    testImagePath = input('Podaj ścieżkę do obrazu testowego: ', 's');

    testImageColor = imread(testImagePath);
    testImageGray = rgb2gray(testImageColor);
    
    % --- Krok 1 i 2: Identyfikacja i Wyrównanie ---
    [nominal, templateImage, templateMask, registeredImage, tform, outputView] = ...
        identifyAndAlign(testImageColor, testImageGray, templateData, MATCH_THRESHOLD, MAX_FEATURES);
    
    % --- Krok 3: Przekształcanie Maski ---
    fprintf('Krok 3: Przekształcanie maski...\n');
    registeredMask = imwarp(templateMask, tform, 'OutputView', outputView, 'Interp', 'nearest');
    
    % --- Krok 4: Detekcja Wad ---
    [defectStats, ~] = detectDefects(registeredImage, templateImage, registeredMask, templateMask, MIN_DEFECT_AREA, THRESHOLD_MULTIPLIER);
    
    % --- Krok 5: Prezentacja Wyników ---
    displayResults(registeredImage, defectStats, nominal, ECCENTRICITY_THRESHOLD, ACCEPTANCE);
    
catch ME 
    if strcmp(ME.identifier, 'MATLAB:load:couldNotFindFile')
        fprintf('\nBŁĄD KRYTYCZNY: Nie znaleziono pliku "templateFeatures.mat"!\n');
        fprintf('Uruchom najpierw "setupTemplates.m" (wersję 3.1).\n');
    elseif strcmp(ME.identifier, 'MATLAB:images:imread:fileDoesNotExist')
        fprintf('\nBŁĄD: Plik obrazu testowego nie istnieje: %s\n', testImagePath);
    else
        fprintf('\nWystąpił nieoczekiwany błąd:\n');
        rethrow(ME);
    end
end


% --- Funkcje Pomocnicze (Lokalne) ---

function [nominal, templateImage, templateMask, registeredImage, tform, outputView] = ...
    identifyAndAlign(testImageColor, testImageGray, templateData, matchThreshold, maxFeatures)
    
    fprintf('Krok 1: Identyfikacja nominału...\n');

    % Detect features but cap the number of strongest to MAX_FEATURES for speed
    pointsTest_all = detectORBFeatures(testImageGray);
    pointsTest = selectStrongest(pointsTest_all, maxFeatures);

    fprintf('   > Znaleziono %d wszystkich cech testowych, wybrano %d najsilniejszych.\n', ...
        pointsTest_all.Count, pointsTest.Count);

    [featuresTest, validPointsTest] = extractFeatures(testImageGray, pointsTest);

    nTemplates = length(templateData);
    matchCounts = zeros(nTemplates, 1);
    allMatchPairs = cell(nTemplates, 1);
    
    % Extract all template features and sizes into simple arrays
    szTest = size(featuresTest, 2);
    fprintf('   > DEBUG: Test descriptor size = %d columns\n', szTest);
    
    templateFeatures = cell(nTemplates, 1);
    templateSizes = zeros(nTemplates, 1);
    
    for i = 1:nTemplates
        templateFeatures{i} = templateData(i).Features;
        templateSizes(i) = size(templateData(i).Features, 2);
        fprintf('   > DEBUG: Template %d (%s) descriptor size = %d columns\n', ...
            i, templateData(i).Name, templateSizes(i));
    end
    
    % Check compatibility and filter to only compatible templates
    compatibleIdx = find(templateSizes == szTest);
    nIncompatible = nTemplates - length(compatibleIdx);
    
    fprintf('   > Compatible templates: %d, Incompatible: %d\n', length(compatibleIdx), nIncompatible);
    if ~isempty(compatibleIdx)
        fprintf('   > Compatible indices: ');
        fprintf('%d ', compatibleIdx);
        fprintf('\n');
    end
    
    if isempty(compatibleIdx)
        error('BRAK KOMPATYBILNYCH SZABLONÓW! Wszystkie szablony mają inny rozmiar deskryptorów niż obraz testowy.');
    end
    
    % Match ONLY against compatible templates (serial processing for reliability)
    for j = 1:length(compatibleIdx)
        i = compatibleIdx(j);
        fprintf('   > DEBUG: Matching with template %d (%s)...\n', i, templateData(i).Name);
        
        % Double-check before matching
        sz1 = size(featuresTest, 2);
        sz2 = size(templateFeatures{i}, 2);
        if sz1 ~= sz2
            fprintf('   > ERROR: Size mismatch at runtime! Test=%d, Template=%d\n', sz1, sz2);
            continue;
        end
        
        indexPairs = matchFeatures(featuresTest, templateFeatures{i}, ...
            'MatchThreshold', 50, 'MaxRatio', matchThreshold, 'Unique', true);
        matchCounts(i) = size(indexPairs, 1);
        allMatchPairs{i} = indexPairs;
        fprintf('   > DEBUG: Found %d matches\n', matchCounts(i));
    end
    
    [maxMatches, bestIndex] = max(matchCounts);
    
    if maxMatches < 10 
        error('Nie udało się rozpoznać nominału. Zbyt niska jakość obrazu lub brak wzorca.');
    end
    
    nominal = templateData(bestIndex).Name;
    templateImage = templateData(bestIndex).Image; 
    templateMask = templateData(bestIndex).Mask; 
    
    fprintf('   > Rozpoznano: %s (Liczba "czystych" cech: %d)\n', nominal, maxMatches);
    fprintf('Krok 2: Wyrównywanie obrazu (rejestracja)...\n');
    
    bestPairs = allMatchPairs{bestIndex};
    matchedPointsTest = validPointsTest(bestPairs(:, 1));
    matchedPointsTemplate = templateData(bestIndex).Points(bestPairs(:, 2));
    
    [tform, ~] = estimateGeometricTransform2D(matchedPointsTest, matchedPointsTemplate, 'projective');
    
    outputView = imref2d(size(templateImage));
    registeredImage = imwarp(testImageColor, tform, 'OutputView', outputView);
    
    fprintf('   > Wyrównywanie zakończone.\n');
end




% --- Wklej to w miejsce STAREJ funkcji detectDefects ---

function [stats, defectMap] = detectDefects(registeredImage, templateImage, registeredMask, templateMask, minArea, ~)
    
    % ============================================================================
    % PODEJŚCIE V11.0: DETEKCJA BRAKUJĄCYCH CZĘŚCI + ANOMALIE WEWNĘTRZNE
    % 
    % KLUCZOWE SPOSTRZEŻENIE:
    % - Zagięcia/rozdarcia = BRAKUJĄCE CZĘŚCI geometrii banknotu
    % - Po rejestracji mamy MASKĘ SZABLONU = idealny kształt
    % - Możemy wykryć gdzie banknut NIE WYPEŁNIA idealnego kształtu!
    %
    % STRATEGIA:
    % 1. Znajdź rzeczywistą maskę banknotu (gdzie są piksele)
    % 2. Porównaj z maską szablonu (gdzie POWINNY być piksele)  
    % 3. Różnica = brakujące części (zagięcia, rozdarcia, dziury)
    % 4. Dodatkowo: anomalie kolorów wewnątrz (plamy, napisy)
    % ============================================================================
    
    fprintf('Krok 4: Wykrywanie wad (v11.0 - Brakujące części + Anomalie)...\n');

    % Parametry
    MIN_MISSING_AREA = 1000;    % Minimalny obszar brakującej części
    MIN_COLOR_ANOMALY = 3000;   % Minimalny obszar anomalii kolorów
    EDGE_MARGIN = 10;           % Margines brzegów
    
    % === 1. WYKRYJ BRAKUJĄCE CZĘŚCI BANKNOTU ===
    fprintf('   > Szukam brakujących części (zagięcia/rozdarcia)...\n');
    
    % Stwórz rzeczywistą maskę banknotu (gdzie są piksele)
    regGray = rgb2gray(registeredImage);
    realMask = regGray > 10; % proste progowanie - gdzie jest cokolwiek
    
    % Wygładź i wypełnij małe dziury
    realMask = imclose(realMask, strel('disk', 5));
    realMask = imfill(realMask, 'holes');
    realMask = bwareaopen(realMask, 500);
    
    % DEBUG: Porównaj powierzchnie
    areaTemplate = sum(templateMask(:));
    areaReal = sum(realMask(:));
    missingPercent = ((areaTemplate - areaReal) / areaTemplate) * 100;
    fprintf('   > Powierzchnia szablonu: %d px, rzeczywista: %d px (brakuje %.2f%%)\n', ...
        areaTemplate, areaReal, missingPercent);
    
    % Gdzie POWINIEN być banknut (maska szablonu) vs gdzie JEST (realMask)
    % Różnica = brakujące części!
    missingParts = templateMask & ~realMask;
    
    % Usuń małe artefakty i brzegi
    missingParts = bwareaopen(missingParts, MIN_MISSING_AREA);
    [h, w] = size(missingParts);
    edgeMask = true(h, w);
    edgeMask(1:EDGE_MARGIN, :) = false;
    edgeMask(end-EDGE_MARGIN+1:end, :) = false;
    edgeMask(:, 1:EDGE_MARGIN) = false;
    edgeMask(:, end-EDGE_MARGIN+1:end) = false;
    missingParts = missingParts & edgeMask;
    
    nMissing = bwconncomp(missingParts).NumObjects;
    fprintf('   > Znaleziono %d brakujących części\n', nMissing);
    
    % === 2. WYKRYJ ANOMALIE KOLORÓW (plamy, napisy) ===
    fprintf('   > Szukam anomalii kolorów (plamy, napisy)...\n');
    lab_reg = rgb2lab(imgaussfilt(registeredImage, 3));
    lab_temp = rgb2lab(imgaussfilt(templateImage, 3));
    
    % Delta E w LAB - TYLKO wewnątrz rzeczywistej maski
    deltaE = sqrt(sum((lab_reg - lab_temp).^2, 3));
    deltaE(~realMask) = 0;
    
    % BARDZO WYSOKIE progowanie (tylko duże różnice)
    if any(deltaE(:) > 0)
        threshold = mean(deltaE(deltaE > 0)) + 3*std(deltaE(deltaE > 0));
    else
        threshold = inf;
    end
    colorAnomalies = deltaE > threshold;
    colorAnomalies = bwareaopen(colorAnomalies, MIN_COLOR_ANOMALY);
    colorAnomalies = imclose(colorAnomalies, strel('disk', 5));
    
    nColor = bwconncomp(colorAnomalies).NumObjects;
    fprintf('   > Znaleziono %d anomalii kolorów\n', nColor);
    
    % === 3. POŁĄCZ BRAKUJĄCE CZĘŚCI + ANOMALIE KOLORÓW ===
    fprintf('   > Łączę wykryte wady...\n');
    defectMap = missingParts | colorAnomalies;

    % === 4. ANALIZA WYKRYTYCH REGIONÓW + KLASYFIKACJA ===
    cc = bwconncomp(defectMap);
    rawStats = regionprops(cc, 'Area', 'BoundingBox', 'Eccentricity', 'PixelIdxList');

    % klasyfikacja typu (rozdarcie vs plama) i wagi (minor/major)
    stats = struct('Area', {}, 'BoundingBox', {}, 'Eccentricity', {}, 'Type', {}, 'Severity', {}, 'AreaPct', {});
    for i = 1:numel(rawStats)
        s = rawStats(i);
        areaPct = double(s.Area) / double(areaTemplate);

        % udział pikseli z map: która dominuje w regionie?
        pix = s.PixelIdxList;
        fracMissing = sum(missingParts(pix)) / numel(pix);
        fracColor = sum(colorAnomalies(pix)) / numel(pix);
        if fracMissing >= fracColor
            typeStr = 'Rozdarcie / Zagniecenie';
        else
            typeStr = 'Plama / Napis';
        end

        % heurystyka powagi uszkodzenia
        isEdgeTouch = false;
        bb = s.BoundingBox;
        x1 = round(max(1, bb(1))); y1 = round(max(1, bb(2)));
        x2 = round(min(size(missingParts,2), bb(1)+bb(3)));
        y2 = round(min(size(missingParts,1), bb(2)+bb(4)));
        if x1 <= EDGE_MARGIN || y1 <= EDGE_MARGIN || ...
           (size(missingParts,2) - x2) <= EDGE_MARGIN || (size(missingParts,1) - y2) <= EDGE_MARGIN
            isEdgeTouch = true;
        end

        majorByArea = areaPct >= 0.015;         % 1.5% obszaru
        majorEdgeTear = (s.Eccentricity > 0.92) && isEdgeTouch && (areaPct >= 0.001);
        if majorByArea || majorEdgeTear
            severityStr = 'major';
        else
            severityStr = 'minor';
        end

        stats(end+1) = struct( ...
            'Area', s.Area, ...
            'BoundingBox', s.BoundingBox, ...
            'Eccentricity', s.Eccentricity, ...
            'Type', typeStr, ...
            'Severity', severityStr, ...
            'AreaPct', areaPct ...
        ); %#ok<AGROW>
    end

    fprintf('   > Analiza zakończona. Znaleziono %d wad.\n', length(stats));
end
  
    


   







function displayResults(image, stats, nominal, eccentricityThreshold, acceptance)
    % Zmienione: pokazuj tylko poważne (major) wady. Logika drobnych wad została usunięta.
    fprintf('Krok 5: Wyświetlanie wyników (tylko poważne wady)...\n');

    figure;
    imshow(image);
    hold on;

    % jeśli nie ma żadnych wykrytych obszarów => akceptuj
    if isempty(stats)
        titleStr = sprintf('Nominał: %s. Status: ZAAKCEPTOWANY', nominal);
        fprintf('\n--- WERDYKT: BANKNOT ZAAKCEPTOWANY (brak wad) ---\n');
        title(titleStr, 'FontSize', 14); hold off; return;
    end

    % Wyfiltruj tylko poważne wady
    isMajor = arrayfun(@(s) isfield(s,'Severity') && strcmpi(s.Severity,'major'), stats);
    majorStats = stats(isMajor);
    nMajor = numel(majorStats);

    % Rysuj tylko poważne wady (czerwone obramowanie)
    for i = 1:length(majorStats)
        bb = majorStats(i).BoundingBox;
        rectangle('Position', bb, 'EdgeColor', 'r', 'LineWidth', 2);

        % etykieta (tylko major)
        if isfield(majorStats(i), 'Type')
            label = majorStats(i).Type;
        else
            if majorStats(i).Eccentricity > eccentricityThreshold
                label = 'Rozdarcie / Zagniecenie';
            else
                label = 'Plama / Napis';
            end
        end
        label = [label ' (poważna)'];
        text(bb(1), bb(2) - 10, label, 'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold');
        fprintf('   > Wada #%d: Typ: %s, Powierzchnia: %d px (%.3f%%), MAJOR\n', ...
            i, label, majorStats(i).Area, 100*majorStats(i).AreaPct);
    end

    % Decyzja: akceptujemy jeśli brak poważnych wad
    if nMajor == 0
        titleStr = sprintf('Nominał: %s. Status: ZAAKCEPTOWANY (brak poważnych wad)', nominal);
        fprintf('\n--- WERDYKT: BANKNOT ZAAKCEPTOWANY (brak poważnych wad) ---\n');
    else
        titleStr = sprintf('Nominał: %s. Status: DO SPRAWDZENIA (poważne: %d)', nominal, nMajor);
        fprintf('\n--- WERDYKT: BANKNOT DO SPRAWDZENIA ---\n');
    end
    title(titleStr, 'FontSize', 14);
    hold off;

    title(titleStr, 'FontSize', 14);
    hold off;
end

% mały pomocnik do czytelnego fprintf
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end