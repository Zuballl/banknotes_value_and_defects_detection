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
    vars = who;
    remove = setdiff(vars, preserveVars);
    clear(remove{:}); clc; close all;
end

fprintf('--- System Wizyjnej Kontroli Banknotów (v3.1 Ograniczenie Cech) ---\n\n');

% Parametry konfiguracyjne systemu
MATCH_THRESHOLD = 0.8; 
MIN_DEFECT_AREA = 150; 
ECCENTRICITY_THRESHOLD = 0.92;
THRESHOLD_MULTIPLIER = 1.8; 
% Kryteria akceptacji jakości banknotu
ACCEPTANCE = struct( ...
    'MINOR_MAX_COUNT', 3, ...            % maksymalna liczba drobnych wad
    'MINOR_TOTAL_AREA_PCT', 0.006, ...   % dopuszczalna łączna powierzchnia drobnych wad
    'MAJOR_AREA_PCT', 0.015, ...         % próg klasyfikacji wady jako poważnej (1.5%)
    'EDGE_MARGIN', 12 ...                % margines brzegowy
);

if ~exist('MAX_FEATURES','var') || isempty(MAX_FEATURES)
    MAX_FEATURES = 10000; 
end
if ~exist('FAST_MODE','var') || isempty(FAST_MODE)
    FAST_MODE = false;
end
if FAST_MODE
    MAX_FEATURES = min(MAX_FEATURES, 3000);
    fprintf('FAST_MODE: ON (MAX_FEATURES=%d)\n', MAX_FEATURES);
end

% Główna procedura weryfikacji
try
    data = load('templateFeatures.mat');
    templateData = data.templateData;
    fprintf('Załadowano bazę %d wzorców banknotów.\n', length(templateData));
    
    % Wczytanie obrazu testowego
    testImagePath = input('Podaj ścieżkę do obrazu testowego: ', 's');

    testImageColor = imread(testImagePath);
    testImageGray = rgb2gray(testImageColor);
    
    % Identyfikacja nominału i rejestracja geometryczna
    [nominal, templateImage, templateMask, registeredImage, tform, outputView] = ...
        identifyAndAlign(testImageColor, testImageGray, templateData, MATCH_THRESHOLD, MAX_FEATURES);
    
    % Transformacja maski wzorcowej
    fprintf('Krok 3: Przekształcanie maski...\n');
    registeredMask = imwarp(templateMask, tform, 'OutputView', outputView, 'Interp', 'nearest');
    
    % Detekcja defektów i uszkodzeń
    [defectStats, ~] = detectDefects(registeredImage, templateImage, registeredMask, templateMask, MIN_DEFECT_AREA, THRESHOLD_MULTIPLIER);
    
    % Wizualizacja wyników analizy
    displayResults(registeredImage, defectStats, nominal, ECCENTRICITY_THRESHOLD, ACCEPTANCE);
    
catch ME 
    if strcmp(ME.identifier, 'MATLAB:load:couldNotFindFile')
        fprintf('\nBŁĄD: Brak pliku "templateFeatures.mat".\n');
        fprintf('Należy wcześniej wykonać skrypt "setupTemplates.m".\n');
    elseif strcmp(ME.identifier, 'MATLAB:images:imread:fileDoesNotExist')
        fprintf('\nBŁĄD: Plik obrazu testowego nie istnieje: %s\n', testImagePath);
    else
        fprintf('\nWystąpił nieoczekiwany błąd:\n');
        rethrow(ME);
    end
end


% Funkcje pomocnicze

function [nominal, templateImage, templateMask, registeredImage, tform, outputView] = ...
    identifyAndAlign(testImageColor, testImageGray, templateData, matchThreshold, maxFeatures)
    
    fprintf('Krok 1: Identyfikacja nominału...\n');

    % Detekcja punktów kluczowych na obrazie testowym
    pointsTest_all = detectORBFeatures(testImageGray);
    pointsTest = selectStrongest(pointsTest_all, maxFeatures);

    fprintf('   > Znaleziono %d wszystkich cech testowych, wybrano %d najsilniejszych.\n', ...
        pointsTest_all.Count, pointsTest.Count);

    [featuresTest, validPointsTest] = extractFeatures(testImageGray, pointsTest);

    nTemplates = length(templateData);
    matchCounts = zeros(nTemplates, 1);
    allMatchPairs = cell(nTemplates, 1);
    
    % Przygotowanie danych wzorcowych do porównania
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
    
    % Weryfikacja kompatybilności deskryptorów
    compatibleIdx = find(templateSizes == szTest);
    nIncompatible = nTemplates - length(compatibleIdx);
    
    fprintf('   > Compatible templates: %d, Incompatible: %d\n', length(compatibleIdx), nIncompatible);
    if ~isempty(compatibleIdx)
        fprintf('   > Compatible indices: ');
        fprintf('%d ', compatibleIdx);
        fprintf('\n');
    end
    
    if isempty(compatibleIdx)
        error('Brak kompatybilnych wzorców - niezgodność wymiarów deskryptorów.');
    end
    
    % Procedura dopasowywania cech
    for j = 1:length(compatibleIdx)
        i = compatibleIdx(j);
        fprintf('   > DEBUG: Matching with template %d (%s)...\n', i, templateData(i).Name);
        
        % Walidacja zgodności wymiarów przed dopasowaniem
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
        error('Niewystarczająca liczba dopasowań - niemożliwa identyfikacja nominału.');
    end
    
    nominal = templateData(bestIndex).Name;
    templateImage = templateData(bestIndex).Image; 
    templateMask = templateData(bestIndex).Mask; 
    
    fprintf('   > Rozpoznano: %s (Liczba "czystych" cech: %d)\n', nominal, maxMatches);
    fprintf('Krok 2: Wyrównywanie obrazu (rejestracja)...\n');
    
    % Estymacja transformacji geometrycznej
    bestPairs = allMatchPairs{bestIndex};
    matchedPointsTest = validPointsTest(bestPairs(:, 1));
    matchedPointsTemplate = templateData(bestIndex).Points(bestPairs(:, 2));
    
    [tform, ~] = estimateGeometricTransform2D(matchedPointsTest, matchedPointsTemplate, 'projective');
    
    outputView = imref2d(size(templateImage));
    registeredImage = imwarp(testImageColor, tform, 'OutputView', outputView);
    
    fprintf('   > Wyrównywanie zakończone.\n');
end


function [stats, defectMap] = detectDefects(registeredImage, templateImage, registeredMask, templateMask, minArea, ~)
    
    % Funkcja realizuje detekcję defektów banknotu poprzez analizę
    % geometryczną (brakujące fragmenty) oraz kolorymetryczną (plamy, napisy).
    % Wykorzystuje porównanie rzeczywistej maski banknotu z idealną maską wzorca
    % oraz analizę różnic kolorystycznych w przestrzeni LAB.
    
    fprintf('Krok 4: Wykrywanie wad (v11.0 - Brakujące części + Anomalie)...\n');

    % Parametry detekcji defektów
    MIN_MISSING_AREA = 1000;    % próg powierzchni dla brakujących fragmentów
    MIN_COLOR_ANOMALY = 3000;   % próg powierzchni dla anomalii kolorystycznych
    EDGE_MARGIN = 10;           % szerokość marginesu brzegowego
    
    % Detekcja brakujących fragmentów (zagięcia, rozdarcia)
    fprintf('   > Szukam brakujących części (zagięcia/rozdarcia)...\n');
    
    % Stwórz rzeczywistą maskę banknotu (gdzie są piksele)
    regGray = rgb2gray(registeredImage);
    realMask = regGray > 10;
    
    % Operacje morfologiczne w celu wygładzenia maski
    realMask = imclose(realMask, strel('disk', 5));
    realMask = imfill(realMask, 'holes');
    realMask = bwareaopen(realMask, 500);
    
    % Analiza porównawcza powierzchni
    areaTemplate = sum(templateMask(:));
    areaReal = sum(realMask(:));
    missingPercent = ((areaTemplate - areaReal) / areaTemplate) * 100;
    fprintf('   > Powierzchnia szablonu: %d px, rzeczywista: %d px (brakuje %.2f%%)\n', ...
        areaTemplate, areaReal, missingPercent);
    
    % Identyfikacja brakujących fragmentów poprzez różnicę logiczną masek
    missingParts = templateMask & ~realMask;
    
    % Filtracja artefaktów i eliminacja obszarów brzegowych
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
    
    % Detekcja anomalii kolorymetrycznych
    fprintf('   > Szukam anomalii kolorów (plamy, napisy)...\n');
    lab_reg = rgb2lab(imgaussfilt(registeredImage, 3));
    lab_temp = rgb2lab(imgaussfilt(templateImage, 3));
    
    % Obliczenie różnic kolorystycznych w przestrzeni LAB
    deltaE = sqrt(sum((lab_reg - lab_temp).^2, 3));
    deltaE(~realMask) = 0;
    
    % Adaptacyjne progowanie różnic kolorystycznych
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
    
    % Agregacja wykrytych defektów
    fprintf('   > Łączę wykryte wady...\n');
    defectMap = missingParts | colorAnomalies;

    % Ekstrakcja cech i klasyfikacja defektów
    cc = bwconncomp(defectMap);
    rawStats = regionprops(cc, 'Area', 'BoundingBox', 'Eccentricity', 'PixelIdxList');

    % klasyfikacja typu (rozdarcie vs plama) i wagi (minor/major)
    stats = struct('Area', {}, 'BoundingBox', {}, 'Eccentricity', {}, 'Type', {}, 'Severity', {}, 'AreaPct', {});
    for i = 1:numel(rawStats)
        s = rawStats(i);
        areaPct = double(s.Area) / double(areaTemplate);

        % Determinacja typu defektu na podstawie składowych
        pix = s.PixelIdxList;
        fracMissing = sum(missingParts(pix)) / numel(pix);
        fracColor = sum(colorAnomalies(pix)) / numel(pix);
        if fracMissing >= fracColor
            typeStr = 'Rozdarcie / Zagniecenie';
        else
            typeStr = 'Plama / Napis';
        end

        % Ocena ciężkości defektu zgodnie z ustalonymi kryteriami
        isEdgeTouch = false;
        bb = s.BoundingBox;
        x1 = round(max(1, bb(1))); y1 = round(max(1, bb(2)));
        x2 = round(min(size(missingParts,2), bb(1)+bb(3)));
        y2 = round(min(size(missingParts,1), bb(2)+bb(4)));
        if x1 <= EDGE_MARGIN || y1 <= EDGE_MARGIN || ...
           (size(missingParts,2) - x2) <= EDGE_MARGIN || (size(missingParts,1) - y2) <= EDGE_MARGIN
            isEdgeTouch = true;
        end

        majorByArea = areaPct >= 0.015;
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
    % Funkcja realizuje wizualizację wyników analizy z zaznaczeniem wykrytych defektów
    fprintf('Krok 5: Wyświetlanie wyników (tylko poważne wady)...\n');

    figure;
    imshow(image);
    hold on;

    % Przypadek braku wykrytych defektów
    if isempty(stats)
        titleStr = sprintf('Nominał: %s. Status: ZAAKCEPTOWANY', nominal);
        fprintf('\n--- WERDYKT: BANKNOT ZAAKCEPTOWANY (brak wad) ---\n');
        title(titleStr, 'FontSize', 14); hold off; return;
    end

    % Selekcja defektów o statusie poważnym
    isMajor = arrayfun(@(s) isfield(s,'Severity') && strcmpi(s.Severity,'major'), stats);
    majorStats = stats(isMajor);
    nMajor = numel(majorStats);

    % Wizualizacja poważnych defektów
    for i = 1:length(majorStats)
        bb = majorStats(i).BoundingBox;
        rectangle('Position', bb, 'EdgeColor', 'r', 'LineWidth', 2);

        % Generacja etykiety tekstowej
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

    % Procedura decyzyjna akceptacji banknotu
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

% Funkcja pomocnicza operatora trójargumentowego
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end