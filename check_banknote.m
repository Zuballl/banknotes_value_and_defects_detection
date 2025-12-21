% check_banknote.m
% Główny skrypt systemu wizyjnej kontroli banknotów.
% Rozpoznaje nominał banknotu i wykrywa defekty (plamy, rozdarcia, zagięcia).
% 
% Działanie:
% 1. Wczytuje obraz testowy
% 2. Identyfikuje nominał przez dopasowanie cech ORB do wzorców
% 3. Wyrównuje obraz testowy do wzorca (rejestracja geometryczna)
% 4. Wykrywa defekty przez analizę geometryczną i kolorystyczną
% 5. Wyświetla wyniki z zaznaczeniem poważnych wad

% Inteligentne czyszczenie workspace - zachowuje wybrane zmienne jeśli istnieją
% To pozwala na szybsze testowanie bez ponownego ustawiania parametrów
preserveVars = {'testImagePath','MAX_FEATURES','FAST_MODE'};
keep = false;
for k = 1:numel(preserveVars)
    if exist(preserveVars{k}, 'var')
        keep = true; break;
    end
end
if ~keep
    clear; clc; close all;  % standardowe czyszczenie
else
    vars = who;  % lista wszystkich zmiennych
    remove = setdiff(vars, preserveVars);  % zmienne do usunięcia
    clear(remove{:}); clc; close all;  % usuń tylko niepotrzebne
end

fprintf('--- System Wizyjnej Kontroli Banknotów (v3.1 Ograniczenie Cech) ---\n\n');

% Parametry konfiguracyjne systemu
% Te wartości można dostrajać w zależności od potrzeb
MATCH_THRESHOLD = 0.8;  % próg dopasowania cech (0-1, wyżej = bardziej restrykcyjne)
MIN_DEFECT_AREA = 150;  % minimalna powierzchnia defektu w pikselach
ECCENTRICITY_THRESHOLD = 0.92;  % próg mimośrodowości dla rozdarć (0-1, bliżej 1 = bardziej wydłużone)
THRESHOLD_MULTIPLIER = 1.8;  % mnożnik progu dla detekcji różnic 
% Kryteria akceptacji jakości banknotu
% Struktura definiująca progi klasyfikacji - można je modyfikować
% w zależności od tego jak rygorystyczny ma być system
ACCEPTANCE = struct( ...
    'MINOR_MAX_COUNT', 3, ...            % maksymalna liczba drobnych wad (minor defects)
    'MINOR_TOTAL_AREA_PCT', 0.006, ...   % dopuszczalna łączna powierzchnia drobnych wad (0.6%)
    'MAJOR_AREA_PCT', 0.015, ...         % próg klasyfikacji wady jako poważnej - 1.5% powierzchni
    'EDGE_MARGIN', 12 ...                % margines brzegowy w pikselach - ignorujemy defekty przy krawędziach
);

% Parametr MAX_FEATURES kontroluje ile cech maksymalnie analizujemy
% Więcej cech = dokładniejsze rozpoznawanie, ale wolniejsze działanie
if ~exist('MAX_FEATURES','var') || isempty(MAX_FEATURES)
    MAX_FEATURES = 10000;  % domyślnie analizujemy do 10000 cech
end

% FAST_MODE - tryb szybki, użyteczny podczas testowania
% Ogranicza liczbę cech dla przyspieszenia obliczeń
if ~exist('FAST_MODE','var') || isempty(FAST_MODE)
    FAST_MODE = false;  % domyślnie wyłączony
end
if FAST_MODE
    MAX_FEATURES = min(MAX_FEATURES, 3000);  % w trybie szybkim max 3000 cech
    fprintf('FAST_MODE: ON (MAX_FEATURES=%d)\n', MAX_FEATURES);
end

% Główna procedura weryfikacji banknotu
% Obsługuje błędy poprzez try-catch
try
    % Wczytujemy bazę wzorców przygotowaną wcześniej przez setupTemplates.m
    data = load('templateFeatures.mat');
    templateData = data.templateData;
    fprintf('Załadowano bazę %d wzorców banknotów.\n', length(templateData));
    
    % Wczytanie obrazu testowego - użytkownik podaje ścieżkę
    testImagePath = input('Podaj ścieżkę do obrazu testowego: ', 's');

    % Ładujemy obraz w kolorze i konwertujemy do skali szarości
    testImageColor = imread(testImagePath);
    testImageGray = rgb2gray(testImageColor);
    
    % Identyfikacja nominału i rejestracja geometryczna
    % Ta funkcja rozpoznaje jaki to banknot i wyrównuje go do wzorca
    [nominal, templateImage, templateMask, registeredImage, tform, outputView] = ...
        identifyAndAlign(testImageColor, testImageGray, templateData, MATCH_THRESHOLD, MAX_FEATURES);
    
    % Transformacja maski wzorcowej
    % Musimy też przekształcić maskę wzorca tak samo jak obraz testowy,
    % żeby móc porównać co jest bankotem a co tłem
    fprintf('Krok 3: Przekształcanie maski...\n');
    registeredMask = imwarp(templateMask, tform, 'OutputView', outputView, 'Interp', 'nearest');
    
    % Detekcja defektów i uszkodzeń
    % Porównujemy wyrównany banknot z wzorcem i szukamy różnic
    [defectStats, ~] = detectDefects(registeredImage, templateImage, registeredMask, templateMask, MIN_DEFECT_AREA, THRESHOLD_MULTIPLIER);
    
    % Wizualizacja wyników analizy
    % Wyświetlamy obraz z zaznaczonymi defektami i werdyktem
    displayResults(registeredImage, defectStats, nominal, ECCENTRICITY_THRESHOLD, ACCEPTANCE);
    
catch ME  % ME = MException (obiekt wyjątku)
    % Obsługa różnych typów błędów z przyjaznymi komunikatami
    if strcmp(ME.identifier, 'MATLAB:load:couldNotFindFile')
        fprintf('\nBŁĄD: Brak pliku "templateFeatures.mat".\n');
        fprintf('Należy wcześniej wykonać skrypt "setupTemplates.m".\n');
    elseif strcmp(ME.identifier, 'MATLAB:images:imread:fileDoesNotExist')
        fprintf('\nBŁĄD: Plik obrazu testowego nie istnieje: %s\n', testImagePath);
    else
        % Jeśli to jakiś inny błąd, po prostu go wyrzucamy dalej
        fprintf('\nWystąpił nieoczekiwany błąd:\n');
        rethrow(ME);
    end
end


%% ========================================================================
% FUNKCJE POMOCNICZE
% ==========================================================================

%% Funkcja identifyAndAlign - rozpoznaje nominał i wyrównuje obraz
% 
% Argumenty wejściowe:
%   testImageColor - kolorowy obraz testowy
%   testImageGray - obraz testowy w skali szarości
%   templateData - struktura z danymi wzorców
%   matchThreshold - próg dopasowania cech
%   maxFeatures - maksymalna liczba cech do analizy
%
% Argumenty wyjściowe:
%   nominal - rozpoznany nominał (np. "10 PLN #1")
%   templateImage - obraz wzorca
%   templateMask - maska wzorca
%   registeredImage - wyrównany obraz testowy
%   tform - transformacja geometryczna
%   outputView - parametry widoku wyjściowego
function [nominal, templateImage, templateMask, registeredImage, tform, outputView] = ...
    identifyAndAlign(testImageColor, testImageGray, templateData, matchThreshold, maxFeatures)
    
    fprintf('Krok 1: Identyfikacja nominału...\n');

    % Detekcja punktów kluczowych na obrazie testowym algorytmem ORB
    % Wykrywamy wszystkie możliwe punkty charakterystyczne
    pointsTest_all = detectORBFeatures(testImageGray);
    % Wybieramy tylko najsilniejsze - lepiej mieć mniej ale lepszych cech
    pointsTest = selectStrongest(pointsTest_all, maxFeatures);

    fprintf('   > Znaleziono %d wszystkich cech testowych, wybrano %d najsilniejszych.\n', ...
        pointsTest_all.Count, pointsTest.Count);

    % Ekstrakcja deskryptorów - wektorów opisujących każdy punkt
    % validPointsTest to punkty dla których udało się wyekstrahować deskryptory
    [featuresTest, validPointsTest] = extractFeatures(testImageGray, pointsTest);

    % Przygotowanie do porównania z wzorcami
    nTemplates = length(templateData);
    matchCounts = zeros(nTemplates, 1);  % tablica na liczbę dopasowań dla każdego wzorca
    allMatchPairs = cell(nTemplates, 1);  % pary dopasowanych punktów
    
    % Przygotowanie danych wzorcowych do porównania
    % Sprawdzamy wymiary deskryptorów - muszą się zgadzać
    szTest = size(featuresTest, 2);
    fprintf('   > DEBUG: Test descriptor size = %d columns\n', szTest);
    
    templateFeatures = cell(nTemplates, 1);
    templateSizes = zeros(nTemplates, 1);
    
    % Zbieramy informacje o wymiarach deskryptorów wszystkich wzorców
    for i = 1:nTemplates
        templateFeatures{i} = templateData(i).Features;
        templateSizes(i) = size(templateData(i).Features, 2);
        fprintf('   > DEBUG: Template %d (%s) descriptor size = %d columns\n', ...
            i, templateData(i).Name, templateSizes(i));
    end
    
    % Weryfikacja kompatybilności deskryptorów
    % Dopasowanie działa tylko jeśli wymiary się zgadzają
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
    
    % Procedura dopasowywania cech - szukamy najlepszego wzorca
    % Dla każdego kompatybilnego wzorca liczymy ile cech się dopasowało
    for j = 1:length(compatibleIdx)
        i = compatibleIdx(j);
        fprintf('   > DEBUG: Matching with template %d (%s)...\n', i, templateData(i).Name);
        
        % Walidacja zgodności wymiarów przed dopasowaniem (podwójna kontrola)
        sz1 = size(featuresTest, 2);
        sz2 = size(templateFeatures{i}, 2);
        if sz1 ~= sz2
            fprintf('   > ERROR: Size mismatch at runtime! Test=%d, Template=%d\n', sz1, sz2);
            continue;
        end
        
        % Dopasowanie cech - funkcja matchFeatures znajduje które punkty
        % z obrazu testowego odpowiadają punktom ze wzorca
        indexPairs = matchFeatures(featuresTest, templateFeatures{i}, ...
            'MatchThreshold', 50, 'MaxRatio', matchThreshold, 'Unique', true);
        matchCounts(i) = size(indexPairs, 1);  % ile się dopasowało
        allMatchPairs{i} = indexPairs;
        fprintf('   > DEBUG: Found %d matches\n', matchCounts(i));
    end
    
    % Wybieramy wzorzec z największą liczbą dopasowań
    [maxMatches, bestIndex] = max(matchCounts);
    
    % Jeśli mamy za mało dopasowań, nie możemy pewnie rozpoznać nominału
    if maxMatches < 10 
        error('Niewystarczająca liczba dopasowań - niemożliwa identyfikacja nominału.');
    end
    
    % Wyciągamy informacje o najlepiej dopasowanym wzorcu
    nominal = templateData(bestIndex).Name;
    templateImage = templateData(bestIndex).Image; 
    templateMask = templateData(bestIndex).Mask; 
    
    fprintf('   > Rozpoznano: %s (Liczba "czystych" cech: %d)\n', nominal, maxMatches);
    fprintf('Krok 2: Wyrównywanie obrazu (rejestracja)...\n');
    
    % Estymacja transformacji geometrycznej
    % Obliczamy jak przekształcić obraz testowy, żeby wyglądał tak jak wzorzec
    bestPairs = allMatchPairs{bestIndex};
    matchedPointsTest = validPointsTest(bestPairs(:, 1));  % punkty z obrazu testowego
    matchedPointsTemplate = templateData(bestIndex).Points(bestPairs(:, 2));  % odpowiadające punkty ze wzorca
    
    % Obliczamy transformację projekcyjną (uwzględnia perspektywę)
    [tform, ~] = estimateGeometricTransform2D(matchedPointsTest, matchedPointsTemplate, 'projective');
    
    % Stosujemy transformację do obrazu testowego
    outputView = imref2d(size(templateImage));  % ustawiamy rozmiar wyjścia jak wzorzec
    registeredImage = imwarp(testImageColor, tform, 'OutputView', outputView);
    
    fprintf('   > Wyrównywanie zakończone.\n');
end


%% Funkcja detectDefects - wykrywa defekty na banknocie
%
% Funkcja realizuje detekcję defektów banknotu poprzez analizę
% geometryczną (brakujące fragmenty) oraz kolorymetryczną (plamy, napisy).
% Wykorzystuje porównanie rzeczywistej maski banknotu z idealną maską wzorca
% oraz analizę różnic kolorystycznych w przestrzeni LAB.
%
% Argumenty wejściowe:
%   registeredImage - wyrównany obraz testowy
%   templateImage - obraz wzorca
%   registeredMask - wyrównana maska testowa
%   templateMask - maska wzorca
%   minArea - minimalna powierzchnia defektu
%
% Argumenty wyjściowe:
%   stats - struktura z informacjami o wykrytych defektach
%   defectMap - mapa binarna wszystkich defektów
function [stats, defectMap] = detectDefects(registeredImage, templateImage, registeredMask, templateMask, minArea, ~)
    
    % Funkcja realizuje detekcję defektów banknotu poprzez analizę
    % geometryczną (brakujące fragmenty) oraz kolorymetryczną (plamy, napisy).
    % Wykorzystuje porównanie rzeczywistej maski banknotu z idealną maską wzorca
    % oraz analizę różnic kolorystycznych w przestrzeni LAB.
    
    fprintf('Krok 4: Wykrywanie wad (v11.0 - Brakujące części + Anomalie)...\n');

    % Parametry detekcji defektów - można je dostrajać
    MIN_MISSING_AREA = 1000;    % próg powierzchni dla brakujących fragmentów (w pikselach)
    MIN_COLOR_ANOMALY = 3000;   % próg powierzchni dla anomalii kolorystycznych (w pikselach)
    EDGE_MARGIN = 10;           % szerokość marginesu brzegowego (ignorujemy brzegi)
    
    % Detekcja brakujących fragmentów (zagięcia, rozdarcia)
    fprintf('   > Szukam brakujących części (zagięcia/rozdarcia)...\n');
    
    % Tworzymy rzeczywistą maskę banknotu - gdzie faktycznie są piksele
    % Piksele ciemniejsze niż 10 uznajemy za tło
    regGray = rgb2gray(registeredImage);
    realMask = regGray > 10;
    
    % Operacje morfologiczne w celu wygładzenia maski
    % Zamykanie (close) - łączy małe przerwy
    realMask = imclose(realMask, strel('disk', 5));
    % Wypełnianie dziur
    realMask = imfill(realMask, 'holes');
    % Usuwanie bardzo małych obiektów (szum)
    realMask = bwareaopen(realMask, 500);
    
    % Analiza porównawcza powierzchni - ile % banknotu brakuje
    areaTemplate = sum(templateMask(:));  % powierzchnia idealnego wzorca
    areaReal = sum(realMask(:));  % powierzchnia rzeczywistego banknotu
    missingPercent = ((areaTemplate - areaReal) / areaTemplate) * 100;
    fprintf('   > Powierzchnia szablonu: %d px, rzeczywista: %d px (brakuje %.2f%%)\n', ...
        areaTemplate, areaReal, missingPercent);
    
    % Identyfikacja brakujących fragmentów poprzez różnicę logiczną masek
    % missingParts = tam gdzie jest wzorzec ALE nie ma rzeczywistości
    missingParts = templateMask & ~realMask;
    
    % Filtracja artefaktów - usuwamy za małe obiekty
    missingParts = bwareaopen(missingParts, MIN_MISSING_AREA);
    
    % Eliminacja obszarów brzegowych - przy krawędziach często są błędy rejestracji
    [h, w] = size(missingParts);
    edgeMask = true(h, w);  % zaczynamy od pełnej maski
    % Zerujemy brzegi - góra, dół, lewo, prawo
    edgeMask(1:EDGE_MARGIN, :) = false;
    edgeMask(end-EDGE_MARGIN+1:end, :) = false;
    edgeMask(:, 1:EDGE_MARGIN) = false;
    edgeMask(:, end-EDGE_MARGIN+1:end) = false;
    missingParts = missingParts & edgeMask;  % stosujemy maskę brzegową
    
    nMissing = bwconncomp(missingParts).NumObjects;  % ile znaleźliśmy brakujących fragmentów
    fprintf('   > Znaleziono %d brakujących części\n', nMissing);
    
    % Detekcja anomalii kolorymetrycznych (plamy, napisy)
    fprintf('   > Szukam anomalii kolorów (plamy, napisy)...\n');
    
    % Konwersja do przestrzeni LAB - lepsza dla porównań kolorystycznych
    % LAB lepiej oddaje różnice widoczne dla ludzkiego oka niż RGB
    % imgaussfilt wygładza obraz - redukuje szum
    lab_reg = rgb2lab(imgaussfilt(registeredImage, 3));
    lab_temp = rgb2lab(imgaussfilt(templateImage, 3));
    
    % Obliczenie różnic kolorystycznych w przestrzeni LAB
    % deltaE (Delta E) to standardowa miara różnicy kolorów
    deltaE = sqrt(sum((lab_reg - lab_temp).^2, 3));
    deltaE(~realMask) = 0;  % zerujemy różnice poza bankotem
    
    % Adaptacyjne progowanie różnic kolorystycznych
    % Używamy średniej + 3 odchylenia standardowego jako progu
    % Dzięki temu próg dostosowuje się do konkretnego obrazu
    if any(deltaE(:) > 0)
        threshold = mean(deltaE(deltaE > 0)) + 3*std(deltaE(deltaE > 0));
    else
        threshold = inf;  % jeśli nie ma różnic, ustawiamy nieosiągalny próg
    end
    colorAnomalies = deltaE > threshold;  % binaryzacja - tworzymy maskę anomalii
    colorAnomalies = bwareaopen(colorAnomalies, MIN_COLOR_ANOMALY);  % usuwamy małe obiekty
    colorAnomalies = imclose(colorAnomalies, strel('disk', 5));  % łączymy bliskie plamy
    
    nColor = bwconncomp(colorAnomalies).NumObjects;
    fprintf('   > Znaleziono %d anomalii kolorów\n', nColor);
    
    % Agregacja wykrytych defektów - łączymy oba typy defektów
    fprintf('   > Łączę wykryte wady...\n');
    defectMap = missingParts | colorAnomalies;  % suma logiczna obu map

    % Ekstrakcja cech i klasyfikacja defektów
    % Analizujemy każdy wykryty defekt osobno
    cc = bwconncomp(defectMap);  % znajdujemy połączone komponenty
    rawStats = regionprops(cc, 'Area', 'BoundingBox', 'Eccentricity', 'PixelIdxList');

    % Klasyfikacja typu (rozdarcie vs plama) i wagi (minor/major)
    % Tworzymy szczegółową strukturę z informacjami o każdym defekcie
    stats = struct('Area', {}, 'BoundingBox', {}, 'Eccentricity', {}, 'Type', {}, 'Severity', {}, 'AreaPct', {});
    for i = 1:numel(rawStats)
        s = rawStats(i);
        areaPct = double(s.Area) / double(areaTemplate);  % jaki % powierzchni banknotu zajmuje defekt

        % Determinacja typu defektu na podstawie składowych
        % Sprawdzamy który typ defektu dominuje w tym obszarze
        pix = s.PixelIdxList;  % indeksy pikseli defektu
        fracMissing = sum(missingParts(pix)) / numel(pix);  % frakcja pikseli z brakujących części
        fracColor = sum(colorAnomalies(pix)) / numel(pix);  % frakcja pikseli z anomalii kolorystycznych
        if fracMissing >= fracColor
            typeStr = 'Rozdarcie / Zagniecenie';
        else
            typeStr = 'Plama / Napis';
        end

        % Ocena ciężkości defektu zgodnie z ustalonymi kryteriami
        % Sprawdzamy czy defekt dotyka brzegu obrazu
        isEdgeTouch = false;
        bb = s.BoundingBox;  % [x, y, width, height]
        x1 = round(max(1, bb(1))); y1 = round(max(1, bb(2)));
        x2 = round(min(size(missingParts,2), bb(1)+bb(3)));
        y2 = round(min(size(missingParts,1), bb(2)+bb(4)));
        if x1 <= EDGE_MARGIN || y1 <= EDGE_MARGIN || ...
           (size(missingParts,2) - x2) <= EDGE_MARGIN || (size(missingParts,1) - y2) <= EDGE_MARGIN
            isEdgeTouch = true;
        end

        % Klasyfikacja jako "major" (poważny) jeśli:
        % 1. Powierzchnia >= 1.5% całego banknotu, LUB
        % 2. Jest to wydłużony defekt przy brzegu (może być rozdarcie)
        majorByArea = areaPct >= 0.015;  % duży obszar
        majorEdgeTear = (s.Eccentricity > 0.92) && isEdgeTouch && (areaPct >= 0.001);  % wydłużone przy brzegu
        if majorByArea || majorEdgeTear
            severityStr = 'major';  % poważny defekt
        else
            severityStr = 'minor';  % drobny defekt
        end

        % Zapisujemy wszystkie informacje o defekcie
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
  
    


   







%% Funkcja displayResults - wizualizuje wyniki analizy
%
% Funkcja realizuje wizualizację wyników analizy z zaznaczeniem wykrytych defektów.
% Pokazuje tylko poważne wady (major) i wydaje werdykt o akceptacji banknotu.
%
% Argumenty wejściowe:
%   image - obraz banknotu do wyświetlenia
%   stats - struktura z informacjami o defektach
%   nominal - rozpoznany nominał
%   eccentricityThreshold - próg mimośrodowości
%   acceptance - struktura z kryteriami akceptacji
function displayResults(image, stats, nominal, eccentricityThreshold, acceptance)
    % Funkcja realizuje wizualizację wyników analizy z zaznaczeniem wykrytych defektów
    fprintf('Krok 5: Wyświetlanie wyników (tylko poważne wady)...\n');

    figure;  % tworzymy nowe okno
    imshow(image);  % wyświetlamy obraz banknotu
    hold on;  % pozwala rysować na obrazie

    % Przypadek braku wykrytych defektów - banknot idealny
    if isempty(stats)
        titleStr = sprintf('Nominał: %s. Status: ZAAKCEPTOWANY', nominal);
        fprintf('\n--- WERDYKT: BANKNOT ZAAKCEPTOWANY (brak wad) ---\n');
        title(titleStr, 'FontSize', 14); hold off; return;
    end

    % Selekcja defektów o statusie poważnym
    % Filtrujemy tylko te defekty które są sklasyfikowane jako "major"
    isMajor = arrayfun(@(s) isfield(s,'Severity') && strcmpi(s.Severity,'major'), stats);
    majorStats = stats(isMajor);
    nMajor = numel(majorStats);

    % Wizualizacja poważnych defektów na obrazie
    for i = 1:length(majorStats)
        bb = majorStats(i).BoundingBox;  % prostokąt okalający defekt
        % Rysujemy czerwony prostokąt wokół defektu
        rectangle('Position', bb, 'EdgeColor', 'r', 'LineWidth', 2);

        % Generacja etykiety tekstowej - typ defektu
        if isfield(majorStats(i), 'Type')
            label = majorStats(i).Type;  % używamy zapisanego typu
        else
            % Fallback - jeśli brak pola Type, klasyfikujemy na podstawie mimośrodowości
            if majorStats(i).Eccentricity > eccentricityThreshold
                label = 'Rozdarcie / Zagniecenie';
            else
                label = 'Plama / Napis';
            end
        end
        label = [label ' (poważna)'];  % dodajemy oznaczenie wagi
        % Wyświetlamy etykietę nad defektem (10px powyżej)
        text(bb(1), bb(2) - 10, label, 'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold');
        % Logujemy informacje o defekcie do konsoli
        fprintf('   > Wada #%d: Typ: %s, Powierzchnia: %d px (%.3f%%), MAJOR\n', ...
            i, label, majorStats(i).Area, 100*majorStats(i).AreaPct);
    end

    % Procedura decyzyjna akceptacji banknotu
    % Jeśli nie ma poważnych wad - akceptujemy banknot
    % Jeśli są poważne wady - wymagana kontrola ręczna
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

%% Funkcja pomocnicza ternary - operator trójargumentowy
% Przydatna do zwięzłych warunków if-else w jednej linii
% Użycie: wynik = ternary(warunek, wartość_jeśli_prawda, wartość_jeśli_fałsz)
function out = ternary(cond, a, b)
    if cond
        out = a;  % zwracamy a jeśli warunek prawdziwy
    else
        out = b;  % zwracamy b jeśli warunek fałszywy
    end
end