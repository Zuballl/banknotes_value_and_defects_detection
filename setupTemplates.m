% setupTemplates.m
% Skrypt przygotowujący bazę wzorców banknotów do dalszej analizy.
% Dla każdego obrazu wzorcowego wykonuje detekcję cech ORB, filtruje je
% i zapisuje do pliku .mat, który później wykorzystuje check_banknote.m

clear; clc; close all;

fprintf('Rozpoczynam analizę wzorców i filtrowanie cech...\n');

% Konfiguracja parametrów
templateFolder = 'templates';  % folder z obrazami wzorcowymi  % folder z obrazami wzorcowymi

% Lista wszystkich plików wzorcowych - po 3 obrazy dla każdego nominału.
% Mamy tu banknoty: 10, 20, 50, 100 i 200 PLN.
templateFiles = {
    '10zl_wzorcowe_1.jpg', '10zl_wzorcowe_2.jpg', '10zl_wzorcowe_3.jpg', ... % banknoty 10 PLN
    '20zl_wzorcowe_1.jpg', '20zl_wzorcowe_2.jpg', '20zl_wzorcowe_3.jpg', ... % banknoty 20 PLN
    '50zl_wzorcowe_1.jpg', '50zl_wzorcowe_2.jpg', '50zl_wzorcowe_3.jpg', ... % banknoty 50 PLN
    '100zl_wzorcowe_1.jpg', '100zl_wzorcowe_2.jpg', '100zl_wzorcowe_3.jpg', ... % banknoty 100 PLN
    '200zl_wzorcowe_1.jpg', '200zl_wzorcowe_2.jpg', '200zl_wzorcowe_3.jpg' ... % banknoty 200 PLN
};
% Nazwy opisowe dla każdego wzorca - ułatwia debugowanie i raportowanie wyników
templateNames = {'10 PLN #1', '10 PLN #2', '10 PLN #3', ...
                 '20 PLN #1', '20 PLN #2', '20 PLN #3', ...
                 '50 PLN #1', '50 PLN #2', '50 PLN #3', ...
                 '100 PLN #1', '100 PLN #2', '100 PLN #3', ...
                 '200 PLN #1', '200 PLN #2', '200 PLN #3'};

% Struktura przechowująca wszystkie dane o wzorcach:
% - Name: nazwa wzorca (np. "10 PLN #1")
% - Image: oryginalny obraz w kolorze
% - Points: wykryte punkty charakterystyczne
% - Features: deskryptory tych punktów (wektory opisujące cechy)
% - Mask: maska logiczna oddzielająca banknot od tła
templateData = struct('Name', {}, 'Image', {}, 'Points', {}, 'Features', {}, 'Mask', {});

% Limit cech do wyekstrahowania - więcej cech = lepsza dokładność, ale wolniejsze działanie
MAX_FEATURES = 10000;

% Przetwarzanie obrazów wzorcowych - główna pętla
fprintf('Liczba plików do przetworzenia: %d\n', length(templateFiles));
for i = 1:length(templateFiles)
    fileName = templateFiles{i};
    fullPath = fullfile(templateFolder, fileName);
    
    % Sprawdzamy czy plik faktycznie istnieje - jeśli nie, pomijamy go
    fprintf('Sprawdzam plik: %s\n', fullPath);
    if ~exist(fullPath, 'file')
        warning('Nie znaleziono pliku wzorca: %s. Pomijam.', fullPath);
        continue;
    end
    
    fprintf('Przetwarzam: %s...\n', fileName);
    
    % Wczytujemy obraz w kolorze i konwertujemy do skali szarości
    % (detekcja cech ORB działa na obrazach w odcieniach szarości)
    imgColor = imread(fullPath);
    imgGray = rgb2gray(imgColor);
    
    % Tworzenie maski w celu separacji banknotu od tła
    % Wykorzystujemy przestrzeń HSV, bo kanał saturacji (nasycenie)
    % dobrze oddziela kolorowy banknot od jednolitego tła
    fprintf('   > Tworzę maskę...\n');
    imgHSV = rgb2hsv(imgColor);
    sChannel = imgHSV(:,:,2);  % wyciągamy kanał saturacji (S)
    level = graythresh(sChannel);  % automatyczne obliczenie progu metodą Otsu
    mask = sChannel > level;  % binaryzacja - tworzymy maskę logiczną
    mask = imfill(mask, 'holes');  % wypełniamy dziury w masce
    mask = bwareaopen(mask, 5000);  % usuwamy małe obiekty (szum) 

    % Detekcja punktów charakterystycznych algorytmem ORB
    % ORB (Oriented FAST and Rotated BRIEF) znajduje charakterystyczne miejsca na obrazie,
    % które są łatwe do rozpoznania (np. narożniki, krawędzie)
    points_all = detectORBFeatures(imgGray);
    
    % Filtracja punktów znajdujących się wewnątrz maski banknotu
    % Chcemy tylko cechy z banknotu, a nie z tła
    fprintf('   > Filtruję cechy...\n');
    locations = points_all.Location;  % pozycje wszystkich punktów
    % Konwertujemy współrzędne na indeksy w macierzy
    indices = sub2ind(size(mask), round(locations(:,2)), round(locations(:,1)));
    is_inside = mask(indices);  % sprawdzamy które punkty leżą w masce
    points_inside = points_all(is_inside);  % wybieramy tylko te wewnątrz

    % Selekcja punktów o najwyższej sile odpowiedzi
    % Jeśli mamy za dużo punktów, wybieramy tylko te "najlepsze" -
    % czyli te z najsilniejszą odpowiedzią detektora
    points_filtered_strong = selectStrongest(points_inside, MAX_FEATURES);

    
    fprintf('   > Znaleziono %d cech na banknocie, wybrano %d najsilniejszych.\n', ...
        points_inside.Count, points_filtered_strong.Count);

    % Ekstrakcja deskryptorów dla wybranych punktów
    % Deskryptor to wektor liczb opisujący otoczenie punktu - dzięki niemu
    % można później rozpoznać ten sam punkt na innym obrazie
    [features_filtered, valid_points_filtered] = extractFeatures(imgGray, points_filtered_strong);
    
    % Zapis danych do struktury wynikowej
    % Każdy element struktury zawiera wszystkie informacje o jednym wzorcu
    templateData(i).Name = templateNames{i};
    templateData(i).Image = imgColor;
    templateData(i).Points = valid_points_filtered;
    templateData(i).Features = features_filtered; 
    templateData(i).Mask = mask;
end

% Zapis przetworzonych danych wzorcowych do pliku MAT
% Ten plik będzie później wczytywany przez check_banknote.m
if ~isempty(templateData)
    save('templateFeatures.mat', 'templateData');
    fprintf('\nPrzetworzono %d wzorców i zapisano cechy do "templateFeatures.mat".\n', length(templateData));
else
    fprintf('\nBŁĄD: Nie przetworzono żadnych wzorców.\n');
end