% --- setupTemplates.m (Wersja 3.1 - Ograniczenie Cech) ---
%
% Cel: Ograniczenie liczby cech do 10 000 najsilniejszych.
% Uruchom: KONIECZNIE PONOWNIE (tylko raz)

clear; clc; close all;

fprintf('Rozpoczynam analizę wzorców i filtrowanie cech...\n');

% --- Konfiguracja ---
templateFolder = 'templates';
templateFiles = {
    '10zl_wzorcowe.jpg', % Użyj tutaj nazwy Twojego pliku wzorcowego
    '20zl_wzorcowe.jpg', 
    '50zl_wzorcowe.jpg', 
    '100zl_wzorcowe.jpg',
    '200zl_wzorcowe.jpg'
};
templateNames = {'10 PLN', '20 PLN', '50 PLN', '100 PLN', '200 PLN'};

templateData = struct('Name', {}, 'Image', {}, 'Points', {}, 'Features', {}, 'Mask', {});
MAX_FEATURES = 10000; % Ustalamy limit 10 000 cech

% --- Przetwarzanie ---
for i = 1:length(templateFiles)
    fileName = templateFiles{i};
    fullPath = fullfile(templateFolder, fileName);
    
    if ~exist(fullPath, 'file')
        warning('Nie znaleziono pliku wzorca: %s. Pomijam.', fullPath);
        continue;
    end
    
    fprintf('Przetwarzam: %s...\n', fileName);
    
    imgColor = imread(fullPath);
    imgGray = rgb2gray(imgColor);
    
    % === Krok 1: Tworzenie MASKI ===
    fprintf('   > Tworzę maskę...\n');
    imgHSV = rgb2hsv(imgColor);
    sChannel = imgHSV(:,:,2); 
    level = graythresh(sChannel);
    mask = sChannel > level;
    mask = imfill(mask, 'holes');
    mask = bwareaopen(mask, 5000); 

    % === Krok 2: Wykrywanie cech ===
    points_all = detectORBFeatures(imgGray);
    
    % === Krok 3: Filtrowanie cech (wewnątrz maski) ===
    fprintf('   > Filtruję cechy...\n');
    locations = points_all.Location;
    indices = sub2ind(size(mask), round(locations(:,2)), round(locations(:,1)));
    is_inside = mask(indices);
    points_inside = points_all(is_inside);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%% TU JEST KLUCZOWA ZMIANA (WERSJA 3.1) %%%
    %
    % Wybierz 10 000 NAJMOCNIEJSZYCH cech spośród tych w masce
    %
    points_filtered_strong = selectStrongest(points_inside, MAX_FEATURES);
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('   > Znaleziono %d cech na banknocie, wybrano %d najsilniejszych.\n', ...
        points_inside.Count, points_filtered_strong.Count);

    % === Krok 4: Ekstrakcja deskryptorów ===
    [features_filtered, valid_points_filtered] = extractFeatures(imgGray, points_filtered_strong);
    
    % === Krok 5: Zapis do struktury ===
    templateData(i).Name = templateNames{i};
    templateData(i).Image = imgColor;
    templateData(i).Points = valid_points_filtered;
    templateData(i).Features = features_filtered; 
    templateData(i).Mask = mask;
end

% --- Zapis do pliku ---
if ~isempty(templateData)
    save('templateFeatures.mat', 'templateData');
    fprintf('\nGotowe! Przetworzono %d wzorców i zapisano "czyste" cechy do "templateFeatures.mat".\n', length(templateData));
else
    fprintf('\nBŁĄD: Nie przetworzono żadnych wzorców.\n');
end