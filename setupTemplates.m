clear; clc; close all;

fprintf('Rozpoczynam analizę wzorców i filtrowanie cech...\n');

% Konfiguracja parametrów
templateFolder = 'templates';
templateFiles = {
    '10zl_wzorcowe_1.jpg', '10zl_wzorcowe_2.jpg', '10zl_wzorcowe_3.jpg', ... % banknoty 10 PLN
    '20zl_wzorcowe_1.jpg', '20zl_wzorcowe_2.jpg', '20zl_wzorcowe_3.jpg', ... % banknoty 20 PLN
    '50zl_wzorcowe_1.jpg', '50zl_wzorcowe_2.jpg', '50zl_wzorcowe_3.jpg', ... % banknoty 50 PLN
    '100zl_wzorcowe_1.jpg', '100zl_wzorcowe_2.jpg', '100zl_wzorcowe_3.jpg', ... % banknoty 100 PLN
    '200zl_wzorcowe_1.jpg', '200zl_wzorcowe_2.jpg', '200zl_wzorcowe_3.jpg' ... % banknoty 200 PLN
};
templateNames = {'10 PLN #1', '10 PLN #2', '10 PLN #3', ...
                 '20 PLN #1', '20 PLN #2', '20 PLN #3', ...
                 '50 PLN #1', '50 PLN #2', '50 PLN #3', ...
                 '100 PLN #1', '100 PLN #2', '100 PLN #3', ...
                 '200 PLN #1', '200 PLN #2', '200 PLN #3'};

templateData = struct('Name', {}, 'Image', {}, 'Points', {}, 'Features', {}, 'Mask', {});
MAX_FEATURES = 10000; % maksymalna liczba ekstrahowanych cech

% Przetwarzanie obrazów wzorcowych
fprintf('Liczba plików do przetworzenia: %d\n', length(templateFiles));
for i = 1:length(templateFiles)
    fileName = templateFiles{i};
    fullPath = fullfile(templateFolder, fileName);
    
    fprintf('Sprawdzam plik: %s\n', fullPath);
    if ~exist(fullPath, 'file')
        warning('Nie znaleziono pliku wzorca: %s. Pomijam.', fullPath);
        continue;
    end
    
    fprintf('Przetwarzam: %s...\n', fileName);
    
    imgColor = imread(fullPath);
    imgGray = rgb2gray(imgColor);
    
    % Tworzenie maski w celu separacji banknotu od tła
    fprintf('   > Tworzę maskę...\n');
    imgHSV = rgb2hsv(imgColor);
    sChannel = imgHSV(:,:,2); 
    level = graythresh(sChannel);
    mask = sChannel > level;
    mask = imfill(mask, 'holes');
    mask = bwareaopen(mask, 5000); 

    % Detekcja punktów charakterystycznych
    points_all = detectORBFeatures(imgGray);
    
    % Filtracja punktów znajdujących się wewnątrz maski banknotu
    fprintf('   > Filtruję cechy...\n');
    locations = points_all.Location;
    indices = sub2ind(size(mask), round(locations(:,2)), round(locations(:,1)));
    is_inside = mask(indices);
    points_inside = points_all(is_inside);

    % Selekcja punktów o najwyższej sile odpowiedzi
    points_filtered_strong = selectStrongest(points_inside, MAX_FEATURES);

    
    fprintf('   > Znaleziono %d cech na banknocie, wybrano %d najsilniejszych.\n', ...
        points_inside.Count, points_filtered_strong.Count);

    % Ekstrakcja deskryptorów dla wybranych punktów
    [features_filtered, valid_points_filtered] = extractFeatures(imgGray, points_filtered_strong);
    
    % Zapis danych do struktury wynikowej
    templateData(i).Name = templateNames{i};
    templateData(i).Image = imgColor;
    templateData(i).Points = valid_points_filtered;
    templateData(i).Features = features_filtered; 
    templateData(i).Mask = mask;
end

% Zapis przetworzonych danych wzorcowych
if ~isempty(templateData)
    save('templateFeatures.mat', 'templateData');
    fprintf('\nPrzetworzono %d wzorców i zapisano cechy do "templateFeatures.mat".\n', length(templateData));
else
    fprintf('\nBŁĄD: Nie przetworzono żadnych wzorców.\n');
end