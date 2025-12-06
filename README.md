# System wizyjnej kontroli banknotów

System służy do automatycznej identyfikacji nominału polskich banknotów oraz detekcji ich uszkodzeń mechanicznych i wizualnych.

## Wymagania

- MATLAB R2020a lub nowszy
- Computer Vision Toolbox
- Image Processing Toolbox

## Struktura projektu

```
├── setupTemplates.m          # Przygotowanie bazy wzorców
├── check_banknote.m          # Główny program weryfikacji
├── templates/                # Katalog z obrazami wzorcowymi
└── test/                     # Katalog z obrazami testowymi
```

## Instrukcja użytkowania

### 1. Przygotowanie bazy wzorców

Przed pierwszym uruchomieniem należy wygenerować bazę cech charakterystycznych:

```matlab
run setupTemplates.m
```

Skrypt przetworzy obrazy wzorcowe banknotów (10, 20, 50, 100, 200 PLN) i zapisze wyekstrahowane cechy do pliku `templateFeatures.mat`.

### 2. Weryfikacja banknotu

Po przygotowaniu bazy wzorców można przystąpić do analizy:

```matlab
run check_banknote.m
```

Program poprosi o podanie ścieżki do obrazu testowego, a następnie przeprowadzi:
- Identyfikację nominału
- Rejestrację geometryczną obrazu
- Detekcję defektów i uszkodzeń

## Algorytm działania

### Etap przygotowania wzorców

1. Wczytanie obrazów referencyjnych
2. Generacja maski w celu separacji banknotu od tła (przestrzeń HSV)
3. Detekcja punktów kluczowych (ORB)
4. Filtracja punktów do obszaru banknotu
5. Ekstrakcja deskryptorów

### Etap weryfikacji

1. **Identyfikacja nominału** - dopasowanie cech testowych do wzorców
2. **Rejestracja geometryczna** - wyrównanie obrazu metodą transformacji rzutowej
3. **Transformacja maski wzorcowej** - dopasowanie do geometrii testowej
4. **Detekcja defektów**:
   - Analiza brakujących fragmentów (zagięcia, rozdarcia)
   - Detekcja anomalii kolorystycznych (plamy, napisy)
5. **Klasyfikacja wad** - ocena typu i ciężkości defektu

## Parametry konfiguracyjne

Kluczowe parametry znajdują się w pliku `check_banknote.m`:

- `MAX_FEATURES` - maksymalna liczba ekstrahowanych cech (domyślnie 10000)
- `MATCH_THRESHOLD` - próg podobieństwa przy dopasowywaniu (0.8)
- `MIN_DEFECT_AREA` - minimalna powierzchnia wykrywanego defektu (150 px)
- `MAJOR_AREA_PCT` - próg klasyfikacji wady jako poważnej (1.5%)

## Uwagi implementacyjne

- System wykorzystuje punkty charakterystyczne ORB ze względu na odporność na rotację i skalowanie
- Detekcja defektów opiera się na porównaniu maski rzeczywistej z idealną maską wzorca
- Anomalie kolorystyczne wykrywane są w przestrzeni LAB (metryka ΔE)
- Klasyfikacja defektów uwzględnia zarówno powierzchnię, jak i położenie względem brzegów

## Struktura wyników

Program wyświetla obraz z zaznaczonymi defektami oraz prezentuje werdykt:
- **ZAAKCEPTOWANY** - brak poważnych wad
- **DO SPRAWDZENIA** - wykryto defekty wymagające inspekcji


