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

## Wykorzystane techniki

### Identyfikacja nominału
- **ORB (Oriented FAST and Rotated BRIEF)** - detekcja i deskrypcja punktów kluczowych
- **Feature matching** - dopasowanie deskryptorów binarnych (Hamming distance)
- **Transformacja rzutowa** - wyrównanie geometryczne (estimateGeometricTransform2D)

### Segmentacja banknotu
- **Progowanie Otsu** - automatyczna binaryzacja kanału saturacji (HSV)
- **Operacje morfologiczne** - zamknięcie, wypełnianie dziur, usuwanie małych obszarów

### Detekcja uszkodzeń mechanicznych
- **Różnica masek logicznych** - porównanie maski rzeczywistej z wzorcową
- **Analiza komponentów spójnych** - identyfikacja brakujących fragmentów
- **Filtracja brzegowa** - eliminacja artefaktów rejestracji

### Detekcja anomalii wizualnych
- **Przestrzeń kolorów LAB** - analiza różnic perceptualnych (metryka ΔE)
- **Filtr Gaussa** - wygładzenie przed porównaniem kolorów
- **Progowanie adaptacyjne** - μ + 3σ dla różnic kolorystycznych

### Klasyfikacja defektów
- **Analiza regionów** - powierzchnia, mimośród, położenie brzegowe
- **Próg dwupoziomowy** - minor (≤1.5%) vs major (>1.5% lub krawędziowe)

## Parametry konfiguracyjne

Kluczowe parametry znajdują się w pliku `check_banknote.m`:

- `MAX_FEATURES` - maksymalna liczba ekstrahowanych cech (domyślnie 10000)
- `MATCH_THRESHOLD` - próg podobieństwa przy dopasowywaniu (0.8)
- `MIN_DEFECT_AREA` - minimalna powierzchnia wykrywanego defektu (150 px)
- `MAJOR_AREA_PCT` - próg klasyfikacji wady jako poważnej (1.5%)

## Wyniki

Program generuje wizualizację z zaznaczonymi defektami i wydaje werdykt:
- **ZAAKCEPTOWANY** - brak poważnych wad
- **DO SPRAWDZENIA** - wykryto defekty wymagające inspekcji


