# Komputerowe Symulacje Cieczy Nieściśliwych
Kod źródłowy oraz wygenerowane dzięki jego wynikom animacje dla wyników z pracy magisterskiej. Jej celem było stworzenie algorytmu na bazie metody MAC (*Marker-And-Cell*) opierającego się na równaniach Naviera-Stokesa do tworzenia symulacji w czasie cieczy nieściśliwych o różnych parametrach i konfiguracjach.

## 📁 Code
Kod źródłowy w C++. Podzielony na klasy:
- ``System`` - główny algorytm i jego ustawienia,
- ``Type`` - typ wyliczeniowy, określa typy komórek na siatce obliczeniowej,
- ``Cell`` - typ opakowujący ``Type``, komórka dla siatki określającej flagi,
- ``Particle`` - klasa reprezentująca cząstkę znaczoną,
- ``Matrix`` - klasa pomocnicza do reprezentacji macierzy.

## 📁 Plots
Pliki generujące wykresy i animacje. 

## 📁 Animations
Zawiera przypadki rozpatrywane i opisywane w pracy w formie gifów. Niektóre przykłady poniżej.

### Problem *lid-driven cavity* 
![Lid-driven cavity problem animation](https://github.com/TMaczek/msc_thesis/blob/main/animations/6_2_anim.gif)

### Fala załamująca się na płytkiej wodzie
![shallow wave breaking](https://github.com/TMaczek/msc_thesis/blob/main/animations/7_2_anim.gif)

### Opuszczona tama (*broken dam problem*) z przeszkodą 
![broken dam with obstacle](https://github.com/TMaczek/msc_thesis/blob/main/animations/7_4_anim.gif)
