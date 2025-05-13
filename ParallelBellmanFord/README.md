# Bellman-Ford Parallel Algorithm

Projekat implementira sekvencijalnu i tri optimizovane paralelne verzije Bellman-Ford algoritma koristeći C++ i OpenMP. Fokus je na efikasnosti, skalabilnosti i tačnosti algoritama za pronalaženje najkraćih puteva u grafovima sa mogućim negativnim težinama grana.

## Sadržaj

- Sekvencijalna Bellman-Ford implementacija
- Paralelna verzija 1 (edge-based sa `#pragma omp`)
- Paralelna verzija 2 (vertex-based sa lokalnim kopijama)
- Paralelna verzija 3 (hibridna sa optimizovanim chunking-om)
- Automatsko testiranje na više veličina grafova
- Export rezultata u CSV format (`bellman_ford_results.csv`)
- Keširanje generisanih grafova sa LRU policy

## Kako pokrenuti

### 1. Instaliraj kompajler sa OpenMP podrškom

Na Linux-u:

```bash
sudo apt update
sudo apt install g++ build-essential
```

Na macOS:

1. Instalirati Xcode command line tools:
   ```bash
   xcode-select --install
   ```
2. (Opciono) Instalirati GCC sa OpenMP podrškom:
   ```bash
   brew install gcc
   ```
   > **Napomena:** Da bi se koristio OpenMP možda je potrebno kompajlirati sa `g++-13`:
   ```bash
   g++-13 -fopenmp -O3 -std=c++17 bellman_ford.cpp -o bellman_ford
   ```
Na Windows-u:

**Opcija 1: MSYS2**  
1. Preuzmi i instaliraj MSYS2 sa [https://www.msys2.org](https://www.msys2.org)  
2. U `MSYS2 MinGW` terminalu, pokreni:
   ```bash
   pacman -Syu
   pacman -S mingw-w64-x86_64-gcc
   ```
3. Kompajliraj:
   ```bash
   g++ -fopenmp -O3 -std=c++17 bellman_ford.cpp -o bellman_ford.exe
   ```

**Opcija 2: Visual Studio**  
- Instaliraj Visual Studio sa komponentom „Desktop development with C++”  
- Kreiraj novi C++ projekat i uključi OpenMP:  
  `Properties > C/C++ > Language > OpenMP Support > Yes`

### 2. Kompajliranje programa

```bash
g++ -fopenmp -O3 -std=c++17 bellman_ford.cpp -o bellman_ford
```

### 3. Pokretanje programa

```bash
./bellman_ford
```

Rezultati će biti sačuvani u fajlu `bellman_ford_results.csv`.

## Struktura projekta

```text
├── ParallelBellmanFord.cpp          # Glavni izvorni kod
├── bellman_ford_results.csv  # Rezultati izvođenja (generiše se)
└── README.md                 # Ovaj fajl
└── Implementacija-ParalelniBellmanFordAlgoritam.pptx #prezentacija
└── Izvještaj_projekat.docx #izvještaj
```

## Primjer rezultata

| Algorithm   | GraphSize | ExecutionTime (s) | Correctness |
|-------------|-----------|-------------------|-------------|
| Sequential  | 500       | 0.012             | true        |
| ParallelV1  | 500       | 0.005             | true        |
| ParallelV2  | 500       | 0.004             | true        |
| ParallelV3  | 500       | 0.003             | true        |

## Tehnologije korištene

- C++
- OpenMP (za paralelizaciju)
- STL (`vector`, `unordered_map`, `chrono`, `fstream`, itd.)

## Napomena

- Program detektuje negativne cikluse u grafovima.
- Koristi se cache za brže testiranje istih veličina grafova (do 5 različitih veličina istovremeno).
- Rezultati se porede radi provjere ispravnosti između paralelnih i sekvencijalnih verzija.


