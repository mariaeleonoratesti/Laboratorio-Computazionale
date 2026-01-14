# XY Model Monte Carlo Simulation (2D)

Questo repository contiene un’implementazione completa del **modello XY bidimensionale** basata su simulazioni **Monte Carlo con algoritmo di Metropolis**, scritta in **Fortran 90** e interfacciata a **Python** tramite `f2py`.  
L’analisi dei risultati, la visualizzazione delle osservabili e lo studio della transizione di **Berezinskii–Kosterlitz–Thouless (BKT)** sono effettuati tramite un **Jupyter Notebook**.

---

## Contenuto del repository

- `xy_fortran.f90`  
  Modulo Fortran 90 che implementa:
  - il modello XY 2D con condizioni periodiche al contorno,
  - l’algoritmo di Metropolis,
  - il calcolo di energia, magnetizzazione e correlazioni,
  - l’identificazione dei vortici,
  - scansioni in temperatura,
  - analisi statistiche (binning, bootstrap, autocorrelazione),
  - stima della temperatura critica BKT.

- `xy_organized.html`
  html del Jupyter Notebook per:
  - compilare e importare il modulo Fortran in Python,
  - eseguire simulazioni Monte Carlo,
  - analizzare le osservabili fisiche,
  - produrre grafici di energia, magnetizzazione, suscettibilità,
  - studiare la densità di vortici e la transizione BKT.

---

## Modello fisico

Il **modello XY bidimensionale** descrive spin classici di modulo unitario, interagenti tramite l’Hamiltoniana H = - Σ_{<i,j>} cos(θ_i - θ_j).

In due dimensioni il modello non presenta ordine magnetico a lungo raggio a temperatura finita (teorema di Mermin–Wagner), ma mostra una **transizione topologica di Berezinskii–Kosterlitz–Thouless**, associata alla dissociazione di coppie vortice–antivortice.

---

## Requisiti

### Fortran
- Compilatore Fortran compatibile con Fortran 90 (es. `gfortran`)

### Python
- Python ≥ 3.8
- Python < 3.12 (per compatibilità con f2py)
- `numpy`
- `matplotlib`
- `jupyter`
- `scipy`

---

## Compilazione del modulo Fortran

Il modulo Fortran può essere compilato come estensione Python usando `f2py` come fatto nelle prime celle del Jupyter Notebook:

```bash
python -m numpy.f2py -c xy_fortran.f90 -m xy_model --opt='-O3 -march=native'
