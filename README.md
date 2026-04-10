# Self-Organised Criticality in the Forest Fire Model

**Individual Project — Complexity Science**
**Submission: April 10, 2026**

---

## Overview

This project applies **complexity science** to the **Drossel–Schwabl forest-fire model** (1992), a canonical system exhibiting **self-organised criticality (SOC)**. The forest fire model captures all three hallmarks of SOC:

| Component | In this model |
|-----------|--------------|
| **Noise** | Lightning strikes (random, independent, prob. `f` per tree) |
| **Avalanche** | Entire connected tree cluster burns instantaneously |
| **Connectivity** | Tree density ρ — controls cluster size and cascade potential |

---

## Repository Structure

```
.
├── simulation.py          # Main Python simulation
├── main.tex               # LaTeX journal manuscript (source)
├── main.pdf               # Compiled PDF report
├── requirements.txt       # Python dependencies
├── figures/               # Generated figures (PDF + PNG)
│   ├── avalanche_distribution.png
│   ├── connectivity_analysis.png
│   ├── density_timeseries.png
│   └── grid_snapshots.png
└── README.md
```

---

## How to Run

### 1. Install dependencies
```bash
pip install -r requirements.txt
```

### 2. Run the simulation
```bash
python simulation.py
```

This will:
- Run the main forest-fire simulation (L=80, 8000 steps)
- Generate all 4 figures in `figures/`
- Print the power-law exponent and R²

Expected output:
```
Avalanche events : ~19,000
Mean density     : ~0.39
Max avalanche    : ~2000
Power-law alpha  : ~0.44
```

### 3. Compile the LaTeX report
```bash
pdflatex main.tex
pdflatex main.tex   # run twice for cross-references
```

---

## Model Parameters

| Parameter | Symbol | Value | Description |
|-----------|--------|-------|-------------|
| Grid size | L | 80 | Square lattice L×L |
| Growth prob. | p | 0.05 | Empty → Tree per step |
| Lightning prob. | f | 0.001 | Tree struck per step |
| Drive ratio | θ = f/p | 0.02 | Controls SOC regime |
| Burn-in | — | 2000 | Steps discarded |
| Production | — | 8000 | Steps analysed |

---

## Key Results

- **Avalanche distribution:** Fat-tailed, consistent with power law N(s) ∝ s^{−α}, α ≈ 0.44
- **Critical density:** System self-organises to ρ ≈ 0.39 without external tuning
- **Connectivity:** Mean avalanche size increases monotonically with tree density

---

## Dependencies

- Python ≥ 3.10
- numpy, matplotlib, scipy

---

## Acknowledgements

Simulation design, code, and report assisted by **Claude** (Anthropic, claude.ai).
All scientific results verified independently.
