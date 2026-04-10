
"""
Self-Organised Criticality in the Forest Fire Model
====================================================
Drossel-Schwabl (1992) forest-fire model with instant fire spread.

States:  EMPTY = 0 | TREE = 1

Algorithm per time step
-----------------------
  1. Each TREE cell is struck by lightning independently with probability f
     -> its connected cluster burns instantly  (avalanche)
  2. Each EMPTY cell grows a tree with probability p  (noise / recovery)

Author : Rayeli Nandy
Date   : 10 April 2026
"""

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from collections import deque
from scipy import stats

EMPTY = 0
TREE  = 1


def init_grid(L, rho=0.6, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    return rng.choice([EMPTY, TREE], size=(L, L),
                      p=[1 - rho, rho]).astype(np.int8)


def bfs_cluster(grid, r0, c0, L):
    cluster, queue = [(r0, c0)], deque([(r0, c0)])
    visited = {(r0, c0)}
    while queue:
        r, c = queue.popleft()
        for dr, dc in ((-1, 0), (1, 0), (0, -1), (0, 1)):
            nr, nc = r + dr, c + dc
            if (0 <= nr < L and 0 <= nc < L
                    and grid[nr, nc] == TREE
                    and (nr, nc) not in visited):
                visited.add((nr, nc))
                cluster.append((nr, nc))
                queue.append((nr, nc))
    return cluster


def step(grid, p, f, L, rng):
    step_avalanches = []
    tree_rows, tree_cols = np.where(grid == TREE)
    if len(tree_rows):
        mask = rng.random(len(tree_rows)) < f
        for r0, c0 in zip(tree_rows[mask], tree_cols[mask]):
            if grid[r0, c0] == TREE:
                cluster = bfs_cluster(grid, r0, c0, L)
                for r, c in cluster:
                    grid[r, c] = EMPTY
                step_avalanches.append(len(cluster))
    empty_rows, empty_cols = np.where(grid == EMPTY)
    if len(empty_rows):
        mask = rng.random(len(empty_rows)) < p
        grid[empty_rows[mask], empty_cols[mask]] = TREE
    return grid, step_avalanches


def run_simulation(L=80, p=0.05, f=0.001, steps=8000, burn_in=2000, seed=42):
    rng  = np.random.default_rng(seed)
    grid = init_grid(L, rho=0.6, rng=rng)
    avalanche_sizes, densities, snapshots = [], [], {}
    snap_at = {burn_in, burn_in + steps // 2, burn_in + steps - 1}
    for t in range(burn_in + steps):
        grid, av_list = step(grid, p, f, L, rng)
        if t >= burn_in:
            densities.append(np.sum(grid == TREE) / L**2)
            avalanche_sizes.extend(av_list)
        if t in snap_at:
            snapshots[t - burn_in] = grid.copy()
    return dict(avalanche_sizes=avalanche_sizes, densities=densities,
                snapshots=snapshots,
                params=dict(L=L, p=p, f=f, steps=steps, burn_in=burn_in))


def plot_avalanche_distribution(sizes, outdir):
    s = np.array([x for x in sizes if x > 0], dtype=float)
    bins = np.logspace(np.log10(max(s.min(), 1)), np.log10(s.max()), 40)
    counts, edges = np.histogram(s, bins=bins)
    centres = (edges[:-1] + edges[1:]) / 2
    mask = counts > 3
    slope, intercept, r, *_ = stats.linregress(
        np.log10(centres[mask]), np.log10(counts[mask]))
    alpha = -slope
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.loglog(centres[mask], counts[mask], "o", ms=5, color="steelblue",
              label="Simulation data")
    ax.loglog(centres[mask], 10**intercept * centres[mask]**slope,
              "--", color="tomato", lw=2,
              label=rf"Power-law fit  $\alpha={alpha:.2f}$  ($R^2={r**2:.3f}$)")
    ax.set_xlabel(r"Avalanche size $s$ (trees burned)", fontsize=12)
    ax.set_ylabel(r"Frequency $N(s)$", fontsize=12)
    ax.set_title("Avalanche Size Distribution", fontsize=13, fontweight="bold")
    ax.legend(fontsize=10)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"avalanche_distribution.{ext}"), dpi=150)
    plt.show()
    plt.close(fig)
    print(f"  alpha = {alpha:.3f}   R2 = {r**2:.4f}")
    return alpha, r**2


def plot_density_timeseries(densities, outdir):
    rho = np.array(densities)
    fig, ax = plt.subplots(figsize=(8, 3.5))
    ax.plot(rho, lw=0.5, color="forestgreen", alpha=0.85, label="Tree density")
    mu = rho.mean()
    ax.axhline(mu, color="darkred", lw=1.8, ls="--",
               label=rf"Mean $\rho={mu:.3f}$")
    ax.fill_between(range(len(rho)), mu - rho.std(), mu + rho.std(),
                    color="green", alpha=0.12, label=r"$\pm1\sigma$")
    ax.set_xlabel("Time step", fontsize=12)
    ax.set_ylabel(r"Tree density $\rho$", fontsize=12)
    ax.set_title("Tree Density at Self-Organised Critical State",
                 fontsize=13, fontweight="bold")
    ax.legend(fontsize=10)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"density_timeseries.{ext}"), dpi=150)
    plt.show()
    plt.close(fig)


def plot_snapshots(snapshots, outdir):
    cmap = mcolors.ListedColormap(["#d4b07a", "#2d6a2d"])
    norm = mcolors.BoundaryNorm([0, 0.5, 1.5], cmap.N)
    keys = sorted(snapshots.keys())
    labels = ["Early (burn-in end)",
              f"Mid  (t={keys[1]:,})",
              f"Late / SOC (t={keys[2]:,})"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    for ax, k, lbl in zip(axes, keys, labels):
        ax.imshow(snapshots[k], cmap=cmap, norm=norm, interpolation="nearest")
        rho = np.mean(snapshots[k] == TREE)
        ax.set_title(f"{lbl}\n$\\rho = {rho:.3f}$", fontsize=11)
        ax.axis("off")
    fig.suptitle("Forest Grid Snapshots  (tan=empty, green=tree)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"grid_snapshots.{ext}"), dpi=150)
    plt.show()
    plt.close(fig)


def plot_connectivity_analysis(outdir, L=60):
    p = 0.05
    f_values = [0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02]
    records = []
    print("  Connectivity sweep ...")
    for f in f_values:
        rng  = np.random.default_rng(7)
        grid = init_grid(L, rho=0.6, rng=rng)
        sizes, rhos = [], []
        for t in range(6000):
            grid, av = step(grid, p, f, L, rng)
            if t >= 1000:
                rhos.append(np.sum(grid == TREE) / L**2)
                sizes.extend(av)
        rho_m = np.mean(rhos)
        s_m   = np.mean(sizes) if sizes else 0
        records.append((f, rho_m, s_m, len(sizes)))
        print(f"    f={f:.4f}  rho={rho_m:.3f}  <s>={s_m:.1f}  n={len(sizes)}")
    rho_v = [r[1] for r in records]
    s_v   = [r[2] for r in records]
    n_v   = [r[3] for r in records]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].plot(rho_v, s_v, "o-", color="steelblue", ms=8, lw=2)
    axes[0].set_xlabel(r"Mean tree density $\rho$  (connectivity)", fontsize=12)
    axes[0].set_ylabel(r"Mean avalanche size $\langle s \rangle$", fontsize=12)
    axes[0].set_title("Connectivity vs. Mean Avalanche Size",
                      fontsize=12, fontweight="bold")
    axes[0].grid(True, alpha=0.3)
    axes[1].plot(rho_v, n_v, "s-", color="tomato", ms=8, lw=2)
    axes[1].set_xlabel(r"Mean tree density $\rho$  (connectivity)", fontsize=12)
    axes[1].set_ylabel("Number of avalanches", fontsize=12)
    axes[1].set_title("Connectivity vs. Avalanche Frequency",
                      fontsize=12, fontweight="bold")
    axes[1].grid(True, alpha=0.3)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"connectivity_analysis.{ext}"), dpi=150)
    plt.show()
    plt.close(fig)


def main():
    outdir = "figures"
    os.makedirs(outdir, exist_ok=True)
    print("=" * 60)
    print("  Forest Fire SOC Simulation")
    print("  Drossel-Schwabl Model (instant spread)")
    print("=" * 60)

    print("\n[1/4] Main simulation (L=80, p=0.05, f=0.001, 8000 steps) ...")
    res = run_simulation(L=80, p=0.05, f=0.001, steps=8000, burn_in=2000, seed=42)
    sizes, densities, snapshots = res["avalanche_sizes"], res["densities"], res["snapshots"]
    print(f"  Avalanche events : {len(sizes)}")
    print(f"  Mean density     : {np.mean(densities):.4f}")
    print(f"  Max avalanche    : {max(sizes)}")

    print("\n[2/4] Avalanche size distribution ...")
    alpha, r2 = plot_avalanche_distribution(sizes, outdir)

    print("\n[3/4] Density time series & snapshots ...")
    plot_density_timeseries(densities, outdir)
    plot_snapshots(snapshots, outdir)

    print("\n[4/4] Connectivity analysis (L=60) ...")
    plot_connectivity_analysis(outdir, L=60)

    print(f"\nAll figures -> '{outdir}/'")
    print("=" * 60)
    print(f"  alpha = {alpha:.3f}  (literature: ~1.1-1.2)")
    print(f"  R2    = {r2:.4f}")
    print("=" * 60)


if __name__ == "__main__":
    main()
