

import os
import numpy as np
import matplotlib.pyplot as plt

# Data exports live in ../data relative to this script, so the script works
# no matter where it is launched from.
DATA_DIR = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
)


def load_stress_csv(path):
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()
    rows = []
    for line in lines:
        parts = line.split()
        if len(parts) >= 3:
            try:
                rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
            except ValueError:
                continue
    arr = np.array(rows)
    if arr.size == 0:
        raise ValueError(f"No numeric data rows found in {path}")
    return arr[:, 0], arr[:, 1], arr[:, 2]   # X, Peel, Shear


def goland_reissner(E, nu, Ea, nua, t, ta, Lo, b, P_total, n=400):
    c = Lo / 2
    Pbar = P_total / b
    Ga = Ea / (2 * (1 + nua))

    u2 = np.sqrt(3 * (1 - nu**2) / 2) * (1 / t) * np.sqrt(Pbar / (t * E))
    k = np.cosh(u2 * c) / (np.cosh(u2 * c) + 2 * np.sqrt(2) * np.sinh(u2 * c))

    beta = np.sqrt(8 * Ga * t / (E * ta))

    gamma = (6 * Ea * t / (E * ta)) ** 0.25
    lam = gamma * (c / t)
    Delta = 0.5 * (np.sin(2 * lam) + np.sinh(2 * lam))
    R1 = np.cosh(lam) * np.sin(lam) + np.sinh(lam) * np.cos(lam)
    R2 = -np.cosh(lam) * np.sin(lam) + np.sinh(lam) * np.cos(lam)
    kp = (k * c / t) * np.sqrt(3 * (1 - nu**2) * Pbar / (t * E))

    x = np.linspace(-c, c, n)

    tau = -(Pbar / (8 * c)) * (
        (beta * c / t) * (1 + 3 * k) * np.cosh(beta * x / t) / np.sinh(beta * c / t)
        + 3 * (1 - k)
    )

    A = ((R2 * lam**2 * k) / 2 + lam * kp) / (np.cosh(lam) * np.cos(lam)) \
        * np.cosh(lam * x / c) * np.cos(lam * x / c)
    B = ((R1 * lam**2 * k) / 2 + lam * kp) / (np.sinh(lam) * np.sin(lam)) \
        * np.sinh(lam * x / c) * np.sin(lam * x / c)
    sigma = (Pbar * t) / (c**2 * Delta) * (A + B)

    x_shifted = x + c   # 0 = overlap start, matches the FEA path coordinate

    return x_shifted, -tau, -sigma



BASE = dict(E=70000, nu=0.33, Ea=2000, nua=0.35, t=1.6, ta=0.2, Lo=12.5, b=25, P_total=2000)

SCENARIOS = [
    ("Baseline (ta=0.2, Lo=12.5, Ea=2000)", "base_1.csv",             {},               None),
    ("ta = 0.1 mm",                          "base_1_t01.csv",         {"ta": 0.1},      "ta"),
    ("ta = 0.4 mm",                          "base_1_t04.csv",         {"ta": 0.4},      "ta"),
    ("Lo = 6.25 mm",                         "base_1_overlap625.csv",  {"Lo": 6.25},     "Lo"),
    ("Lo = 25 mm",                           "base_1_overlap250.csv",  {"Lo": 25.0},     "Lo"),
    ("Lo = 50 mm",                           "layover_50.csv",         {"Lo": 50.0},     "Lo"),
    ("Ea = 500 MPa",                         "base_1_E500.csv",        {"Ea": 500},      "Ea"),
    ("Ea = 5000 MPa",                        "base_1_E5000.csv",       {"Ea": 5000},     "Ea"),
]

GROUP_TITLE = {"ta": "Adhesive thickness sweep", "Lo": "Overlap length sweep", "Ea": "Adhesive modulus sweep"}


results = []
for label, fname, overrides, group in SCENARIOS:
    path = os.path.join(DATA_DIR, fname)
    fea_x, fea_peel, fea_shear = load_stress_csv(path)

    params = {**BASE, **overrides}
    an_x, an_tau, an_sigma = goland_reissner(**params)

    an_shear_at_fea = np.interp(fea_x, an_x, an_tau)
    an_peel_at_fea = np.interp(fea_x, an_x, an_sigma)

    shear_rmse = np.sqrt(np.mean((fea_shear - an_shear_at_fea) ** 2))
    peel_rmse = np.sqrt(np.mean((fea_peel - an_peel_at_fea) ** 2))

    fea_shear_peak = fea_shear[np.argmax(np.abs(fea_shear))]
    an_shear_peak = an_tau[np.argmax(np.abs(an_tau))]
    fea_peel_peak = fea_peel[np.argmax(np.abs(fea_peel))]
    an_peel_peak = an_sigma[np.argmax(np.abs(an_sigma))]

    results.append(dict(
        label=label, group=group, params=params,
        fea_x=fea_x, fea_peel=fea_peel, fea_shear=fea_shear,
        an_x=an_x, an_tau=an_tau, an_sigma=an_sigma,
        shear_rmse=shear_rmse, peel_rmse=peel_rmse,
        fea_shear_peak=fea_shear_peak, an_shear_peak=an_shear_peak,
        fea_peel_peak=fea_peel_peak, an_peel_peak=an_peel_peak,
    ))

by_label = {r["label"]: r for r in results}
baseline = results[0]


print(f"{'Scenario':<38}{'Shear pk FEA':>13}{'Shear pk GR':>13}{'%diff':>8}"
      f"{'Peel pk FEA':>13}{'Peel pk GR':>13}{'%diff':>8}{'RMSE(sh)':>10}{'RMSE(pl)':>10}")
print("-" * 128)
for r in results:
    sh_pct = 100 * (r["fea_shear_peak"] - r["an_shear_peak"]) / r["an_shear_peak"]
    pl_pct = 100 * (r["fea_peel_peak"] - r["an_peel_peak"]) / r["an_peel_peak"]
    print(f"{r['label']:<38}{r['fea_shear_peak']:>13.3f}{r['an_shear_peak']:>13.3f}{sh_pct:>7.1f}%"
          f"{r['fea_peel_peak']:>13.3f}{r['an_peel_peak']:>13.3f}{pl_pct:>7.1f}%"
          f"{r['shear_rmse']:>10.3f}{r['peel_rmse']:>10.3f}")

print("\nNote on the peel comparison: Goland & Reissner's original (1944) peel-stress\n"
      "derivation has documented internal inconsistencies (a plane-stress/plane-strain\n"
      "mismatch later corrected by other authors, e.g. Tsai et al.). For this joint's\n"
      "thin/stiff adhesive combination, the classical formula predicts a much smaller\n"
      "edge peel value than FEA shows — a known, citable limitation of the classical\n"
      "theory, not a modeling error. Discuss this explicitly rather than treating it\n"
      "as a validation failure.")


colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]

for group in ["ta", "Lo", "Ea"]:
    members = [baseline] + [r for r in results if r["group"] == group]
    fig, (ax_sh, ax_pl) = plt.subplots(1, 2, figsize=(12, 5))
    for i, r in enumerate(members):
        c = colors[i % len(colors)]
        ax_sh.plot(r["an_x"], r["an_tau"], "-", color=c, label=f'{r["label"]} (GR)')
        ax_sh.plot(r["fea_x"], r["fea_shear"], "--", color=c, label=f'{r["label"]} (FEA)')
        ax_pl.plot(r["an_x"], r["an_sigma"], "-", color=c, label=f'{r["label"]} (GR)')
        ax_pl.plot(r["fea_x"], r["fea_peel"], "--", color=c, label=f'{r["label"]} (FEA)')

    ax_sh.set_title(f"{GROUP_TITLE[group]} — shear stress")
    ax_sh.set_xlabel("Distance along overlap (mm)")
    ax_sh.set_ylabel("Shear stress (MPa)")
    ax_sh.legend(fontsize=8)

    ax_pl.set_title(f"{GROUP_TITLE[group]} — peel stress")
    ax_pl.set_xlabel("Distance along overlap (mm)")
    ax_pl.set_ylabel("Peel stress (MPa)")
    ax_pl.legend(fontsize=8)

    plt.tight_layout()
    outname = f"SLJ_{group}_sweep_comparison.png"
    plt.savefig(outname, dpi=200)
    print(f"Saved {outname}")

plt.show()


MESH_RUNS = [
    ("2 elements", "change_mesh_mesh2.csv"),
    ("4 elements (baseline)", "base_1.csv"),
    ("8 elements", "change_mesh_mesh8.csv"),
]

print("\nMesh convergence (baseline geometry, varying through-thickness element count):")
print(f"{'Mesh':<24}{'Peak shear':>12}{'Peak peel':>12}{'%chg shear':>12}{'%chg peel':>12}")
print("-" * 72)
mesh_results = []
for label, fname in MESH_RUNS:
    x, peel, shear = load_stress_csv(os.path.join(DATA_DIR, fname))
    sh_peak = shear[np.argmax(np.abs(shear))]
    pl_peak = peel[np.argmax(np.abs(peel))]
    mesh_results.append((label, sh_peak, pl_peak))

ref_sh, ref_pl = mesh_results[0][1], mesh_results[0][2]  # 2-element run as reference
for i, (label, sh_peak, pl_peak) in enumerate(mesh_results):
    if i == 0:
        print(f"{label:<24}{sh_peak:>12.3f}{pl_peak:>12.3f}{'--':>12}{'--':>12}")
    else:
        prev_sh, prev_pl = mesh_results[i-1][1], mesh_results[i-1][2]
        sh_chg = 100 * (sh_peak - prev_sh) / prev_sh
        pl_chg = 100 * (pl_peak - prev_pl) / prev_pl
        print(f"{label:<24}{sh_peak:>12.3f}{pl_peak:>12.3f}{sh_chg:>11.2f}%{pl_chg:>11.2f}%")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
labels = [m[0] for m in mesh_results]
sh_vals = [m[1] for m in mesh_results]
pl_vals = [m[2] for m in mesh_results]
ax1.plot(labels, sh_vals, "o-", color="#d62728")
ax1.set_ylabel("Peak shear stress (MPa)")
ax1.set_title("Mesh convergence — shear")
ax1.tick_params(axis='x', rotation=20)
ax2.plot(labels, pl_vals, "o-", color="#1f77b4")
ax2.set_ylabel("Peak peel stress (MPa)")
ax2.set_title("Mesh convergence — peel")
ax2.tick_params(axis='x', rotation=20)
plt.tight_layout()
plt.savefig("SLJ_mesh_convergence.png", dpi=200)
print("Saved SLJ_mesh_convergence.png")
plt.show()
