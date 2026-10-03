# 02 · Composite Shell Buckling

Buckling and post-buckling analysis of a **laminated composite plate/shell**
in Abaqus (eigenvalue buckling + Riks arc-length), cross-checked against a
from-scratch **classical laminated-plate (CLT) buckling calculator**.

---

## Physics / method

* **Abaqus eigenvalue (Buckle) analyses** give the critical load and the first
  six buckling mode shapes for each laminate.
* **Abaqus Riks (arc-length) analysis** traces the non-linear
  load–displacement equilibrium path past the limit point, revealing
  post-buckling behaviour.
* **`classical_buckling_calculator.py`** computes the closed-form critical
  in-plane load `Nx,cr` for a simply-supported specially-orthotropic laminate
  (Whitney / Jones formula), built on a full Classical Lamination Theory
  **A–B–D** matrix implementation (reduced stiffness → `Qbar` → A, B, D).

Three parametric studies are performed by the calculator:

| Study | Variable | Stacking sequence |
|---|---|---|
| **A** | Lay-up angle θ | `[±θ]₂ₛ`, 8 plies, t = 1 mm |
| **B** | Number of plies | `[±45]ₙ`, t_total = 1 mm |
| **C** | Thickness (b/t) | quasi-isotropic `[0/45/−45/90]ₛ` |

---

## Folder layout

```
02-composite-shell-buckling/
├── models/                    composite_shell_buckling.cae / .jnl
├── scripts/                   classical_buckling_calculator.py
├── reports/                   Composite_Shell_Buckling_Project_Report.docx
├── draft/                     2.docx (earlier draft)
├── runs-main/                 Baseline jobs: Job-Buckle, Job-Riks (inputs + logs)
├── runs-all/                  Full sweep: angles A0/A15/A45/A90, ply counts B-4/8/16, thickness C-Thin/C-thicker025
├── images/
│   ├── mode-shapes/           First 6 mode shapes per run + Riks frames
│   └── (top level)            Summary plots
└── media/                     job-riks.mp4 — Riks post-buckling animation
```

---

## Script

`scripts/classical_buckling_calculator.py`

```bash
cd 02-composite-shell-buckling/scripts
python classical_buckling_calculator.py      # needs numpy
```

Prints `Nx,cr` tables for studies A, B and C and notes the `D16 = D26 = 0`
specially-orthotropic assumption used by the closed-form formula.

---

## Results

* Critical load increases with lay-up angle up to θ ≈ 45° then falls.
* Finer ply subdivision (Study B) drives `D16/D26` towards zero, i.e. the
  laminate becomes closer to specially orthotropic, so the classical formula
  becomes more accurate.
* Thin plates buckle at lower load; the `b/t` trend from Study C matches
  classical theory.
* Riks analysis captures the unstable post-buckling branch missed by the
  linear eigenvalue solution (`media/job-riks.mp4`).

---

## How to reuse

Open `models/composite_shell_buckling.cae` for the full model. Every sweep
input deck is in `runs-all/` and the baseline in `runs-main/`. Re-submit the
`.inp` files to regenerate the `.odb` result files (not committed — see the
root `.gitignore`).
