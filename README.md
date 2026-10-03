# Abaqus Structural Analysis

A collection of **structural & solid-mechanics** studies carried out with
**Abaqus** (FEA) and cross-validated with **classical analytical theory**
(closed-form / series solutions) written in **Python** and **MATLAB**.

Each project pairs a finite-element model with a hand-derived analytical
benchmark, does its own **parametric** and/or **mesh-convergence** study, and
ships the model, scripts, data and a written report.

---

## Projects

| # | Project | Type | Key methods | Tools |
|---|---------|------|-------------|-------|
| 01 | [Single Lap Joint — Adhesive Bond](01-single-lap-joint-adhesive) | FEA + analytical | Adhesive shear/peel stresses, Goland & Reissner, parametric & mesh studies | Abaqus, Python |
| 02 | [Composite Shell Buckling](02-composite-shell-buckling) | FEA + analytical | Eigenvalue buckling, Riks (arc-length) post-buckling, Classical Lamination Theory (A–B–D / Whitney-Jones `Nx,cr`) | Abaqus, Python |
| 03 | [Plate Deflection — Levy Method](03-plate-deflection-levy) | Analytical | Kirchhoff thin-plate theory, Levy/Fourier sine series, arbitrary-precision solve | Python, MATLAB |
| 04 | [Membrane Shell of Revolution](04-membrane-shell-revolution) | Analytical | Membrane theory for shells of revolution, paraboloidal dome under self-weight | Python, MATLAB |

> **01 and 02** are **Abaqus finite-element** studies (models + result scripts).
> **03 and 04** are **analytical/theory** studies intended for comparison with
> finite-element results.

---

## Highlights

* **01 — Single Lap Joint.** Baseline model plus a full **parametric sweep**
  over adhesive thickness (0.1/0.2/0.4 mm), overlap length (6.25–50 mm) and
  adhesive modulus (500/2000/5000 MPa), a 3-level **mesh-convergence** study,
  and a Python script that overlays FEA path stresses on the Goland & Reissner
  closed-form solution.
* **02 — Composite Shell Buckling.** Baseline **Buckle** + **Riks** runs, a
  lay-up **angle sweep** (0/15/45/90°), a **ply-count sweep** (4/8/16), a
  **b/t thickness sweep**, mode-shape images and a post-buckling **animation**,
  all backed by a from-scratch CLT buckling calculator.
* **03 — Plate Deflection.** Levy-method solution for a C–F–S–S rectangular
  plate under a point load, in both Python (`mpmath`) and MATLAB (`vpa`).
* **04 — Membrane Dome.** Membrane forces/stresses in a paraboloidal dome under
  self-weight (steel and concrete variants), in Python and MATLAB.

---

## Repository structure

```
abaqus-structural-analysis/
├── README.md                       ← you are here
├── .gitignore                      ← excludes large/regenerable solver output
├── 01-single-lap-joint-adhesive/   ← project 01 (+ its own README)
├── 02-composite-shell-buckling/    ← project 02 (+ its own README)
├── 03-plate-deflection-levy/       ← project 03 (+ its own README)
├── 04-membrane-shell-revolution/   ← project 04 (+ its own README)
├── docs/                           ← combined written reports (projects 03 & 04)
└── misc/folder3-duplicates/        ← preserved extra copies (see misc/README.md)
```

Every project folder contains its **own `README.md`** describing the theory,
the folder contents and how to run the scripts.

---

## Getting started

### Abaqus finite-element projects (01, 02)

Open the committed model file in **Abaqus/CAE**:

```
01-single-lap-joint-adhesive/models/single_lap_joint.cae
02-composite-shell-buckling/models/composite_shell_buckling.cae
```

To re-run a job, submit the corresponding input deck, e.g.

```bash
abaqus job=Job-Riks  input=02-composite-shell-buckling/runs-all/Job-Riks.inp
abaqus viewer odb=Job-Riks
```

The `.inp` decks for every parametric/sweep case are kept in
`.../abaqus-jobs/` (project 01) and `.../runs-all/` + `.../runs-main/`
(project 02).

### Analytical theory projects (03, 04)

```bash
cd 03-plate-deflection-levy/python   && python project1_plate_deflection.py
cd 04-membrane-shell-revolution/python && python project2_membrane_dome_steel.py
```

The MATLAB counterparts (`.m`) need **MATLAB** (project 03 also needs the
*Symbolic Math Toolbox* for high-precision `vpa`).

---

## Requirements

| Tool | Used by |
|---|---|
| **Abaqus** (CAE / Standard) | 01, 02 |
| **Python 3** + `numpy`, `pandas`, `matplotlib` | 01, 02, 03, 04 |
| **`mpmath`** | 03 (arbitrary-precision plate solve) |
| **MATLAB** (+ Symbolic Math Toolbox for 03) | 03, 04 |

---

## A note on large files

Large, **binary, fully-regenerable** solver outputs are intentionally **not**
committed (see [`.gitignore`](.gitignore)) — `*.odb`, `*.fil`, `*.msg`,
`*.dat`, `*.sta`, `*.com`, `*.prt`, `*.log`, `*.ipm`, etc. They are the bulk
of the raw data but can be recreated at any time by running the committed
`.cae` / `.inp` files.

Committed (kept): `.cae`, `.jnl`, `.inp`, `.py`, `.m`, `.csv`, `.png`,
`.docx` and `.mp4` — i.e. models, inputs, scripts, extracted data, figures
and reports.

---

## Author

**Reza Nadimi** — <rezaxr14@gmail.com>

---

## License

Released for academic / portfolio reference. See individual project reports
for any course-specific terms.

