# 01 · Single Lap Joint — Adhesive Bond Analysis

Numerical study of an **adhesive-bonded single-lap joint (SLJ)** using Abaqus
finite-element analysis, validated against the classical **Goland & Reissner**
closed-form solution for the shear and peel stress distribution along the
overlap. Includes a full **parametric study** and a **mesh-convergence study**.

---

## Physics / method

An adhesive layer is a stress-transfer medium: the load is carried by
**shear** in the middle of the overlap and by **peel** at the ends of the
overlap, where both stresses peak. The joint is modelled in Abaqus and the
stress resultants are extracted along a path through the bond line. The same
quantities are predicted analytically with the Goland & Reissner (1944)
formulae for a lap joint, and the two are compared.

**Baseline parameters used in the study**

| Parameter | Symbol | Value |
|---|---|---|
| Adherend modulus | E | 70 000 MPa (aluminium) |
| Adherend Poisson ratio | ν | 0.33 |
| Adherend thickness | t | 1.6 mm |
| Adhesive modulus | Eₐ | 2 000 MPa |
| Adhesive Poisson ratio | νₐ | 0.35 |
| Adhesive thickness | tₐ | 0.2 mm |
| Overlap length | L₀ | 12.5 mm |
| Joint width | b | 25 mm |
| Applied load | P | 2 000 N |

Sweeps: adhesive thickness (tₐ = 0.1 / 0.2 / 0.4 mm), overlap length
(L₀ = 6.25 / 12.5 / 25 / 50 mm) and adhesive modulus (Eₐ = 500 / 2 000 / 5 000 MPa).

---

## Folder layout

```
01-single-lap-joint-adhesive/
├── models/                       Baseline Abaqus model (single_lap_joint.cae / .jnl)
├── images/                       Baseline set-up figures (assembly, loads, BCs, coupling, ...)
├── reports/                      SLJ_Project_Report_Persian_FINAL_v2.docx
├── mesh-convergence/
│   ├── models/                   1.cae / 1lo50.cae
│   ├── images/                   S & U contour plots (mesh2 / mesh8 / Lo50)
│   └── abaqus-jobs/              Job inputs (.inp) + path-stress CSVs
└── parametric-study/
    ├── models/                   Parametric Abaqus model
    ├── data/                     Path stress exports for every sweep case
    ├── scripts/                  slj_goland_reissner_compare.py
    ├── images/                   S & U contour plots for every sweep case
    └── abaqus-jobs/              Parametric job inputs (.inp)
```

---

## Scripts

`parametric-study/scripts/slj_goland_reissner_compare.py`

* Loads every path-stress export (`*.csv`) from `../data`.
* Evaluates the **Goland & Reissner** shear and peel stresses for the
  matching parameter set.
* Produces overlay plots (GR vs FEA) for each sweep and a
  mesh-convergence table/figure.

```bash
cd 01-single-lap-joint-adhesive/parametric-study/scripts
python slj_goland_reissner_compare.py      # needs numpy + matplotlib
```

> Note: the peel comparison has a **known limitation** — the original
> Goland & Reissner peel derivation carries a documented plane-stress /
> plane-strain inconsistency (later corrected by Tsai et al.), so the
> classical peel value is much lower than FEA for this thin/stiff adhesive.
> This is discussed explicitly in the script and report.

---

## Results

* Shear-stress peaks are predicted well by the analytical model; peel is
  under-predicted by classical theory (see note above).
* Mesh convergence shows the peak stresses stabilise as the number of
  through-thickness elements in the adhesive increases (2 → 4 → 8 elements).
* Longer overlaps and thicker/softer adhesives reduce the peak stresses.

---

## How to reuse

Open `models/single_lap_joint.cae` in Abaqus/CAE to inspect or re-run the
model. Regenerate the `.odb` by submitting the committed `.inp` files in
`mesh-convergence/abaqus-jobs` and `parametric-study/abaqus-jobs`
(solver outputs are not stored in this repository — see the root `.gitignore`).
