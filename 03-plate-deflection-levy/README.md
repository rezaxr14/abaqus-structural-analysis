# 03 · Rectangular Plate Deflection (Levy Method)

Analytical (Levy) solution for the **deflection of a rectangular thin plate
under a concentrated point load**, with a mixed **clamped / free / simply-
supported** boundary. Implemented in both **Python** and **MATLAB**.

---

## Physics / method

Classical **Kirchhoff thin-plate theory**. The boundary conditions are:

| Edge | Condition |
|---|---|
| x = 0 | **Clamped** (C) |
| x = a | **Free** (F) |
| y = 0 | Simply supported (S) |
| y = b | Simply supported (S) |

Because both edges at y = 0 and y = b are simply supported, the deflection is
expanded as a **Fourier sine series** in y (Levy's method). Each term satisfies
a fourth-order ODE in x with repeated characteristic roots; the domain is split
at the load line into two regions and the 8 unknowns per mode are found from
clamping + continuity + load-jump + free-edge conditions.

The solve is carried out in **arbitrary precision** (Python `mpmath`, MATLAB
`vpa`) because `cosh`/`sinh` of large arguments at higher Fourier modes span an
enormous dynamic range and would lose precision or overflow in double
precision.

Two material variants are provided: **aluminium** and **titanium**.

---

## Folder layout

```
03-plate-deflection-levy/
├── python/     project1_plate_deflection.py           (mpmath implementation)
├── matlab/     project1_plate_deflection_titanium.m   (Symbolic Math Toolbox / vpa)
├── data/       plate_deflection_results*.csv          (deflection at the load point)
├── images/     plate_deflected_shape.png              (3-D deflected surface)
└── reports/    (combined project-3 & 4 reports are in /docs)
```

---

## Scripts

```bash
cd 03-plate-deflection-levy/python
python project1_plate_deflection.py     # needs numpy, pandas, matplotlib, mpmath
```

```matlab
% MATLAB — requires the Symbolic Math Toolbox
cd 03-plate-deflection-levy/matlab
project1_plate_deflection_titanium
```

Both versions print a convergence table (deflection vs number of Fourier
terms), the deflection for the required load cases (P, 2P, 4P) and write
`plate_deflection_results.csv`. A placeholder is provided to paste finite-
element results for a direct analytical-vs-FE comparison.

---

## Results

* Convergence with the number of Fourier/Levy terms is rapid.
* Load-case deflections scale linearly with P, as expected for thin-plate
  theory.
* The maximum deflection is checked against the plate thickness to confirm
  thin-plate validity. Titanium (stiffer) deflects less than aluminium.

---

## How to reuse

Edit the input block at the top of either script (geometry, material, load,
number of modes) to study a different plate; the method is general.
