"""
    x = 0  : CLAMPED   (C)
    x = a  : FREE      (F)
    y = 0  : SIMPLY SUPPORTED (S)
    y = b  : SIMPLY SUPPORTED (S)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  
import mpmath as mp

mp.mp.dps = 300   


a_plate = 1.2    
b_plate = 0.8      
h       = 0.010    

material = "Aluminum"
E   = 70.0e9       
nu  = 0.33         

x0 = a_plate / 2   
y0 = b_plate / 2  

P_base   = 1000.0   
load_cases = [1 * P_base, 2 * P_base, 4 * P_base] 

n_modes = 41        
grid_n  = 41        


D = E * h**3 / (12 * (1 - nu**2))     
print(f"Material: {material}  (E = {E/1e9:.1f} GPa, nu = {nu})")
print(f"Plate: a = {a_plate} m, b = {b_plate} m, h = {h*1000:.1f} mm")
print(f"Flexural rigidity D = {D:.3f} N.m")
print(f"Load point: (x0, y0) = ({x0}, {y0}) m\n")


def phi(x, beta):
    c, s = mp.cosh(beta * x), mp.sinh(beta * x)
    return [c, x * s, s, x * c]


def dphi(x, beta):
    c, s = mp.cosh(beta * x), mp.sinh(beta * x)
    return [beta * s, s + beta * x * c, beta * c, c + beta * x * s]


def d2phi(x, beta):
    c, s = mp.cosh(beta * x), mp.sinh(beta * x)
    return [beta**2 * c, 2*beta*c + beta**2*x*s, beta**2 * s, 2*beta*s + beta**2*x*c]


def d3phi(x, beta):
    c, s = mp.cosh(beta * x), mp.sinh(beta * x)
    return [beta**3 * s, 3*beta**2*s + beta**3*x*c, beta**3 * c, 3*beta**2*c + beta**3*x*s]



def solve_mode(m, a, b, x0, y0, nu, D):
    beta = mp.pi * m / b
    Fm = (2.0 / b) * mp.sin(m * mp.pi * y0 / b)   

    M = mp.zeros(8, 8)
    rhs = mp.zeros(8, 1)

    # Row 1: clamped, Y_I(0) = 0
    row = phi(0, beta)
    for j in range(4):
        M[0, j] = row[j]

    # Row 2: clamped, Y_I'(0) = 0
    row = dphi(0, beta)
    for j in range(4):
        M[1, j] = row[j]

    # Row 3: continuity of w at x0
    rI, rII = phi(x0, beta), phi(x0, beta)
    for j in range(4):
        M[2, j] = rI[j]
        M[2, 4 + j] = -rII[j]

    # Row 4: continuity of w' at x0
    rI, rII = dphi(x0, beta), dphi(x0, beta)
    for j in range(4):
        M[3, j] = rI[j]
        M[3, 4 + j] = -rII[j]

    # Row 5: continuity of w'' at x0
    rI, rII = d2phi(x0, beta), d2phi(x0, beta)
    for j in range(4):
        M[4, j] = rI[j]
        M[4, 4 + j] = -rII[j]

    # Row 6: jump in w''' at x0  ->  D*(Y_II''' - Y_I''') = Fm
    rI, rII = d3phi(x0, beta), d3phi(x0, beta)
    for j in range(4):
        M[5, j] = -rI[j]
        M[5, 4 + j] = rII[j]
    rhs[5, 0] = Fm / D

    # Row 7: free edge moment,  Y_II''(a) - nu*beta^2*Y_II(a) = 0
    r2, r0 = d2phi(a, beta), phi(a, beta)
    for j in range(4):
        M[6, 4 + j] = r2[j] - nu * beta**2 * r0[j]

    # Row 8: free edge Kirchhoff shear, Y_II'''(a) - (2-nu)*beta^2*Y_II'(a) = 0
    r3, r1 = d3phi(a, beta), dphi(a, beta)
    for j in range(4):
        M[7, 4 + j] = r3[j] - (2 - nu) * beta**2 * r1[j]

    sol = mp.lu_solve(M, rhs)
    A = [sol[i, 0] for i in range(4)]
    B = [sol[4 + i, 0] for i in range(4)]
    return beta, A, B


def Y_at(x, x0, beta, A, B):
    """Evaluate the (unit-load) Y_m(x) piecewise function."""
    coeffs = A if x <= x0 else B
    vals = phi(x, beta)
    return sum(c * v for c, v in zip(coeffs, vals))



modes = []   # list of (m, beta, A, B)
for m in range(1, n_modes + 1):
    beta, A, B = solve_mode(m, a_plate, b_plate, x0, y0, nu, D)
    modes.append((m, beta, A, B))


def deflection_unit(x, y):
    """Unit-load (P=1 N) deflection w(x,y), summing all solved modes."""
    total = mp.mpf(0)
    for m, beta, A, B in modes:
        total += Y_at(x, x0, beta, A, B) * mp.sin(beta * y)
    return total


print("Convergence of w(x0,y0) [unit load, P=1 N] with number of Fourier terms:")
running = mp.mpf(0)
check_at = [5, 10, 15, 20, 25, 30, 35, 41]
k = 0
for i, (m, beta, A, B) in enumerate(modes, start=1):
    running += Y_at(x0, x0, beta, A, B) * mp.sin(beta * y0)
    if i in check_at:
        print(f"  n_modes = {i:3d}   w = {float(running)*1000:.6f} mm")
w_unit_x0y0 = float(running)   # deflection at the load point, per unit N of load
print()


results = []
for P in load_cases:
    w_mm = w_unit_x0y0 * P * 1000.0     # convert m -> mm
    results.append({
        'Load_case': f"{P/P_base:.0f}P",
        'P_N': P,
        'w_analytical_mm': w_mm,
        'w_FE_mm': np.nan,          # <-- fill in from Abaqus at the same point
        'diff_percent': np.nan,
    })

df = pd.DataFrame(results)
print("Deflection under the point load - analytical (Levy method) vs FE:")
print(df.to_string(index=False))
df.to_csv('plate_deflection_results.csv', index=False)


xs = np.linspace(0, a_plate, grid_n)
ys = np.linspace(0, b_plate, grid_n)
W = np.zeros((grid_n, grid_n))

# reuse the pre-solved modes; evaluate Y_m at each grid x once per mode
for m, beta, A, B in modes:
    Y_vals = np.array([float(Y_at(xv, x0, beta, A, B)) for xv in xs])
    sin_vals = np.sin(float(beta) * ys)
    W += np.outer(Y_vals, sin_vals)

W_mm = W * P_base * 1000.0   # deflection [mm] for the base load case P_base

Xg, Yg = np.meshgrid(xs, ys, indexing='ij')

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(Xg, Yg, -W_mm, cmap='viridis', edgecolor='none')
# (plotted as -W so downward deflection visually points down)
ax.set_xlabel('x (m)  [0=clamped, a=free]')
ax.set_ylabel('y (m)  [0,b = simply supported]')
ax.set_zlabel('Deflection (mm)')
ax.set_title(f'Deflected shape, P = {P_base:.0f} N ({material} plate)')
fig.colorbar(surf, shrink=0.6, label='w (mm)')
fig.tight_layout()
fig.savefig('plate_deflected_shape.png', dpi=150)
plt.close(fig)

print(f"\nMax deflection on the plotted grid: {W_mm.max():.4f} mm "
      f"(thickness = {h*1000:.1f} mm - check this stays small vs. thickness "
      f"for thin-plate theory to remain valid)")


w_FE_mm = None
if w_FE_mm is not None:
    df['w_FE_mm'] = w_FE_mm
    df['diff_percent'] = 100 * (df['w_FE_mm'] - df['w_analytical_mm']) / df['w_analytical_mm']
    print("\n--- FE comparison ---")
    print(df.to_string(index=False))
    df.to_csv('plate_deflection_results.csv', index=False)
