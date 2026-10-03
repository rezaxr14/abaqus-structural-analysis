
import numpy as np

# ---------------------------------------------------------------------------
# 1. Lamina properties (typical UD carbon/epoxy - EDIT to match your material)
# ---------------------------------------------------------------------------
E1 = 135000.0   # MPa
E2 = 10000.0    # MPa
G12 = 5000.0    # MPa
NU12 = 0.30

# ---------------------------------------------------------------------------
# 2. Core CLT functions
# ---------------------------------------------------------------------------


def reduced_stiffness(E1, E2, nu12, G12):
    """In-plane reduced stiffnesses Q11, Q22, Q12, Q66 for a UD lamina."""
    nu21 = nu12 * E2 / E1
    denom = 1.0 - nu12 * nu21
    Q11 = E1 / denom
    Q22 = E2 / denom
    Q12 = nu12 * E2 / denom
    Q66 = G12
    return Q11, Q22, Q12, Q66


def Qbar(Q11, Q22, Q12, Q66, theta_deg):
    """Transformed reduced stiffness matrix for a ply at angle theta_deg
    (measured from the laminate/loading X-axis)."""
    th = np.radians(theta_deg)
    m, n = np.cos(th), np.sin(th)
    Q11b = Q11 * m**4 + 2 * (Q12 + 2 * Q66) * m**2 * n**2 + Q22 * n**4
    Q22b = Q11 * n**4 + 2 * (Q12 + 2 * Q66) * m**2 * n**2 + Q22 * m**4
    Q12b = (Q11 + Q22 - 4 * Q66) * m**2 * n**2 + Q12 * (m**4 + n**4)
    Q66b = (Q11 + Q22 - 2 * Q12 - 2 * Q66) * m**2 * n**2 + Q66 * (m**4 + n**4)
    Q16b = (Q11 - Q12 - 2 * Q66) * m**3 * n - (Q22 - Q12 - 2 * Q66) * m * n**3
    Q26b = (Q11 - Q12 - 2 * Q66) * m * n**3 - (Q22 - Q12 - 2 * Q66) * m**3 * n
    return np.array([[Q11b, Q12b, Q16b],
                      [Q12b, Q22b, Q26b],
                      [Q16b, Q26b, Q66b]])


def abd_matrix(angles_deg, t_ply, E1=E1, E2=E2, nu12=NU12, G12=G12):
    """Full A, B, D laminate stiffness matrices for a laminate with
    equal-thickness plies, given bottom-to-top ply angles."""
    Q11, Q22, Q12, Q66 = reduced_stiffness(E1, E2, nu12, G12)
    n_ply = len(angles_deg)
    t_total = n_ply * t_ply
    z = [-t_total / 2.0 + i * t_ply for i in range(n_ply + 1)]
    A = np.zeros((3, 3))
    B = np.zeros((3, 3))
    D = np.zeros((3, 3))
    for k, th in enumerate(angles_deg):
        Qb = Qbar(Q11, Q22, Q12, Q66, th)
        A += Qb * (z[k + 1] - z[k])
        B += 0.5 * Qb * (z[k + 1] ** 2 - z[k] ** 2)
        D += (1.0 / 3.0) * Qb * (z[k + 1] ** 3 - z[k] ** 3)
    return A, B, D, t_total


def classical_Nx_cr(D, a, b, m_max=15, n=1):
    """Whitney/Jones closed-form Nx,cr for a simply-supported specially-
    orthotropic plate, minimized over half-wave number m (n=1 fixed)."""
    D11, D22, D12, D66 = D[0, 0], D[1, 1], D[0, 1], D[2, 2]
    best_m, best_N = None, None
    for m in range(1, m_max + 1):
        N = (np.pi ** 2 *
             (D11 * (m / a) ** 4 + 2 * (D12 + 2 * D66) * (m / a) ** 2 * (n / b) ** 2
              + D22 * (n / b) ** 4) / (m / a) ** 2)
        if best_N is None or N < best_N:
            best_N, best_m = N, m
    return best_N, best_m


def pm_angle_stack(theta, n_total):
    """Symmetric, balanced [+theta/-theta] stack with n_total plies (even),
    alternating sign every single ply -- the closest approach to
    specially-orthotropic behavior at a given ply count."""
    half = n_total // 2
    half_stack = [theta if i % 2 == 0 else -theta for i in range(half)]
    return half_stack + half_stack[::-1]


# ---------------------------------------------------------------------------
# 3. Baseline geometry
# ---------------------------------------------------------------------------
a = 300.0   # mm, loaded direction
b = 150.0   # mm

# ---------------------------------------------------------------------------
# 4. Study A: layup angle sweep, [+/-theta]_2s, 8 plies, t_total = 1 mm
# ---------------------------------------------------------------------------
print("Study A: layup angle sweep  [+/-theta]_2s, 8 plies, t_total = 1.0 mm")
print(f"{'theta (deg)':>12} {'Nx,cr classical (N/mm)':>24} {'critical m':>11}")
for theta in [0, 15, 30, 45, 60, 75, 90]:
    angles = [theta] * 8 if theta in (0, 90) else pm_angle_stack(theta, 8)
    A, B, D, t = abd_matrix(angles, t_ply=0.125)
    Ncr, m = classical_Nx_cr(D, a, b)
    print(f"{theta:12d} {Ncr:24.2f} {m:11d}")

# ---------------------------------------------------------------------------
# 5. Study B: number of plies, [+/-45] finely alternating, t_total = 1 mm fixed
# ---------------------------------------------------------------------------
print("\nStudy B: number of layers, [+/-45] finely alternating, t_total = 1.0 mm")
print(f"{'N plies':>8} {'t_ply (mm)':>12} {'Nx,cr classical (N/mm)':>24} {'critical m':>11}")
for N in [4, 8, 16, 32]:
    t_ply = 1.0 / N
    angles = pm_angle_stack(45, N)
    A, B, D, t = abd_matrix(angles, t_ply=t_ply)
    Ncr, m = classical_Nx_cr(D, a, b)
    print(f"{N:8d} {t_ply:12.5f} {Ncr:24.2f} {m:11d}")

# ---------------------------------------------------------------------------
# 6. Study C: b/t sweep, quasi-isotropic [0/45/-45/90]s, a=300 b=150 fixed
# ---------------------------------------------------------------------------
print("\nStudy C: b/t sweep, quasi-isotropic [0/45/-45/90]_s, a=300 b=150 fixed")
angles_qi = [0, 45, -45, 90, 90, -45, 45, 0]
print(f"{'t_total (mm)':>13} {'t_ply (mm)':>12} {'b/t':>8} "
      f"{'Nx,cr classical (N/mm)':>24} {'critical m':>11}")
for t_total in [0.5, 1.0, 2.0, 4.0, 8.0]:
    t_ply = t_total / 8.0
    A, B, D, t = abd_matrix(angles_qi, t_ply=t_ply)
    Ncr, m = classical_Nx_cr(D, a, b)
    print(f"{t_total:13.2f} {t_ply:12.4f} {b / t_total:8.2f} {Ncr:24.2f} {m:11d}")

# ---------------------------------------------------------------------------
# 7. Reminder about D16/D26 (bend-twist coupling)
# ---------------------------------------------------------------------------
print("\nNote: the Nx,cr formula above assumes D16 = D26 = 0 (specially")
print("orthotropic). To check how good that assumption is for any case,")
print("inspect D[0,2] and D[1,2] from abd_matrix() and compare their")
print("magnitude to D11/D22/D66 -- this is exactly what shrinks in Study B.")
