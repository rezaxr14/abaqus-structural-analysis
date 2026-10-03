"""

analytically for the meridian  f(r) = a*r^2.


    N_phi   = - W(r) / (2*pi*r*sin(phi))
    W(r)    = total self-weight of the shell above the parallel at r

    N_phi/R1 + N_theta/R2 = -Z ,      Z = q0*cos(phi)
    N_theta = -R2*q0*cos(phi) - N_phi*(R2/R1)

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


R      = 10.0     
H      = 4.0      
t      = 0.02    
gamma  = 77.0    
E      = 200e6    
nu     = 0.30     
r_edge_ratio = 0.90                                            
r_eval_ratio = 0.60                       
n_pts = 500             


a  = H / R**2      # parabola coefficient:  z = a*r^2
q0 = gamma * t      # self-weight per unit AREA of the mid-surface

print(f"Parabola coefficient   a  = {a:.5f}  (1/m)")
print(f"Self-weight per area   q0 = {q0:.4f} kN/m^2\n")



def shell_geometry(r, a):
    """Return phi, sin(phi), cos(phi), R1, R2 for meridian f(r)=a*r^2."""
    fp   = 2 * a * r                      # f'(r) = tan(phi)
    phi  = np.arctan(fp)
    sphi = fp / np.sqrt(1 + fp**2)
    cphi = 1.0 / np.sqrt(1 + fp**2)
    R1   = (1 + fp**2)**1.5 / (2 * a)
    R2   = np.divide(r, sphi, out=np.zeros_like(r), where=sphi != 0)
    return phi, sphi, cphi, R1, R2


def weight_above(r, a, q0):
    """W(r): closed-form self-weight of the shell above radius r."""
    return (np.pi * q0 / (6 * a**2)) * ((1 + 4 * a**2 * r**2)**1.5 - 1)


def membrane_forces(r, a, q0):
    """N_phi, N_theta at radius r (r may be an array; r=0 handled below)."""
    phi, sphi, cphi, R1, R2 = shell_geometry(r, a)
    W = weight_above(r, a, q0)

    denom = 2 * np.pi * r * sphi
    N_phi = np.divide(-W, denom, out=np.zeros_like(r), where=denom != 0)
    N_theta = -R2 * q0 * cphi - N_phi * (R2 / R1)
    return phi, N_phi, N_theta, R1, R2


# analytic apex limit (r -> 0), used as the value at r=0 and as a sanity check
N_phi0 = -q0 / (4 * a)
N_theta0 = N_phi0                     # N_phi = N_theta at the apex
print(f"Apex values (r=0):   N_phi = N_theta = {N_phi0:.4f} kN/m\n")


r_min = 1e-4 * R                       # avoid the 0/0 singularity at the apex
r = np.linspace(r_min, r_edge_ratio * R, n_pts)

phi, N_phi, N_theta, R1, R2 = membrane_forces(r, a, q0)


sigma_phi = N_phi / t         # [kN/m^2 = kPa]
sigma_theta = N_theta / t     # [kN/m^2 = kPa]


fig, axes = plt.subplots(2, 1, figsize=(7, 8))

axes[0].plot(r / R, N_phi, 'b-', linewidth=1.6, label=r'$N_\phi$ (meridional)')
axes[0].plot(r / R, N_theta, 'r-', linewidth=1.6, label=r'$N_\theta$ (hoop)')
axes[0].set_xlabel('r / R')
axes[0].set_ylabel('Membrane force (kN/m)')
axes[0].set_title('Membrane forces along the meridian - steel dome, self-weight only')
axes[0].legend()
axes[0].grid(True)

axes[1].plot(np.degrees(phi), N_phi, 'b-', linewidth=1.6, label=r'$N_\phi$')
axes[1].plot(np.degrees(phi), N_theta, 'r-', linewidth=1.6, label=r'$N_\theta$')
axes[1].set_xlabel(r'$\phi$ (deg)')
axes[1].set_ylabel('Membrane force (kN/m)')
axes[1].legend()
axes[1].grid(True)

fig.tight_layout()
fig.savefig('membrane_forces_steel.png', dpi=150)
plt.close(fig)


ratios_report = np.arange(0, r_edge_ratio + 1e-9, 0.1)
r_report = ratios_report * R

phi_r, Nphi_r, Ntheta_r, R1_r, R2_r = membrane_forces(r_report, a, q0)

# overwrite the apex row (index 0, r=0) with the analytic limit
Nphi_r[0] = N_phi0
Ntheta_r[0] = N_theta0

sigphi_r = Nphi_r / t
sigtheta_r = Ntheta_r / t

df = pd.DataFrame({
    'r_over_R':        ratios_report,
    'r_m':             r_report,
    'phi_deg':         np.degrees(phi_r),
    'N_phi_kNm':       Nphi_r,
    'N_theta_kNm':     Ntheta_r,
    'sigma_phi_kPa':   sigphi_r,
    'sigma_theta_kPa': sigtheta_r,
})

pd.set_option('display.float_format', lambda x: f'{x:0.4f}')
print("Membrane theory results along the meridian (steel dome):")
print(df.to_string(index=False))

df.to_csv('membrane_theory_results_steel.csv', index=False)

r_eval = r_eval_ratio * R
phi_eval, N_phi_eval, N_theta_eval, R1_eval, R2_eval = membrane_forces(
    np.array([r_eval]), a, q0)
phi_eval = phi_eval[0]
N_phi_eval = N_phi_eval[0]
N_theta_eval = N_theta_eval[0]

sigma_phi_eval = N_phi_eval / t
sigma_theta_eval = N_theta_eval / t

print(f"\n--- Comparison point (away from support edge, r/R = {r_eval_ratio:.2f}) ---")
print(f"r = {r_eval:.3f} m ,  phi = {np.degrees(phi_eval):.2f} deg")
print(f"N_phi   (membrane theory) = {N_phi_eval:8.3f} kN/m  ->  sigma_phi   = {sigma_phi_eval:8.2f} kPa")
print(f"N_theta (membrane theory) = {N_theta_eval:8.3f} kN/m  ->  sigma_theta = {sigma_theta_eval:8.2f} kPa")
