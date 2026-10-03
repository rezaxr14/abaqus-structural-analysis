
%           analytically for the meridian  f(r) = a*r^2.
%     Vertical equilibrium of the cap above radius r:
%         N_phi   = - W(r) / (2*pi*r*sin(phi))
%         W(r)    = total self-weight of the shell above the parallel at r
%     Equilibrium normal to the shell surface:
%         N_phi/R1 + N_theta/R2 = -Z ,      Z = q0*cos(phi)
%         N_theta = -R2*q0*cos(phi) - N_phi*(R2/R1)

clear; clc; close all;


R      = 5.0;     
H      = 8.0;       
t      = 0.15;      
gamma  = 23.544;       
r_edge_ratio = 0.90;                                             
r_eval_ratio = 0.60;                        
n_pts  = 500;         

a  = H / R^2;        
q0 = gamma * t;      
fprintf('Parabola coefficient   a  = %.5f  (1/m)\n', a);
fprintf('Self-weight per area   q0 = %.4f kN/m^2\n\n', q0);

r_min = 1e-4 * R;  
r     = linspace(r_min, r_edge_ratio*R, n_pts).';

%     f(r)  = a*r^2
%     f'(r) = 2*a*r = tan(phi)

fp   = 2*a*r;                           
phi  = atan(fp);                     
sphi = fp ./ sqrt(1 + fp.^2);             % sin(phi)
cphi = 1  ./ sqrt(1 + fp.^2);             % cos(phi)

R1 = (1 + fp.^2).^1.5 / (2*a);            % meridional radius of curvature
R2 = r ./ sphi;                           % second principal radius


W = (pi*q0/(6*a^2)) .* ( (1 + 4*a^2*r.^2).^1.5 - 1 );


N_phi   = - W ./ (2*pi*r.*sphi);
N_theta = - R2.*q0.*cphi - N_phi.*(R2./R1);

% analytic apex limit (r -> 0), used only as a reference / sanity check
N_phi0   = -q0/(4*a);
N_theta0 = N_phi0;                        % N_phi = N_theta at the apex
fprintf('Apex values (r=0):   N_phi = N_theta = %.4f kN/m\n\n', N_phi0);

sigma_phi   = N_phi   / t;     % [kN/m^2 = kPa]
sigma_theta = N_theta / t;     % [kN/m^2 = kPa]

figure('Color','w');
subplot(2,1,1)
plot(r/R, N_phi, 'b-', 'LineWidth', 1.6); hold on
plot(r/R, N_theta, 'r-', 'LineWidth', 1.6);
xlabel('r / R'); ylabel('Membrane force (kN/m)');
legend('N_\phi (meridional)','N_\theta (hoop)','Location','best');
title('Membrane forces along the meridian - self-weight only'); grid on

subplot(2,1,2)
plot(rad2deg(phi), N_phi, 'b-', 'LineWidth', 1.6); hold on
plot(rad2deg(phi), N_theta, 'r-', 'LineWidth', 1.6);
xlabel('\phi  (deg)'); ylabel('Membrane force (kN/m)');
legend('N_\phi','N_\theta','Location','best'); grid on

ratios_report = (0:0.1:r_edge_ratio).';
r_report      = ratios_report * R;

phi_r  = atan(2*a*r_report);
sphi_r = sin(phi_r);
cphi_r = cos(phi_r);
R1_r   = (1 + (2*a*r_report).^2).^1.5 / (2*a);
R2_r   = r_report ./ max(sphi_r, eps);
W_r    = (pi*q0/(6*a^2)) .* ((1 + 4*a^2*r_report.^2).^1.5 - 1);

Nphi_r      = -W_r ./ (2*pi*max(r_report,eps).*max(sphi_r,eps));
Nphi_r(1)   = N_phi0;                      % use analytic apex value
Ntheta_r    = -R2_r.*q0.*cphi_r - Nphi_r.*(R2_r./R1_r);
Ntheta_r(1) = N_theta0;

sigphi_r   = Nphi_r   / t;
sigtheta_r = Ntheta_r / t;

T = table(ratios_report, r_report, rad2deg(phi_r), Nphi_r, Ntheta_r, ...
          sigphi_r, sigtheta_r, ...
    'VariableNames', {'r_over_R','r_m','phi_deg','N_phi_kNm','N_theta_kNm', ...
                       'sigma_phi_kPa','sigma_theta_kPa'});

disp('Membrane theory results along the meridian:');
disp(T)

writetable(T, 'membrane_theory_results.csv');

r_eval    = r_eval_ratio * R;
phi_eval  = atan(2*a*r_eval);
sphi_eval = sin(phi_eval);
cphi_eval = cos(phi_eval);
R1_eval   = (1 + (2*a*r_eval)^2)^1.5 / (2*a);
R2_eval   = r_eval / sphi_eval;
W_eval    = (pi*q0/(6*a^2)) * ((1 + 4*a^2*r_eval^2)^1.5 - 1);

Nphi_eval     = -W_eval / (2*pi*r_eval*sphi_eval);
Ntheta_eval   = -R2_eval*q0*cphi_eval - Nphi_eval*(R2_eval/R1_eval);
sigphi_eval   = Nphi_eval   / t;
sigtheta_eval = Ntheta_eval / t;

fprintf('\n--- Comparison point (away from support edge, r/R = %.2f) ---\n', r_eval_ratio);
fprintf('r = %.3f m ,  phi = %.2f deg\n', r_eval, rad2deg(phi_eval));
fprintf('N_phi   (membrane theory) = %8.3f kN/m  ->  sigma_phi   = %8.2f kPa\n', Nphi_eval, sigphi_eval);
fprintf('N_theta (membrane theory) = %8.3f kN/m  ->  sigma_theta = %8.2f kPa\n', Ntheta_eval, sigtheta_eval);

% --- paste your FE (finite element) results here to compare ------------
N_phi_FE   = [];   % <-- fill in from SAP2000 / ANSYS / ABAQUS
N_theta_FE = [];
if ~isempty(N_phi_FE) && ~isempty(N_theta_FE)
    fprintf('\n--- FE comparison ---\n');
    fprintf('N_phi:   FE = %.3f kN/m , difference = %.2f %%\n', ...
        N_phi_FE, 100*(N_phi_FE - Nphi_eval)/Nphi_eval);
    fprintf('N_theta: FE = %.3f kN/m , difference = %.2f %%\n', ...
        N_theta_FE, 100*(N_theta_FE - Ntheta_eval)/Ntheta_eval);
end
