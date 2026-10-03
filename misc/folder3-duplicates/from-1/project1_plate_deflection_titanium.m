
%      x = 0  : CLAMPED               (C)
%      x = a  : FREE                  (F)
%      y = 0  : SIMPLY SUPPORTED      (S)
%      y = b  : SIMPLY SUPPORTED      (S)
%      w(x,y) = sum_{m=1}^{N}  Y_m(x) * sin(m*pi*y/b)

%      q_m(x) = (2P/b) * sin(m*pi*y0/b) * delta(x-x0)
%
%      Y_m'''' - 2*beta_m^2*Y_m'' + beta_m^4*Y_m = q_m(x)/D

%      Y_m(x) = A1*cosh(bx) + A2*x*sinh(bx) + A3*sinh(bx) + A4*x*cosh(bx)

%      Clamped  @ x=0  : Y(0)=0 , Y'(0)=0                     
%      Continuity @ x0 : w, w', w'' continuous               
%      Load jump @ x0  : D*(Y_II'''-Y_I''') = (2P/b)*sin(m pi y0/b)  
%      Free @ x=a      : M_x=0 , V_x (Kirchhoff shear)=0       


clear; clc; close all;
digits(300);   % 300 significant decimal digits 


a_plate = 1.5;        % plate length, clamped (x=0) to free (x=a) edge   
b_plate = 1.0;        % plate width, between the two simply supported edges 
h       = 0.015;      % plate thickness                              

material = 'Titanium alloy (Ti-6Al-4V)';
E   = 110e9;          % elastic modulus                     
nu  = 0.31;           % Poisson's ratio

x0 = a_plate/2;       % load location (x) - assumed at plate centre   
y0 = b_plate/2;       % load location (y) - assumed at plate centre    

P_base     = 1500;                    % base load magnitude            
load_cases = [1 2 4] * P_base;        % P, 2P, 4P (assignment requirement)

n_modes = 20;          % number of Fourier (Levy) terms 
grid_n  = 21;          % grid resolution for the deflected-shape plot only

D = E*h^3/(12*(1-nu^2));      % flexural rigidity of the plate          [N.m]
fprintf('Material: %s (E = %.1f GPa, nu = %.2f)\n', material, E/1e9, nu);
fprintf('Plate: a = %.2f m, b = %.2f m, h = %.1f mm\n', a_plate, b_plate, h*1000);
fprintf('Flexural rigidity D = %.3f N.m\n', D);
fprintf('Load point: (x0, y0) = (%.3f, %.3f) m\n\n', x0, y0);

av = vpa(a_plate); bv = vpa(b_plate);
x0v = vpa(x0);      y0v = vpa(y0);
nuv = vpa(nu);       Dv = vpa(D);

beta_all = cell(n_modes,1);
A_all    = cell(n_modes,1);
B_all    = cell(n_modes,1);

for m = 1:n_modes
    [beta, Avec, Bvec] = solve_mode(m, av, bv, x0v, y0v, nuv, Dv);
    beta_all{m} = beta;
    A_all{m}    = Avec;
    B_all{m}    = Bvec;
end

fprintf('Convergence of w(x0,y0) [unit load, P=1 N] with number of Fourier terms:\n');
running = vpa(0);
check_at = [3 5 8 10 15 20];
for m = 1:n_modes
    running = running + Y_at(x0v, x0v, beta_all{m}, A_all{m}, B_all{m}) * sin(beta_all{m}*y0v);
    if any(m == check_at)
        fprintf('  n_modes = %2d   w = %.6f mm\n', m, double(running)*1000);
    end
end
w_unit_x0y0 = double(running);    % deflection at load point, per unit N of load
fprintf('\n');

Load_case       = strings(numel(load_cases),1);
P_N              = zeros(numel(load_cases),1);
w_analytical_mm  = zeros(numel(load_cases),1);
w_FE_mm          = nan(numel(load_cases),1);     % <-- fill in from Abaqus

for i = 1:numel(load_cases)
    P = load_cases(i);
    Load_case(i)      = sprintf('%.0fP', P/P_base);
    P_N(i)            = P;
    w_analytical_mm(i) = w_unit_x0y0 * P * 1000;   % m -> mm
end

T = table(Load_case, P_N, w_analytical_mm, w_FE_mm);
disp('Deflection under the point load - analytical (Levy method) vs FE:');
disp(T)
writetable(T, 'plate_deflection_results_titanium.csv');

xs = linspace(0, a_plate, grid_n);
ys = linspace(0, b_plate, grid_n);
W  = zeros(grid_n, grid_n);

for m = 1:n_modes
    beta = beta_all{m}; Avec = A_all{m}; Bvec = B_all{m};
    Yx = zeros(1, grid_n);
    for i = 1:grid_n
        Yx(i) = double(Y_at(vpa(xs(i)), x0v, beta, Avec, Bvec));
    end
    sinY = sin(double(beta) * ys);
    W = W + Yx.' * sinY;             % outer product, accumulate mode contribution
end

W_mm = W * P_base * 1000;            % deflection [mm] for the base load case

[Xg, Yg] = ndgrid(xs, ys);
figure('Color','w');
surf(Xg, Yg, -W_mm);                 % plotted as -W so downward deflection points down
shading interp; colorbar
xlabel('x (m)  [0 = clamped, a = free]');
ylabel('y (m)  [0,b = simply supported]');
zlabel('Deflection (mm)');
title(sprintf('Deflected shape, P = %.0f N (%s plate)', P_base, material));

fprintf('\nMax deflection on the plotted grid: %.4f mm (thickness = %.1f mm - ', ...
    max(W_mm(:)), h*1000);
fprintf('check this stays small vs. thickness for thin-plate theory to remain valid)\n');

function v = phi_fun(x, beta)
    c = cosh(beta*x); s = sinh(beta*x);
    v = [c, x*s, s, x*c];
end

function v = dphi_fun(x, beta)
    c = cosh(beta*x); s = sinh(beta*x);
    v = [beta*s, s + beta*x*c, beta*c, c + beta*x*s];
end

function v = d2phi_fun(x, beta)
    c = cosh(beta*x); s = sinh(beta*x);
    v = [beta^2*c, 2*beta*c + beta^2*x*s, beta^2*s, 2*beta*s + beta^2*x*c];
end

function v = d3phi_fun(x, beta)
    c = cosh(beta*x); s = sinh(beta*x);
    v = [beta^3*s, 3*beta^2*s + beta^3*x*c, beta^3*c, 3*beta^2*c + beta^3*x*s];
end

function [beta, Avec, Bvec] = solve_mode(m, a, b, x0, y0, nu, D)
    beta = vpa(pi)*m/b;
    Fm   = (2/b) * sin(vpa(m)*vpa(pi)*y0/b);   % forcing amplitude, unit load P=1N

    M   = vpa(zeros(8,8));
    rhs = vpa(zeros(8,1));

    % Row 1: clamped, Y_I(0) = 0
    row = phi_fun(vpa(0), beta);
    M(1,1:4) = row;

    % Row 2: clamped, Y_I'(0) = 0
    row = dphi_fun(vpa(0), beta);
    M(2,1:4) = row;

    % Row 3: continuity of w at x0
    rI = phi_fun(x0, beta);
    M(3,1:4) = rI;  M(3,5:8) = -rI;

    % Row 4: continuity of w' at x0
    rI = dphi_fun(x0, beta);
    M(4,1:4) = rI;  M(4,5:8) = -rI;

    % Row 5: continuity of w'' at x0
    rI = d2phi_fun(x0, beta);
    M(5,1:4) = rI;  M(5,5:8) = -rI;

    % Row 6: jump in w''' at x0  ->  D*(Y_II'''-Y_I''') = Fm
    rI = d3phi_fun(x0, beta);
    M(6,1:4) = -rI; M(6,5:8) = rI;
    rhs(6) = Fm/D;

    % Row 7: free edge moment, Y_II''(a) - nu*beta^2*Y_II(a) = 0
    r2 = d2phi_fun(a, beta); r0 = phi_fun(a, beta);
    M(7,5:8) = r2 - nu*beta^2*r0;

    % Row 8: free edge Kirchhoff shear, Y_II'''(a) - (2-nu)*beta^2*Y_II'(a) = 0
    r3 = d3phi_fun(a, beta); r1 = dphi_fun(a, beta);
    M(8,5:8) = r3 - (2-nu)*beta^2*r1;

    sol  = M \ rhs;
    Avec = sol(1:4);
    Bvec = sol(5:8);
end

function val = Y_at(x, x0, beta, Avec, Bvec)
    % Evaluate the (unit-load) Y_m(x) piecewise function.
    if x <= x0
        coeffs = Avec;
    else
        coeffs = Bvec;
    end
    v = phi_fun(x, beta);
    val = sum(coeffs.' .* v);
end
