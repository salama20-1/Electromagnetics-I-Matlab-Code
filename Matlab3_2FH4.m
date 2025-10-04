% Constants
epsilon_0 = 8.854e-12; % Permittivity of free space
k = 1 / (4 * pi * epsilon_0); % Coulomb's constant

% Point Charges
Q1 = 8e-9; % Charge in Coulombs
Q2 = 8e-9;
P1 = [0, 1, 1]; % Position of first charge
P2 = [0, -1, 1]; % Position of second charge
P0 = [0, 0, 0]; % Field point

% Compute electric field due to Q1
r1 = P0 - P1;
r1_mag = norm(r1);
E1 = k * Q1 / r1_mag^2 * (r1 / r1_mag);

% Compute electric field due to Q2
r2 = P0 - P2;
r2_mag = norm(r2);
E2 = k * Q2 / r2_mag^2 * (r2 / r2_mag);

% Line Charge
rho_L = 4e-9; % Line charge density (C/m)
L1 = [7, 0, 0]; % Start of line
L2 = [0, 7, 0]; % End of line

% Numerical integration over the line charge
N = 1000; % Number of segments
t_vals = linspace(0, 1, N); % Parameter t from 0 to 1
E_line = [0, 0, 0];

for i = 1:N
    t = t_vals(i);
    L_point = (1 - t) * L1 + t * L2; % Parametric equation
    r = P0 - L_point;
    r_mag = norm(r);
    dL = norm(L2 - L1) / N;
    dE = k * rho_L * dL / r_mag^2 * (r / r_mag);
    E_line = E_line + dE;
end

% Total Electric Field
E_total = E1 + E2 + E_line;

% Display results
disp('Electric field at (0,0,0) (V/m):');
disp(E_total);
