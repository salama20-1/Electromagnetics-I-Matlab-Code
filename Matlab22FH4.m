V = 0;
S1 = 0; 
S2 = 0; 
S3 = 0;
S4 = 0; 
S5 = 0; 
S6 = 0; 
phi = (pi/4); % initializing phi to lower boundary
theta = (pi/4); % initializing theta to lower boundary
num_of_r_steps = 500; % initializing r discretization
num_of_phi_steps = 500; % initializing phi discretization
num_of_theta_steps = 500; % initializing theta discretization
dr = (2-0)/num_of_r_steps; % r increment
dphi = ((pi/2)-(pi/4))/num_of_phi_steps; % phi increment
dtheta = ((pi/2)-(pi/4))/num_of_theta_steps; % theta increment

% calculating volume
for k = 1:num_of_r_steps
    r = k * dr; % Properly update r inside loop
    for j = 1:num_of_theta_steps
        theta = (pi/4) + (j-1) * dtheta; % Explicit update
        for i = 1:num_of_phi_steps
            phi = (pi/4) + (i-1) * dphi; % Explicit update
            V = V + (r^2) * sin(phi) * dr * dtheta * dphi;
        end
    end
end
fprintf('Volume= %.4f\n', V); % Display volume

% calculating S1 and S2
r1 = 0;
r2 = 2;
S1 = 0;
S2 = 0;
for k = 1:num_of_phi_steps
    phi = (pi/4) + (k-1) * dphi; % Explicit phi update
    for i = 1:num_of_theta_steps
        theta = (pi/4) + (i-1) * dtheta; % Explicit theta update
        S1 = S1 + (r1^2) * sin(phi) * dphi * dtheta;
        S2 = S2 + (r2^2) * sin(phi) * dphi * dtheta;
    end
end

% calculating S3 and S4
S3 = 0;
S4 = 0;
for j = 1:num_of_r_steps
    r = j * dr; % Ensure r is correctly updated
    for i = 1:num_of_theta_steps
        theta = (pi/4) + (i-1) * dtheta; % Explicit update
        S3 = S3 + r * sin(pi/2) * dtheta * dr; % phi1 = pi/2
        S4 = S4 + r * sin(pi/4) * dtheta * dr; % phi2 = pi/4
    end
end

% calculating S5 and S6
S5 = 0;
for k = 1:num_of_phi_steps
    phi = (pi/4) + (k-1) * dphi; % Explicit phi update
    for j = 1:num_of_r_steps
        r = j * dr; % Proper r update
        S5 = S5 + r * dphi * dr;
    end
end
S6 = S5; % Symmetry

S = S1 + S2 + S3 + S4 + S5 + S6;
fprintf('Surface area = %.4f\n', S);
