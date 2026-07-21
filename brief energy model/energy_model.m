%% Organoid Energy Minimization: Parameter Change Study
% =====================================================================
% Modified: ECM pressure using neo-Hookean hyperelastic model
% Model Assumptions:
% 1. 3D spherical geometry with a central lumen.
% 2. Cells are arranged in a shell structure.
% 3. Energy contributions include: Apical/Lateral/Basal surface tensions,
%    Volume elasticity, ECM elastic resistance (hyperelastic), and Hydrostatic pressure.
% 4. Simplify the process using constant hydrostatic pressure
% ======================================================================

clear; 
clc; 
close all;

% =====================================================================
% 1. Global Parameter Configuration
% =====================================================================

% Fixed Biomechanical Parameters
fixed_params.P_cell = 200;          % Intracellular pressure (Pa)
fixed_params.miu = 0.4;             % Poisson's ratio (for ECM)
fixed_params.A = 5e-23;             % Non-adhesive modulus
fixed_params.d = 1e-6;              % Membrane thickness
fixed_params.V0 = 2*sqrt(3)*1e-15;  % Cell volume (m^3)
fixed_params.A_basal_0 = sqrt(3)/2 * (8e-6)^2; % Reference basal area
fixed_params.E_m = 2000;            % Membrane modulus (Pa)
fixed_params.k_1 = 0.01;            % Stiffness coeff
fixed_params.E_ecm = 1300;          % ECM Young's Modulus (Pa)
fixed_params.R0_ecm = 20e-6;        % Reference ECM radius (initial organoid size)

% Search Space Definition
search_space.N_vec = linspace(10, 500, 200);   % Cell number range
search_space.R_vec = linspace(20e-6, 200e-6, 3000); % Inner radius range

% =====================================================================
% 2. Simulation Cases
% =====================================================================

%% Case 1: Change Basal Adhesion (Gamma) - Variable N
fprintf('--- Running Case 1: Change Gamma (Variable N) ---\n');
params_c1 = fixed_params;
params_c1.lambda = 1e-3;       
params_c1.gamma_l = -1e-3;        % Changed alpha to gamma_l
change_var_c1 = struct('name', 'gamma', 'values', -1.5e-3:1e-4:1.5e-3);

[res1, val1] = compute_organoid_equilibrium(params_c1, search_space, change_var_c1, []);

figure('Name', 'Case 1: Gamma Change', 'Position', [100, 300, 1500, 500]);
plot_organoid_results(val1, res1, 'Basal Adhesion (\gamma) [J/m^2]', 'Case 1: \gamma Change (Var N)');


%% Case 2: Change Apical Tension (Lambda) - Variable N
fprintf('--- Running Case 2: Change Lambda (Variable N) ---\n');
params_c2 = fixed_params;
params_c2.gamma = -5e-4;       
params_c2.gamma_l = -1e-3;        % Changed alpha to gamma_l
change_var_c2 = struct('name', 'lambda', 'values', 1e-4:1e-4:3e-3);

[res2, val2] = compute_organoid_equilibrium(params_c2, search_space, change_var_c2, []);

figure('Name', 'Case 2: Lambda Change', 'Position', [100, 300, 1500, 500]);
plot_organoid_results(val2, res2, 'Apical Tension (\lambda) [N/m]', 'Case 2: \lambda Change (Var N)');


%% Case 3: Change Lateral Tension (Gamma_l) - Variable N
fprintf('--- Running Case 3: Change Gamma_l (Variable N) ---\n');
params_c3 = fixed_params;
params_c3.lambda = 1e-3;       
params_c3.gamma = -5e-4;       
change_var_c3 = struct('name', 'gamma_l', 'values', -2e-3:2e-4:2e-3);  % Changed alpha to gamma_l

[res3, val3] = compute_organoid_equilibrium(params_c3, search_space, change_var_c3, []);

figure('Name', 'Case 3: Gamma_l Change', 'Position', [100, 300, 1500, 500]);
plot_organoid_results(val3, res3, 'Lateral Tension (\gamma_l) [N/m]', 'Case 3: \gamma_l Change (Var N)');


% %% Case 4: FIXED N, Change Lateral Tension (Gamma_l)
% fprintf('--- Running Case 4: Fixed N=259, Change Gamma_l ---\n');
% N_fixed_val = 259;
% 
% params_c4 = fixed_params;
% params_c4.lambda = 1e-3;       
% params_c4.gamma = -5e-4;       
% change_var_c4 = struct('name', 'gamma_l', 'values', -5e-3:5e-4:5e-3);  % Changed alpha to gamma_l
% 
% [res4, val4] = compute_organoid_equilibrium(params_c4, search_space, change_var_c4, N_fixed_val);
% 
% figure('Name', 'Case 4: Fixed N Change', 'Position', [100, 300, 1500, 500]);
% plot_organoid_results(val4, res4, 'Lateral Tension (\gamma_l) [N/m]', ['Case 4: \gamma_l Change (Fixed N=' num2str(N_fixed_val) ')']);
% 
% 
% %% Case 5: FIXED N, Change Basal Adhesion (Gamma)
% fprintf('--- Running Case 5: Fixed N=259, Change Gamma ---\n');
% 
% params_c5 = fixed_params;
% params_c5.lambda = 1e-3;       
% params_c5.gamma_l = -1e-3;        % Changed alpha to gamma_l
% change_var_c5 = struct('name', 'gamma', 'values', -4e-3:2e-4:4e-3); 
% 
% [res5, val5] = compute_organoid_equilibrium(params_c5, search_space, change_var_c5, N_fixed_val);
% 
% figure('Name', 'Case 5: Fixed N Change', 'Position', [100, 300, 1500, 500]);
% plot_organoid_results(val5, res5, 'Basal Adhesion (\gamma) [J/m^2]', ['Case 5: \gamma Change (Fixed N=' num2str(N_fixed_val) ')']);
% 
% fprintf('All simulations completed.\n');
% 

% =====================================================================
% 3. Functions
% =====================================================================

function [results, change_values] = compute_organoid_equilibrium(fixed_params, search_space, change_var, fixed_N)
    % COMPUTE_ORGANOID_EQUILIBRIUM Calculates energy minima for parameter changes.
    % Modified to use neo-Hookean hyperelastic model for ECM pressure.
    %
    % Inputs:
    %   fixed_params : Struct of constant parameters.
    %   search_space : Struct with 'N_vec' and 'R_vec'.
    %   change_var   : Struct with 'name' and 'values' for the changing parameter.
    %   fixed_N      : Scalar to fix cell count, or [] for variable N.
    
    change_values = change_var.values;
    N_steps = length(change_values);
    
    % Pre-allocate results
    results.R_b = zeros(1, N_steps); 
    results.R_a = zeros(1, N_steps); 
    results.R_out = zeros(1, N_steps); 
    results.R_in = zeros(1, N_steps);  
    results.H = zeros(1, N_steps);     
    results.Curvature = zeros(1, N_steps); 
    results.Energy = zeros(1, N_steps);    
    results.N_cells = zeros(1, N_steps);   

    % Construct Grid based on whether N is fixed
    if isempty(fixed_N)
        [N_grid, R_grid] = meshgrid(search_space.N_vec, search_space.R_vec);
    else
        % Create a grid where N is constant
        N_vec_fixed = ones(size(search_space.R_vec)) * fixed_N;
        [N_grid, R_grid] = meshgrid(N_vec_fixed, search_space.R_vec);
    end
    [N_rows, N_cols] = size(N_grid);
    
    for k = 1:N_steps
        % Update changing parameter
        params = fixed_params;
        paramName = change_var.name;
        params.(paramName) = change_values(k);
        can = 0.3;
        % Extract params
        P_cell = params.P_cell; 
        miu = params.miu; 
        A = params.A; 
        d = params.d;
        V0 = params.V0; 
        A_basal_0 = params.A_basal_0; 
        E_m = params.E_m;
        lambda = params.lambda; 
        gamma = params.gamma; 
        gamma_l = params.gamma_l;  % Changed from alpha to gamma_l
        E_ecm = params.E_ecm; 
        R0_ecm = params.R0_ecm;  % Reference ECM radius
        
        % ECM shear modulus (for neo-Hookean)
        G_ecm = E_ecm / (2 * (1 + miu));

        F_energy = zeros(N_rows, N_cols);

        for j = 1:N_cols
            for i = 1:N_rows
                N_curr = N_grid(i,j);
                R_curr = R_grid(i,j);
                
                % Geometry Calculations
                r_out = ((3*N_curr*V0)/(4*pi) + R_curr^3)^(1/3);
                factor = sqrt(8*pi / (sqrt(3)*N_curr));
                x_curr = factor * r_out; % Basal edge
                y_curr = factor * R_curr; % Apical edge
                
                % Refined Height (Frustum constraint)
                Hy_curr = V0 / (x_curr^2 + x_curr*y_curr + y_curr^2);
                
                % Areas
                A_basal_curr = (sqrt(3)/2) * x_curr^2;
                A_apical_curr = (sqrt(3)/2) * y_curr^2;
                A_lateral_curr = sqrt(3) * (x_curr + y_curr) * Hy_curr;
                
                % =========================================================
                % MODIFIED: ECM Pressure using neo-Hookean hyperelastic model
                % =========================================================
                % Calculate stretch ratio (lambda = R/R0)
                stretch_ratio = r_out / R0_ecm;
                
                % Ensure physical validity (stretch ratio >= 1 for expansion)
                if stretch_ratio < 1
                    stretch_ratio = 1;  % No compression (organoid smaller than reference)
                end
                
                % % neo-Hookean pressure for spherical cavity expansion
                % % P = (G/2) * (5 - 4/lambda - 1/lambda^4)
                % % where lambda = R/R0 is the stretch ratio
                % 
                %     P_ecm = (G_ecm / 2) * (5 - 4/stretch_ratio - 1/stretch_ratio^4);
                % 
                % 
                % % Total pressure balance
                % % For energy calculation, we need the work done by ECM
                % T_press = (P_ecm^2) / (16 * G_ecm);

                  % Pressure term
                K_press(i,j) = (P_cell + can.*sqrt(abs(R_curr))).^2 / (24*G_ecm);
                K_press(i,i) = 1;
                
                % =========================================================
                % Total Free Energy
                % =========================================================
                E_val = lambda * A_apical_curr ...
                      + gamma * A_basal_curr ...
                      + gamma_l * A_lateral_curr ...  % Changed from alpha to gamma_l
                      + 2*A / (x_curr * y_curr) ...
                      + A / (Hy_curr^2) ...
                      + K_press(i,j)*r_out.*(x_curr.^6)./(y_curr.^4);  % ECM elastic energy
                
                % Physical Constraints (Penalty)
                if Hy_curr < 6e-6 || Hy_curr > 20e-6 || ...
                   x_curr < 6e-6 || x_curr > 20e-6 || ...
                   y_curr < 6e-6 || y_curr > 20e-6 || ...
                   stretch_ratio > 5  % Prevent unrealistically large expansions
                    E_val = NaN;
                end
                
                F_energy(i,j) = E_val;
            end
        end
        
        % Find Minimum
        [min_val, idx] = min(F_energy(:));
        [idx_r, idx_c] = ind2sub(size(F_energy), idx);
        
        if isnan(min_val)
            results.R_b(k) = NaN; 
            results.R_a(k) = NaN; 
            continue;
        end
        
        % Retrieve Optimal Geometry
        N_opt = N_grid(idx_r, idx_c);
        R_opt = R_grid(idx_r, idx_c);
        r_out_opt = ((3*N_opt*V0)/(4*pi) + R_opt^3)^(1/3);
        factor_opt = sqrt(8*pi / (sqrt(3)*N_opt));
        x_opt = factor_opt * r_out_opt;
        y_opt = factor_opt * R_opt;
        Hy_opt = V0 / (x_opt^2 + x_opt*y_opt + y_opt^2);
        
        % Store Results
        results.R_b(k) = x_opt;
        results.R_a(k) = y_opt;
        results.R_out(k) = r_out_opt;
        results.R_in(k) = R_opt;
        results.H(k) = Hy_opt;
        results.Curvature(k) = (x_opt - y_opt) / (y_opt * Hy_opt);
        results.Energy(k) = min_val;
        results.N_cells(k) = N_opt;
    end
end

function plot_organoid_results(change_vals, results, xlabel_str, title_prefix)
    % Standard plot with three subplots: Basal, Apical, and Outer Radius
    valid_idx = ~isnan(results.R_b);
    x_plot = change_vals(valid_idx);
    
    if isempty(x_plot), return; end

    subplot(1, 3, 1);
    scatter(x_plot, results.R_b(valid_idx), 40, 'b', 'filled');
    xlabel(xlabel_str, 'FontSize', 12); 
    ylabel('Basal Edge (R_b) [m]', 'FontSize', 12);
    title('Basal Geometry', 'FontSize', 12); 
    grid on;

    subplot(1, 3, 2);
    scatter(x_plot, results.R_a(valid_idx), 40, 'r', 'filled');
    xlabel(xlabel_str, 'FontSize', 12); 
    ylabel('Apical Edge (R_a) [m]', 'FontSize', 12);
    title('Apical Geometry', 'FontSize', 12); 
    grid on;

    subplot(1, 3, 3);
    scatter(x_plot, results.R_out(valid_idx), 40, 'k', 'filled');
    xlabel(xlabel_str, 'FontSize', 12); 
    ylabel('Outer Radius (R_{out}) [m]', 'FontSize', 12);
    title('Organoid Size', 'FontSize', 12); 
    grid on;

    sgtitle([title_prefix], 'FontSize', 16, 'FontWeight', 'bold');
    set(gcf, 'Color', 'w');
end