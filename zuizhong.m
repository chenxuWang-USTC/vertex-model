%% Organoid Energy Minimization: Parameter Change Study
% =====================================================================
% using a biref energy model,
% We can quickly obtain the approximate radius corresponding to the number of organoid cells, 
% as well as the general growth trend of the organoid.
% Model Assumptions:
% 1. 3D spherical geometry with a central lumen.
% 2. Cells are arranged in a shell structure.
% 3. Energy contributions include: Apical/Lateral/Basal surface tensions,
%    Volume elasticity, ECM elastic resistance, and Hydrostatic pressure.
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
fixed_params.miu = 0.4;             % Poisson's ratio
fixed_params.A = 5e-23;             % Non-adhesive modulus
fixed_params.d = 1e-6;              % Membrane thickness
fixed_params.V0 = 2*sqrt(3)*1e-15;  % Cell volume (m^3)
fixed_params.A_basal_0 = sqrt(3)/2 * (8e-6)^2; % Reference basal area
fixed_params.E_m = 2000;            % Membrane modulus (Pa)
fixed_params.k_1 = 0.01;            % Stiffness coeff
fixed_params.E_ecm = 1300;          % ECM Modulus (Pa)

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
params_c1.alpha = 1e-3;        
change_var_c1 = struct('name', 'gamma', 'values', -1.5e-3:1e-4:1.5e-3);% Adjustable range to obtain a single parameter

[res1, val1] = compute_organoid_equilibrium(params_c1, search_space, change_var_c1, []);

figure('Name', 'Case 1: Gamma Change', 'Position', [100, 300, 1500, 500]);
plot_organoid_results(val1, res1, 'Basal Adhesion (\gamma) [J/m^2]', 'Case 1: \gamma Change (Var N)');


%% Case 2: Change Apical Tension (Lambda) - Variable N
fprintf('--- Running Case 2: Change Lambda (Variable N) ---\n');
params_c2 = fixed_params;
params_c2.gamma = -5e-4;       
params_c2.alpha = 1e-3;        
change_var_c2 = struct('name', 'lambda', 'values', 1e-4:1e-4:3e-3);

[res2, val2] = compute_organoid_equilibrium(params_c2, search_space, change_var_c2, []);

figure('Name', 'Case 2: Lambda Change', 'Position', [100, 300, 1500, 500]);
plot_organoid_results(val2, res2, 'Apical Tension (\lambda) [N/m]', 'Case 2: \lambda Change (Var N)');


%% Case 3: Change Lateral Tension (Alpha) - Variable N
fprintf('--- Running Case 3: Change Alpha (Variable N) ---\n');
params_c3 = fixed_params;
params_c3.lambda = 1e-3;       
params_c3.gamma = -5e-4;       
% Narrowed range to avoid physical constraints violations
change_var_c3 = struct('name', 'alpha', 'values', -2e-3:2e-4:2e-3); 

[res3, val3] = compute_organoid_equilibrium(params_c3, search_space, change_var_c3, []);

figure('Name', 'Case 3: Alpha Change', 'Position', [100, 300, 1500, 500]);
plot_organoid_results(val3, res3, 'Lateral Tension (\alpha) [N/m]', 'Case 3: \alpha Change (Var N)');


%% Case 4: FIXED N, Change Lateral Tension (Alpha)
fprintf('--- Running Case 4: Fixed N=259, Change Alpha ---\n');
N_fixed_val = 259; % Define the fixed cell number

params_c4 = fixed_params;
params_c4.lambda = 1e-3;       
params_c4.gamma = -5e-4;       
change_var_c4 = struct('name', 'alpha', 'values', -5e-3:5e-4:5e-3); 

% Pass N_fixed_val to the function
[res4, val4] = compute_organoid_equilibrium(params_c4, search_space, change_var_c4, N_fixed_val);

figure('Name', 'Case 4: Fixed N Change', 'Position', [100, 300, 1500, 500]);
plot_organoid_results(val4, res4, 'Lateral Tension (\alpha) [N/m]', ['Case 4: \alpha Change (Fixed N=' num2str(N_fixed_val) ')']);


%% Case 5: FIXED N, Change Basal Adhesion (Gamma)
fprintf('--- Running Case 5: Fixed N=259, Change Gamma ---\n');

params_c5 = fixed_params;
params_c5.lambda = 1e-3;       
params_c5.alpha = -1e-3;        
change_var_c5 = struct('name', 'gamma', 'values', -4e-3:2e-4:4e-3); 

[res5, val5] = compute_organoid_equilibrium(params_c5, search_space, change_var_c5, N_fixed_val);

figure('Name', 'Case 5: Fixed N Change', 'Position', [100, 300, 1500, 500]);
plot_organoid_results(val5, res5, 'Basal Adhesion (\gamma) [J/m^2]', ['Case 5: \gamma Change (Fixed N=' num2str(N_fixed_val) ')']);

fprintf('All simulations completed.\n');


% =====================================================================
% 3. Functions
% =====================================================================

function [results, change_values] = compute_organoid_equilibrium(fixed_params, search_space, change_var, fixed_N)
    % COMPUTE_ORGANOID_EQUILIBRIUM Calculates energy minima for parameter changes.
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
        
        % Extract params
        P_cell = params.P_cell; miu = params.miu; A = params.A; d = params.d;
        V0 = params.V0; A_basal_0 = params.A_basal_0; E_m = params.E_m;
        lambda = params.lambda; gamma = params.gamma; alpha = params.alpha;
        E_ecm = params.E_ecm; 
        G = E_ecm / (2*(1+miu));

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
                
                % ECM Pressure & Energy Term
                P_ecm = P_cell - (2*d*E_m*(A_basal_curr/A_basal_0 - 1)/r_out);
                T_press = (P_ecm^2) / (16*G);
                
                % Total Free Energy
                E_val = lambda * A_apical_curr ...
                      + gamma * A_basal_curr ...
                      - alpha * A_lateral_curr ...
                      + 2*A / (x_curr * y_curr) ...
                      + A / (Hy_curr^2) ...
                      + T_press * r_out * A_basal_curr;
                
                % Physical Constraints (Penalty)
                if Hy_curr < 6e-6 || Hy_curr > 20e-6 || ...
                   x_curr < 6e-6 || x_curr > 20e-6 || ...
                   y_curr < 6e-6 || y_curr > 20e-6
                    E_val = NaN;
                end
                
                F_energy(i,j) = E_val;
            end
        end
        
        % Find Minimum
        [min_val, idx] = min(F_energy(:));
        [idx_r, idx_c] = ind2sub(size(F_energy), idx);
        
        if isnan(min_val)
            results.R_b(k) = NaN; results.R_a(k) = NaN; 
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
    % Standard plot for Variable-N cases
    valid_idx = ~isnan(results.R_b);
    x_plot = change_vals(valid_idx);
    
    if isempty(x_plot), return; end

    subplot(1, 3, 1);
    scatter(x_plot, results.R_b(valid_idx), 40, 'b', 'filled');
    xlabel(xlabel_str, 'FontSize', 12); ylabel('Basal Edge (R_b) [m]', 'FontSize', 12);
    title('Basal Geometry', 'FontSize', 12); grid on;

    subplot(1, 3, 2);
    scatter(x_plot, results.R_a(valid_idx), 40, 'r', 'filled');
    xlabel(xlabel_str, 'FontSize', 12); ylabel('Apical Edge (R_a) [m]', 'FontSize', 12);
    title('Apical Geometry', 'FontSize', 12); grid on;

    subplot(1, 3, 3);
    scatter(x_plot, results.R_out(valid_idx), 40, 'k', 'filled');
    xlabel(xlabel_str, 'FontSize', 12); ylabel('Outer Radius (R_{out}) [m]', 'FontSize', 12);
    title('Organoid Size', 'FontSize', 12); grid on;

    sgtitle([title_prefix], 'FontSize', 16, 'FontWeight', 'bold');
    set(gcf, 'Color', 'w');
end

% function plot_fixed_n_results(change_vals, results, xlabel_str, title_prefix)
%     % Specialized plot for Fixed-N cases
%     valid_idx = ~isnan(results.R_in);
%     x_plot = change_vals(valid_idx);
% 
%     if isempty(x_plot)
%         warning('No valid data for %s', title_prefix);
%         return;
%     end
% 
%     % Plot 1: Inner vs Outer Radius
%     subplot(1, 3, 1);
%     plot(x_plot, results.R_in(valid_idx)*1e6, 'b-', 'LineWidth', 2); hold on;
%     plot(x_plot, results.R_out(valid_idx)*1e6, 'r-', 'LineWidth', 2);
%     xlabel(xlabel_str, 'FontSize', 12); ylabel('Radius (\mu m)', 'FontSize', 12);
%     legend('Inner (Lumen)', 'Outer', 'Location', 'best');
%     title('Radii Evolution', 'FontSize', 12); grid on;
% 
%     % Plot 2: Cell Height
%     subplot(1, 3, 2);
%     scatter(x_plot, results.H(valid_idx)*1e6, 40, 'g', 'filled');
%     xlabel(xlabel_str, 'FontSize', 12); ylabel('Cell Height (\mu m)', 'FontSize', 12);
%     title('Cell Shape Response', 'FontSize', 12); grid on;
% 
%     % Plot 3: Energy
%     subplot(1, 3, 3);
%     plot(x_plot, results.Energy(valid_idx), 'k-', 'LineWidth', 2);
%     xlabel(xlabel_str, 'FontSize', 12); ylabel('Min Energy (J)', 'FontSize', 12);
%     title('System Stability', 'FontSize', 12); grid on;
% 
%     sgtitle([title_prefix], 'FontSize', 16, 'FontWeight', 'bold');
%     set(gcf, 'Color', 'w');
% end