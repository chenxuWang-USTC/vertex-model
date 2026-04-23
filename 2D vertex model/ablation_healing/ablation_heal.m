%% SIMULATE_ANNULAR_TRAPEZOID_LANGEVIN_ABLATION_HEALING
% Simulates organoid ablation and healing dynamics using a Langevin approach.
% Features: Cell ablation, topological healing, and mechanical recovery.

clear; 
clc; 
close all;

%% ========== 1. Parameter Initialization ==========
n_initial = 40;                 % Initial number of cells
n = n_initial;                  % Current number of cells
dt = 0.25;                      % Time step
num_steps = 5000;               % Max simulation steps
mu = 0.4;                       % Damping coefficient / Viscosity

% Surface Tensions (N/m)
gamma_a = 2e-3;                 % Apical
gamma_b = -1e-3;                % Basal
gamma_l = -3e-3;                % Lateral

% Mechanical Parameters
A = 5e-23;                      % Gaussian non-adhesive modulus
K_area = 1e9;                   % Area constraint stiffness
kappa_bend = 1e-12;             % Bending stiffness
recover_step = 2000;            % Steps for full pressure recovery after healing

% Pack parameters into struct
params.gamma_a = gamma_a;
params.gamma_b = gamma_b;
params.gamma_l = gamma_l;
params.K_area = K_area;
params.A = A;
params.kappa_bend = kappa_bend;

% Visualization & Ablation Settings
vis_interval = 5;
enable_ablation = true;
ablation_step = 1;
ablated_cell_id = 10;
healing_enabled = true;
healing_threshold = 1.0e-7;     % Distance threshold for healing
is_healed_permanently = false;
healed_step = 0;
is_ablated = false(n_initial, 1);

%% ========== 2. Initial Geometry ==========
L = intital(); % Ensure 'intital.m' exists in path
r_a = L.r_a;
r_b = L.r_b;
A0 = L.A0;
params.A0 = A0;

% Symmetrize initial positions to reduce numerical noise
fprintf('=== Symmetrizing Initial Geometry ===\n');
R_a_mean = mean(sqrt(sum(r_a.^2, 2)));
R_b_mean = mean(sqrt(sum(r_b.^2, 2)));
n = length(r_a);
theta = linspace(0, 2*pi, n+1)'; 
theta(end) = [];
r_a = R_a_mean * [cos(theta), sin(theta)];
r_b = R_b_mean * [cos(theta), sin(theta)];
x = [r_a(:), r_b(:)];


%% ========== 3. Visualization Setup ==========
figure('Position', [50, 50, 1500, 750], 'Name', 'Organoid Ablation & Healing');

% --- Left Panel: Geometry ---
subplot(1, 2, 1); 
hold on; 
grid on; 
axis equal;
h_inner = plot([r_a(:,1); r_a(1,1)], [r_a(:,2); r_a(1,2)], 'ro-', 'LineWidth', 2);
h_outer = plot([r_b(:,1); r_b(1,1)], [r_b(:,2); r_b(1,2)], 'bo-', 'LineWidth', 2);
h_side = zeros(n, 1);
for i = 1:n
    h_side(i) = plot([r_a(i,1), r_b(i,1)], [r_a(i,2), r_b(i,2)], 'g--', 'LineWidth', 1.2);
end

% Force vectors (initialized to zero)
% g_force_a = quiver(r_a(:,1), r_a(:,2), zeros(n,1), zeros(n,1), 0, 'y', 'AutoScale', 'off');
% g_force_b = quiver(r_b(:,1), r_b(:,2), zeros(n,1), zeros(n,1), 0, 'g', 'AutoScale', 'off');

title('Organoid Ablation and Healing', 'FontSize', 24, 'FontWeight', 'bold');
xlim([-1.5e-4, 1.5e-4]); 
ylim([-1.5e-4, 1.5e-4]);
xlabel('X (m)'); 
ylabel('Y (m)');
set(gca, 'FontSize', 14);
status_text = text(0.05, 0.90, 'Status: Initial', 'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold');

% --- Right Panel: Radius Evolution ---
subplot(1, 2, 2); 
hold on; 
grid on;
h_radius_avg = plot(NaN, NaN, 'b-', 'LineWidth', 2, 'DisplayName', 'Avg Radius');
xlabel('Step'); 
ylabel('Radius (m)');
title('Radius Evolution', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best');
xlim([1, num_steps]);
ylim([5e-5, 8e-5]);
set(gca, 'FontSize', 14);

% Logging Arrays
radius_log = zeros(num_steps, 1);

%% ========== 4. Video Initialization ==========
videoFileName = 'ablation_healing_fixed.mp4';
if exist(videoFileName, 'file') %delete before results
    delete(videoFileName);
end
video = VideoWriter(videoFileName, 'MPEG-4');
video.FrameRate = 60;
open(video);

%% ========== 5. Main Simulation Loop ==========
for step = 1:num_steps
    
    % --- 5.1 Ablation Trigger ---
    if enable_ablation && step == ablation_step && ~is_healed_permanently
        if ablated_cell_id <= n
            is_ablated(ablated_cell_id) = true;
            fprintf('*** [Step %d] Cell %d Ablated ***\n', step, ablated_cell_id);
            set(status_text, 'String', 'Status: Ablated', 'Color', 'r');
        end
    end

    % --- 5.2 Healing Detection & Topology Update ---
    if healing_enabled && any(is_ablated) && ~is_healed_permanently
        idx = find(is_ablated, 1);
        next_idx = mod(idx, n) + 1;
        
        dist_inner = norm(r_a(next_idx, :) - r_a(idx, :));
        dist_outer = norm(r_b(next_idx, :) - r_b(idx, :));
        
        % Heal if either inner or outer vertices converge
        if dist_inner < healing_threshold || dist_outer < healing_threshold
            fprintf(' [Step %d] Healing Triggered\n', step);
            
            % Remove ablated cell from arrays
            if idx == n
                r_a = r_a(1:end-1, :); 
                r_b = r_b(1:end-1, :);
            else
                r_a = [r_a(1:idx-1, :); r_a(idx+1:end, :)];
                r_b = [r_b(1:idx-1, :); r_b(idx+1:end, :)];
            end
            
            n = n - 1;
            is_healed_permanently = true;
            healed_step = step;
            is_ablated = false(n, 1); % Reset ablation mask
            
            % Rebuild Visualization Handles for new N
            subplot(1, 2, 1);
            delete(h_side); 
            % delete(g_force_a); 
            % delete(g_force_b);
            h_side = zeros(n, 1);
            for i = 1:n
                h_side(i) = plot([r_a(i,1), r_b(i,1)], [r_a(i,2), r_b(i,2)], 'g--', 'LineWidth', 1.2);
            end
            g_force_a = quiver(r_a(:,1), r_a(:,2), zeros(n,1), zeros(n,1), 0, 'y', 'AutoScale', 'off');
            g_force_b = quiver(r_b(:,1), r_b(:,2), zeros(n,1), zeros(n,1), 0, 'g', 'AutoScale', 'off');
            
            set(status_text, 'String', 'Status: Healed', 'Color', 'b');
            drawnow;
        end
    end

    %% --- 5.3 Force Calculations ---
    P_ecm = zeros(n,1); P_water = zeros(n,1); A_eff = zeros(n,1);
    v_a = zeros(n,2); v_b = zeros(n,2); v_side = zeros(n,2);
    
    for i = 1:n
        next = mod(i, n) + 1;
        v_a(i,:) = r_a(next,:) - r_a(i,:);
        v_b(i,:) = r_b(next,:) - r_b(i,:);
        v_side(i,:) = r_b(i,:) - r_a(i,:);
        
        % Pressure Models
        P_ecm(i) = 700/(77e-6)*norm(r_b(i,:)) - 400;
        P_water_raw = 1000*sqrt(66e-6)/sqrt(norm(r_a(i,:))) - 700;
        
        % Water Pressure Recovery Logic
        if ~is_healed_permanently
            P_water(i) = 0;
        else
            progress = min(1, (step - healed_step) / recover_step);
            P_water(i) = P_water_raw * progress;
        end
        
        % Effective Area Modulus (Zero if ablated)
        A_eff(i) = A * (~is_ablated(i));
    end
    
    % Geometric Quantities
    l_a = sqrt(sum(v_a.^2, 2)); l_b = sqrt(sum(v_b.^2, 2)); l_side = sqrt(sum(v_side.^2, 2));
    e_a = v_a ./ l_a; e_b = v_b ./ l_b; e_side = v_side ./ l_side;

    %% --- 5.4 Bending Forces ---
    F_bend_a = zeros(n, 2); F_bend_b = zeros(n, 2); E_bend = 0;
    target_cos = 0; bend_coeff = 2 * kappa_bend;
    
    for i = 1:n
        if is_ablated(i)
            continue;
        end
        next_idx = mod(i, n) + 1;
        
        % Vectors for angles
        u_in = r_a(next_idx,:) - r_a(i,:); u_out = r_b(next_idx,:) - r_b(i,:);
        w_in = r_b(i,:) - r_a(i,:); w_next = r_b(next_idx,:) - r_a(next_idx,:);
        
        len_uin = norm(u_in); len_uout = norm(u_out);
        len_win = norm(w_in); len_wnext = norm(w_next);
        
        if min([len_uin, len_uout, len_win, len_wnext]) < 1e-12
            continue; 
        end
        
        e_uin = u_in/len_uin; e_uout = u_out/len_uout;
        e_win = w_in/len_win; e_wnext = w_next/len_wnext;
        
        % Calculate 4 corner angles and gradients (Simplified for brevity)
        % Note: Full gradient derivation retained from original logic
        % Angle 1: Inner Start
        cos_t1 = max(-1, min(1, dot(e_win, e_uin)));
        E_bend = E_bend + 0.5*bend_coeff*(cos_t1 - target_cos)^2;
        if abs(cos_t1 - target_cos) > 1e-8
            f = -bend_coeff*(cos_t1 - target_cos);
            g_win = (e_uin - cos_t1*e_win)/len_win; g_uin = (e_win - cos_t1*e_uin)/len_uin;
            F_bend_b(i,:) = F_bend_b(i,:) + f*g_win; F_bend_a(i,:) = F_bend_a(i,:) - f*g_win;
            F_bend_a(next_idx,:) = F_bend_a(next_idx,:) + f*g_uin; F_bend_a(i,:) = F_bend_a(i,:) - f*g_uin;
        end
        % (Other 3 angles follow similar pattern - kept concise here)
        % Angle 2: Outer Start
        cos_t2 = max(-1, min(1, dot(-e_win, e_uout)));
        E_bend = E_bend + 0.5*bend_coeff*(cos_t2 - target_cos)^2;
        if abs(cos_t2 - target_cos) > 1e-8
            f = -bend_coeff*(cos_t2 - target_cos);
            g_win = (e_uout - cos_t2*(-e_win))/len_win; g_uout = ((-e_win) - cos_t2*e_uout)/len_uout;
            F_bend_a(i,:) = F_bend_a(i,:) + f*g_win; F_bend_b(i,:) = F_bend_b(i,:) - f*g_win;
            F_bend_b(next_idx,:) = F_bend_b(next_idx,:) + f*g_uout; F_bend_b(i,:) = F_bend_b(i,:) - f*g_uout;
        end
        % Angle 3: Inner End
        cos_t3 = max(-1, min(1, dot(-e_uin, e_wnext)));
        E_bend = E_bend + 0.5*bend_coeff*(cos_t3 - target_cos)^2;
        if abs(cos_t3 - target_cos) > 1e-8
            f = -bend_coeff*(cos_t3 - target_cos);
            g_uin = (e_wnext - cos_t3*(-e_uin))/len_uin; g_wnext = ((-e_uin) - cos_t3*e_wnext)/len_wnext;
            F_bend_a(i,:) = F_bend_a(i,:) + f*g_uin; F_bend_a(next_idx,:) = F_bend_a(next_idx,:) - f*g_uin;
            F_bend_b(next_idx,:) = F_bend_b(next_idx,:) + f*g_wnext; F_bend_a(next_idx,:) = F_bend_a(next_idx,:) - f*g_wnext;
        end
        % Angle 4: Outer End
        cos_t4 = max(-1, min(1, dot(-e_uout, e_wnext)));
        E_bend = E_bend + 0.5*bend_coeff*(cos_t4 - target_cos)^2;
        if abs(cos_t4 - target_cos) > 1e-8
            f = -bend_coeff*(cos_t4 - target_cos);
            g_uout = (e_wnext - cos_t4*(-e_uout))/len_uout; g_wnext = ((-e_uout) - cos_t4*e_wnext)/len_wnext;
            F_bend_b(i,:) = F_bend_b(i,:) + f*g_uout; F_bend_b(next_idx,:) = F_bend_b(next_idx,:) - f*g_uout;
            F_bend_b(next_idx,:) = F_bend_b(next_idx,:) + f*g_wnext; F_bend_a(next_idx,:) = F_bend_a(next_idx,:) - f*g_wnext;
        end
    end

    %% --- 5.5 Surface Tension Forces ---
    grad_U_a = zeros(n, 2); grad_U_b = zeros(n, 2);
    grad_side_a = zeros(n, 2); grad_side_b = zeros(n, 2);
    
    for i = 1:n
        next = mod(i, n) + 1;
        % Effective tensions (Zero if ablated)
        g_a_eff = gamma_a * (~is_ablated(i));
        g_b_eff = gamma_b * (~is_ablated(i));
        g_l_eff = gamma_l * (~is_ablated(i));

        % Apical/Basal Gradients
        grad_U_a(i,:) = grad_U_a(i,:) + sqrt(3)*g_a_eff*l_a(i)*(-e_a(i,:));
        grad_U_a(next,:) = grad_U_a(next,:) + sqrt(3)*g_a_eff*l_a(i)*e_a(i,:);
        grad_U_b(i,:) = grad_U_b(i,:) + sqrt(3)*g_b_eff*l_b(i)*(-e_b(i,:));
        grad_U_b(next,:) = grad_U_b(next,:) + sqrt(3)*g_b_eff*l_b(i)*e_b(i,:);

        % Lateral Gradients
        grad_side_a(i,:) = grad_side_a(i,:) - e_a(i,:)*g_l_eff/sqrt(12)*(l_side(i)+l_side(next)) ...
                         - e_side(i,:)*g_l_eff/sqrt(12)*(l_a(i)+l_b(i));
        grad_side_a(next,:) = grad_side_a(next,:) + e_a(i,:)*g_l_eff/sqrt(12)*(l_side(i)+l_side(next)) ...
                            - e_side(next,:)*g_l_eff/sqrt(12)*(l_a(i)+l_b(i));
        grad_side_b(i,:) = grad_side_b(i,:) - e_b(i,:)*g_l_eff/sqrt(12)*(l_side(i)+l_side(next)) ...
                         + e_side(i,:)*g_l_eff/sqrt(12)*(l_a(i)+l_b(i));
        grad_side_b(next,:) = grad_side_b(next,:) + e_b(i,:)*g_l_eff/sqrt(12)*(l_side(i)+l_side(next)) ...
                            + e_side(next,:)*g_l_eff/sqrt(12)*(l_a(i)+l_b(i));
    end
    F_surface_a = -grad_U_a; F_surface_b = -grad_U_b;

    %% --- 5.6 Pressure Forces ---
    F_water = zeros(n, 2); F_ecm = zeros(n, 2);
    for i = 1:n
        prev = mod(i-2, n) + 1;
        % Inner Pressure Force
        v_prev = v_a(prev,:); v_i = v_a(i,:);
        t1 = v_i/norm(v_i); t2 = v_prev/norm(v_prev);
        n1 = [t1(2), -t1(1)]; e1 = n1/norm(n1);
        n2 = [t2(2), -t2(1)]; e2 = n2/norm(n2);
        F_water(i,:) = P_water(i)*sqrt(3)/2 * (l_a(prev)^2*e2 + l_a(i)^2*e1);
        
        % Outer Pressure Force
        v_prev = v_b(prev,:); v_i = v_b(i,:);
        t1 = v_i/norm(v_i); t2 = v_prev/norm(v_prev);
        n1 = [-t1(2), t1(1)]; e1 = n1/norm(n1);
        n2 = [-t2(2), t2(1)]; e2 = n2/norm(n2);
        F_ecm(i,:) = P_ecm(i)*sqrt(3)/2 * (l_b(prev)^2*e2 + l_b(i)^2*e1);
    end

    %% --- 5.7 Area Constraint Forces ---
    F_area_a = zeros(n, 2); F_area_b = zeros(n, 2);
    for i = 1:n
        next = mod(i, n) + 1;
        P1=r_a(i,:); P2=r_b(i,:); P3=r_b(next,:); P4=r_a(next,:);
        A_curr = quad_area(P1, P2, P3, P4);
        dA = A_curr - A0;
        
        % Gradients of Area wrt vertices
        dA_dP1 = -0.5 * [-(P2(2)-P4(2)), P2(1)-P4(1)];
        dA_dP2 = -0.5 * [-(P3(2)-P1(2)), P3(1)-P1(1)];
        dA_dP3 = -0.5 * [-(P4(2)-P2(2)), P4(1)-P2(1)];
        dA_dP4 = -0.5 * [-(P1(2)-P3(2)), P1(1)-P3(1)];
        
        K_eff = K_area * (~is_ablated(i));
        F_area_a(i,:) = F_area_a(i,:) - K_eff*dA*dA_dP1;
        F_area_b(i,:) = F_area_b(i,:) - K_eff*dA*dA_dP2;
        F_area_b(next,:) = F_area_b(next,:) - K_eff*dA*dA_dP3;
        F_area_a(next,:) = F_area_a(next,:) - K_eff*dA*dA_dP4;
    end

    %% --- 5.8 Gaussian Curvature Forces ---
    F_gauss_a = zeros(n,2); F_gauss_b = zeros(n,2);
    for i = 1:n
        next = mod(i, n) + 1; prev = mod(i-2, n) + 1;
        % Simplified Gaussian term calculation
        term_a = -(A_eff(i)*2*(-e_a(i,:))/(-l_b(i)*l_a(i)^2) + A_eff(i)*(-e_side(i,:))/(-l_side(i)^2*l_side(next))...
            + A_eff(prev)*2*(e_a(prev,:))/(-l_b(prev)*l_a(prev)^2) + A_eff(prev)*(-e_side(i,:))/(-l_side(i)^2*l_side(prev)));
        term_b = -(A_eff(i)*2*(-e_b(i,:))/(-l_a(i)*l_b(i)^2) + A_eff(i)*(e_side(i,:))/(-l_side(i)^2*l_side(next))...
            + A_eff(prev)*2*(e_b(prev,:))/(-l_a(prev)*l_b(prev)^2) + A_eff(prev)*(e_side(i,:))/(-l_side(i)^2*l_side(prev)));
        
        F_gauss_a(i,:) = term_a;
        F_gauss_b(i,:) = term_b;
    end

    %% --- 5.9 Total Force & Dynamics ---
    vel_a = zeros(n,2); vel_b = zeros(n,2);
    for i = 1:n
        F_tot_a = F_area_a(i,:) + F_surface_a(i,:) - grad_side_a(i,:) + F_water(i,:) + F_bend_a(i,:) + F_gauss_a(i,:);
        F_tot_b = F_area_b(i,:) + F_surface_b(i,:) - grad_side_b(i,:) + F_ecm(i,:) + F_bend_b(i,:) + F_gauss_b(i,:);
        
        vel_a(i,:) = F_tot_a / mu;
        vel_b(i,:) = F_tot_b / mu;
    end

    % Update Positions
    r_a = r_a + vel_a * dt;
    r_b = r_b + vel_b * dt;

    %% --- 5.10 Convergence Check ---
    max_vel = max([max(sqrt(sum(vel_a.^2, 2))); max(sqrt(sum(vel_b.^2, 2)))]);
    if max_vel < 5e-10 && step > 100
        fprintf('Converged at step %d\n', step);
        break;
    end

    %% --- 5.11 Data Logging ---
    R_a_curr = mean(sqrt(sum(r_a.^2, 2)));
    R_b_curr = mean(sqrt(sum(r_b.^2, 2)));
    radius_log(step) = (R_a_curr + R_b_curr) / 2;

    %% --- 5.12 Visualization Update ---
    if mod(step, vis_interval) == 0 || step == 1
        % Update Geometry Plots
        if ~is_healed_permanently && any(is_ablated)
            idx = find(is_ablated, 1);
            % Plot with gap for ablated cell
            xa = [r_a(1:idx,1); NaN; r_a(idx+1:end,1); r_a(1,1)];
            ya = [r_a(1:idx,2); NaN; r_a(idx+1:end,2); r_a(1,2)];
            set(h_inner, 'XData', xa, 'YData', ya);
            xb = [r_b(1:idx,1); NaN; r_b(idx+1:end,1); r_b(1,1)];
            yb = [r_b(1:idx,2); NaN; r_b(idx+1:end,2); r_b(1,2)];
            set(h_outer, 'XData', xb, 'YData', yb);
        else
            set(h_inner, 'XData', [r_a(:,1); r_a(1,1)], 'YData', [r_a(:,2); r_a(1,2)]);
            set(h_outer, 'XData', [r_b(:,1); r_b(1,1)], 'YData', [r_b(:,2); r_b(1,2)]);
        end
        
        % Update Side Lines
        for i = 1:min(length(h_side), n)
            set(h_side(i), 'XData', [r_a(i,1), r_b(i,1)], 'YData', [r_a(i,2), r_b(i,2)]);
            if is_ablated(i) || is_ablated(mod(i-2,n)+1)
                set(h_side(i), 'Color', 'k', 'LineWidth', 2.5);
            else
                set(h_side(i), 'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.2);
            end
        end
        
        % Update Radius Plot
        subplot(1, 2, 2);
        set(h_radius_avg, 'XData', 1:step, 'YData', radius_log(1:step));
        
        % Update Force Vectors
        % set(g_force_a, 'UData', F_gauss_a(:,1)*1000, 'VData', F_gauss_a(:,2)*1000);
        % set(g_force_b, 'UData', F_gauss_b(:,1)*1000, 'VData', F_gauss_b(:,2)*1000);
        
        drawnow limitrate;
        writeVideo(video, getframe(gcf));
    end
end

close(video);
fprintf('Simulation Finished. Video saved.\n');