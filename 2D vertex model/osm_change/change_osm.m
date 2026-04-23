clear; clc;
% SIMULATE_ANNULAR_TRAPEZOID_LANGEVIN
% Simulates annular trapezoidal structure evolution using viscous dynamics:
% mu * dr/dt = -grad(U) + F_water + F_ecm + F_area

%% ========== 1. Parameter Initialization ==========
n = 40;                    % Number of cells
dt = 0.25;                 % Time step
num_steps = 10000;         % Total steps
mu = 0.4;                  % Viscous damping

% Surface Tensions (N/m) ,we can change these tension 
gamma_a = 2e-3;            % Apical (Inner)
gamma_b = -1e-3;           % Basal (Outer)
gamma_l = -3e-3;           % Lateral (Side)

% Mechanical Parameters
A_gauss = 5e-23;           % Gaussian curvature parameter
K_area = 1e9;              % Area constraint stiffness
kappa_bend = 1e-12;        % Bending stiffness

% Pack parameters
params.gamma_a = gamma_a; 
params.gamma_b = gamma_b; 
params.gamma_l = gamma_l;
params.K_area = K_area; 
params.A = A_gauss; 
params.kappa_bend = kappa_bend;

vis_interval = 10;       % Visualization update frequency
W_water = 0; 
W_ecm = 0;  % Work terms

%% ========== 2. Initial Geometry & Symmetrization ==========
L = intital();             % Load initial geometry from helper function
r_a = L.r_a; 
r_b = L.r_b; 
A0 = L.A0; 
params.A0 = A0;

fprintf('=== Symmetrizing Initial Geometry ===\n');
% Calculate mean radii and enforce uniform angular distribution
R_a_initial = mean(sqrt(sum(r_a.^2, 2)));
R_b_initial = mean(sqrt(sum(r_b.^2, 2)));
R_a_mean = mean(R_a_initial);
R_b_mean = mean(R_b_initial);
n = length(r_a);
theta = linspace(0, 2*pi, n+1)'; 
theta(end) = [];  
r_a = R_a_mean * [cos(theta), sin(theta)];
r_b = R_b_mean * [cos(theta), sin(theta)];

% Calculate initial total area for normalization
total_A0 = polyarea(r_b(:,1), r_b(:,2));

%% ========== 3. Pre-compute Geometric Quantities ==========
v_a = zeros(n,2); 
v_b = zeros(n,2); 
v_side = zeros(n,2);
l_a = zeros(n,1); 
l_b = zeros(n,1); 
l_side = zeros(n,1);
e_a = zeros(n,2); 
e_b = zeros(n,2); 
e_side = zeros(n,2);

for i = 1:n
    next_idx = mod(i, n) + 1;
    v_a(i,:) = r_a(next_idx,:) - r_a(i,:);
    v_b(i,:) = r_b(next_idx,:) - r_b(i,:);
    v_side(i,:) = r_b(i,:) - r_a(i,:);
    
    l_a(i) = norm(v_a(i,:));
    l_b(i) = norm(v_b(i,:));
    l_side(i) = norm(v_side(i,:));
    
    e_a(i,:) = v_a(i,:) / l_a(i);
    e_b(i,:) = v_b(i,:) / l_b(i);
    e_side(i,:) = v_side(i,:) / l_side(i);
end

x_init = [r_a(:); r_b(:)];
E_prev = compute_total_energy(x_init, n, params);

%% ========== 4. Visualization Setup ==========
figure('Position', [100, 100, 1500, 500]);

% Left: Geometry
subplot(1,3,1); hold on; grid on; axis equal;
h_inner = plot([r_a(:,1); r_a(1,1)], [r_a(:,2); r_a(1,2)], 'ro-', 'LineWidth', 2);
h_outer = plot([r_b(:,1); r_b(1,1)], [r_b(:,2); r_b(1,2)], 'bo-', 'LineWidth', 2);
h_side = zeros(n,1);
for i = 1:n
    h_side(i) = plot([r_a(i,1), r_b(i,1)], [r_a(i,2), r_b(i,2)], 'g--', 'LineWidth', 1.2);
end
g_force_a = quiver(r_a(:,1), r_a(:,2), zeros(n,1), zeros(n,1), 0, 'y', 'AutoScale', 'off');
g_force_b = quiver(r_b(:,1), r_b(:,2), zeros(n,1), zeros(n,1), 0, 'g', 'AutoScale', 'off');
title('Organoid Geometry'); 
xlabel('x (m)');
ylabel('y (m)');
xlim([-1.5e-4, 1.5e-4]); 
ylim([-1.5e-4, 1.5e-4]);

% Middle: Energy
subplot(1,3,2); hold on; grid on;
h_E = plot(NaN, NaN, 'b-', 'LineWidth', 2);
xlabel('Step'); 
ylabel('Energy (J)'); 
title('Energy Convergence');
legend('Total Energy'); 
xlim([0, num_steps]);

% Right: Normalized Area
subplot(1,3,3); hold on; grid on;
h_area = plot(NaN, NaN, 'm-', 'LineWidth', 2);
yline(1, 'k--', 'Initial Area');
xlabel('Step'); 
ylabel('Normalized Area'); 
title('Area Evolution');
xlim([0, num_steps]); 
ylim([0.5, 1.5]);

% Step Counter Annotation
time_text = annotation('textbox', [0.35, 0.98, 0.3, 0.03], ...
    'String', 'Step: 0', 'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'BackgroundColor', 'w', 'EdgeColor', 'none');

%% ========== 5. Video Initialization ==========
videoFileName = 'organoid_evolution.mp4';
if exist(videoFileName, 'file'), delete(videoFileName); end
video = VideoWriter(videoFileName, 'MPEG-4');
video.FrameRate = 60;
open(video);

energy_log = zeros(num_steps, 1);
area_log = zeros(num_steps, 1);

%% ========== 6. Main Simulation Loop ==========
for step = 1:num_steps
    
    % --- 6.1 Update Geometry & Pressures ---
    P_ecm = zeros(n,1); P_water = zeros(n,1);
    for i = 1:n
        next_idx = mod(i, n) + 1;
        v_a(i,:) = r_a(next_idx,:) - r_a(i,:);
        v_b(i,:) = r_b(next_idx,:) - r_b(i,:);
        v_side(i,:) = r_b(i,:) - r_a(i,:);
        
        % Pressure Models
        P_ecm(i) = 500/(77e-6)*norm(r_b(i,:)) - 400;
        P_water(i) = 800*sqrt(66e-6)/sqrt(norm(r_a(i,:))) - 700;
        
        % Transient Pressure Increase (First 10k steps logic from original)
        if step < 10000
             P_water(i) = P_water(i) + 100;% Hypotonic condition
             % P_water(i) = P_water(i) + 200;% Hypotonic condition (stronger)
             % P_water(i) = P_water(i) - 200;% Hyperosmotic condition
             % P_water(i) = P_water(i) - 500;% Hyperosmotic condition (stronger)
        end
        
        % Update Lengths & Unit Vectors
        l_a(i) = norm(v_a(i,:));
        l_b(i) = norm(v_b(i,:));
        l_side(i) = norm(v_side(i,:));
        e_a(i,:) = v_a(i,:) / l_a(i);
        e_b(i,:) = v_b(i,:) / l_b(i);
        e_side(i,:) = v_side(i,:) / l_side(i);
    end

    % --- 6.2 Bending Forces ---
    F_bend_a = zeros(n,2); 
    F_bend_b = zeros(n,2);
    target_cos = 0; 
    bend_coeff = 2 * kappa_bend;
    
    for i = 1:n
        next_idx = mod(i, n) + 1;
        r_ai=r_a(i,:); 
        r_an=r_a(next_idx,:); 
        r_bi=r_b(i,:); 
        r_bn=r_b(next_idx,:);
        
        % Vectors
        u_in = r_an-r_ai; 
        u_out = r_bn-r_bi;
        w_i = r_bi-r_ai; 
        w_next = r_bn-r_an;
        
        len_uin=norm(u_in); 
        len_uout=norm(u_out); 
        len_win=norm(w_i); 
        len_wnext=norm(w_next);
        if min([len_uin,len_uout,len_win,len_wnext]) < 1e-12, continue; end
        
        e_uin=u_in/len_uin; 
        e_uout=u_out/len_uout; 
        e_win=w_i/len_win; 
        e_wnext=w_next/len_wnext;
        
        % Corner 1 (Inner Start)
        c1 = max(-1,min(1,dot(e_win,e_uin)));
        if abs(c1) > 1e-8
            f = -bend_coeff*c1;
            gw = (e_uin-c1*e_win)/len_win; 
            gu = (e_win-c1*e_uin)/len_uin;
            F_bend_b(i,:) = F_bend_b(i,:) + f*gw; 
            F_bend_a(i,:) = F_bend_a(i,:) - f*gw;
            F_bend_a(next_idx,:) = F_bend_a(next_idx,:) + f*gu; 
            F_bend_a(i,:) = F_bend_a(i,:) - f*gu;
        end
        % Corner 2 (Outer Start)
        c2 = max(-1,min(1,dot(-e_win,e_uout)));
        if abs(c2) > 1e-8
            f = -bend_coeff*c2;
            gw = (e_uout-c2*(-e_win))/len_win; 
            gu = ((-e_win)-c2*e_uout)/len_uout;
            F_bend_a(i,:) = F_bend_a(i,:) + f*gw; 
            F_bend_b(i,:) = F_bend_b(i,:) - f*gw;
            F_bend_b(next_idx,:) = F_bend_b(next_idx,:) + f*gu; 
            F_bend_b(i,:) = F_bend_b(i,:) - f*gu;
        end
        % Corner 3 (Inner End)
        c3 = max(-1,min(1,dot(-e_uin,e_wnext)));
        if abs(c3) > 1e-8
            f = -bend_coeff*c3;
            gu = (e_wnext-c3*(-e_uin))/len_uin; 
            gn = ((-e_uin)-c3*e_wnext)/len_wnext;
            F_bend_a(i,:) = F_bend_a(i,:) + f*gu; 
            F_bend_a(next_idx,:) = F_bend_a(next_idx,:) - f*gu;
            F_bend_b(next_idx,:) = F_bend_b(next_idx,:) + f*gn;
            F_bend_a(next_idx,:) = F_bend_a(next_idx,:) - f*gn;
        end
        % Corner 4 (Outer End)
        c4 = max(-1,min(1,dot(-e_uout,e_wnext)));
        if abs(c4) > 1e-8
            f = -bend_coeff*c4;
            gu = (e_wnext-c4*(-e_uout))/len_uout; 
            gn = ((-e_uout)-c4*e_wnext)/len_wnext;
            F_bend_b(i,:) = F_bend_b(i,:) + f*gu; 
            F_bend_b(next_idx,:) = F_bend_b(next_idx,:) - f*gu;
            F_bend_b(next_idx,:) = F_bend_b(next_idx,:) + f*gn;
            F_bend_a(next_idx,:) = F_bend_a(next_idx,:) - f*gn;
        end
    end

    % --- 6.3 Surface Tension Gradients ---
    grad_U_a = zeros(n,2); 
    grad_U_b = zeros(n,2);
    grad_side_a = zeros(n,2); 
    grad_side_b = zeros(n,2);
    
    for i = 1:n
        next_idx = mod(i, n) + 1;
        % Apical
        grad_U_a(i,:) = grad_U_a(i,:) + sqrt(3)*gamma_a*l_a(i)*(-e_a(i,:));
        grad_U_a(next_idx,:) = grad_U_a(next_idx,:) + sqrt(3)*gamma_a*l_a(i)*e_a(i,:);
        % Basal
        grad_U_b(i,:) = grad_U_b(i,:) + sqrt(3)*gamma_b*l_b(i)*(-e_b(i,:));
        grad_U_b(next_idx,:) = grad_U_b(next_idx,:) + sqrt(3)*gamma_b*l_b(i)*e_b(i,:);
        % Lateral
        gs_i = gamma_l/sqrt(12);
        grad_side_a(i,:) = grad_side_a(i,:) - e_a(i,:)*gs_i*(l_side(i)+l_side(next_idx)) - e_side(i,:)*gs_i*(l_a(i)+l_b(i));
        grad_side_a(next_idx,:) = grad_side_a(next_idx,:) + e_a(i,:)*gs_i*(l_side(i)+l_side(next_idx)) - e_side(next_idx,:)*gs_i*(l_a(i)+l_b(i));
        grad_side_b(i,:) = grad_side_b(i,:) - e_b(i,:)*gs_i*(l_side(i)+l_side(next_idx)) + e_side(i,:)*gs_i*(l_a(i)+l_b(i));
        grad_side_b(next_idx,:) = grad_side_b(next_idx,:) + e_b(i,:)*gs_i*(l_side(i)+l_side(next_idx)) + e_side(next_idx,:)*gs_i*(l_a(i)+l_b(i));
    end
    F_surf_a = -grad_U_a; 
    F_surf_b = -grad_U_b;

    % --- 6.4 Pressure Forces ---
    F_water = zeros(n,2); 
    F_ecm = zeros(n,2);
    for i = 1:n
        prev_idx = mod(i-2, n) + 1;
        % Inner (Outward)
        t1 = v_a(i,:)/l_a(i); 
        t2 = v_a(prev_idx,:)/l_a(prev_idx);
        n1 = [t1(2),-t1(1)]; 
        e1 = n1/norm(n1); 
        n2 = [t2(2),-t2(1)]; 
        e2 = n2/norm(n2);
        F_water(i,:) = P_water(i)*sqrt(3)/2 * (l_a(prev_idx)^2*e2 + l_a(i)^2*e1);
        % Outer (Inward)
        t1 = v_b(i,:)/l_b(i); 
        t2 = v_b(prev_idx,:)/l_b(prev_idx);
        n1 = [-t1(2),t1(1)]; 
        e1 = n1/norm(n1); 
        n2 = [-t2(2),t2(1)]; 
        e2 = n2/norm(n2);
        F_ecm(i,:) = P_ecm(i)*sqrt(3)/2 * (l_b(prev_idx)^2*e2 + l_b(i)^2*e1);
    end

    % --- 6.5 Area Constraint Forces ---
    F_area_a = zeros(n,2); F_area_b = zeros(n,2);
    for i = 1:n
        next_idx = mod(i, n) + 1;
        P1=r_a(i,:); 
        P2=r_b(i,:); 
        P3=r_b(next_idx,:); 
        P4=r_a(next_idx,:);
        A_curr = quad_area(P1,P2,P3,P4);
        dA = A_curr - A0;
        dA_dP1 = -0.5*[-(P2(2)-P4(2)), P2(1)-P4(1)];
        dA_dP2 = -0.5*[-(P3(2)-P1(2)), P3(1)-P1(1)];
        dA_dP3 = -0.5*[-(P4(2)-P2(2)), P4(1)-P2(1)];
        dA_dP4 = -0.5*[-(P1(2)-P3(2)), P1(1)-P3(1)];
        
        F_area_a(i,:) = F_area_a(i,:) - K_area*dA*dA_dP1;
        F_area_b(i,:) = F_area_b(i,:) - K_area*dA*dA_dP2;
        F_area_b(next_idx,:) = F_area_b(next_idx,:) - K_area*dA*dA_dP3;
        F_area_a(next_idx,:) = F_area_a(next_idx,:) - K_area*dA*dA_dP4;
    end

    % --- 6.6 Gaussian Curvature Forces ---
    F_gauss_a = zeros(n,2); F_gauss_b = zeros(n,2);
    for i = 1:n
        next_idx = mod(i, n) + 1; prev_idx = mod(i-2, n) + 1;
        F_gauss_a(i,:) = -(A_gauss*2*(-e_a(i,:))/(-l_b(i)*l_a(i)^2) + A_gauss*(-e_side(i,:))/(-l_side(i)^2*l_side(next_idx))...
            +A_gauss*2*(e_a(prev_idx,:))/(-l_b(prev_idx)*l_a(prev_idx)^2) + A_gauss*(-e_side(i,:))/(-l_side(i)^2*l_side(prev_idx)));
        F_gauss_b(i,:) = -(A_gauss*2*(-e_b(i,:))/(-l_a(i)*l_b(i)^2) + A_gauss*(e_side(i,:))/(-l_side(i)^2*l_side(next_idx))...
            +A_gauss*2*(e_b(prev_idx,:))/(-l_a(prev_idx)*l_b(prev_idx)^2) + A_gauss*(e_side(i,:))/(-l_side(i)^2*l_side(prev_idx)));
    end

    % --- 6.7 Dynamics & Position Update ---
    vel_a = zeros(n,2); vel_b = zeros(n,2);
    for i = 1:n
        F_tot_a = F_area_a(i,:) + F_surf_a(i,:) - grad_side_a(i,:) + F_water(i,:) + F_bend_a(i,:) + F_gauss_a(i,:);
        F_tot_b = F_area_b(i,:) + F_surf_b(i,:) - grad_side_b(i,:) + F_ecm(i,:) + F_bend_b(i,:) + F_gauss_b(i,:);
        
        vel_a(i,:) = F_tot_a / mu;
        vel_b(i,:) = F_tot_b / mu;
        
        r_a(i,:) = r_a(i,:) + vel_a(i,:)*dt;
        r_b(i,:) = r_b(i,:) + vel_b(i,:)*dt;
    end

    % --- 6.8 Energy & Area Logging ---
    x_new = [r_a(:); r_b(:)];
    E_curr = compute_total_energy(x_new, n, params) - W_ecm - W_water; 
    dE = E_curr - E_prev;
    energy_log(step) = E_curr;
    E_prev = E_curr;
    
    % Calculate Normalized Area
    total_area = polyarea(r_b(:,1), r_b(:,2));
    area_log(step) = total_area / total_A0;

    % --- 6.9 Convergence Check ---
    max_vel = max([max(sqrt(sum(vel_a.^2,2))); max(sqrt(sum(vel_b.^2,2)))]);
    if max_vel < 1e-10 && step > 100
        fprintf('Converged at step %d\n', step);
        break;
    end

    % --- 6.10 Visualization ---
    if mod(step, vis_interval) == 0 || step == 1
        set(h_inner, 'XData', [r_a(:,1); r_a(1,1)], 'YData', [r_a(:,2); r_a(1,2)]);
        set(h_outer, 'XData', [r_b(:,1); r_b(1,1)], 'YData', [r_b(:,2); r_b(1,2)]);
        for i = 1:n
            set(h_side(i), 'XData', [r_a(i,1), r_b(i,1)], 'YData', [r_a(i,2), r_b(i,2)]);
        end
        
        set(h_E, 'XData', 1:step, 'YData', energy_log(1:step));
        set(h_area, 'XData', 1:step, 'YData', area_log(1:step));
        
        % Update Force Arrows (Scaled for visibility)
        set(g_force_a, 'UData', F_bend_a(:,1)*1000, 'VData', F_bend_a(:,2)*1000);
        set(g_force_b, 'UData', F_bend_b(:,1)*1000, 'VData', F_bend_b(:,2)*1000);
        
        set(time_text, 'String', sprintf('Step: %d', step));
        
        drawnow limitrate;
        writeVideo(video, getframe(gcf));
        
        % Console Output
        status = 'OK'; if dE > 0, status = 'WARN: E Incr'; end
        fprintf('%-6d | %.4e | %.4e | %-10s |\n', step, E_curr, dE, status);
    end
end

close(video);
disp(['Animation saved: ', videoFileName]);
disp('Simulation Complete!');