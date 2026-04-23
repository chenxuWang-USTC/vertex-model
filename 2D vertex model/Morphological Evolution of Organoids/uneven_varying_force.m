%% SIMULATE_ANNULAR_TRAPEZOID_LANGEVIN_MULTI_CASE
% Based on viscous dynamics equation: mu * dr/dt = -grad(U) + F_water + F_ecm + F_vol


clear;
clc;

    %% ========== 1. Base Parameters (Common to all Cases) ==========
    n = 40;                    % Number of trapezoidal units (cells)
    dt = 0.2;                  % Time step
    num_steps = 10000;          % Total simulation steps
    mu = 0.4;                  % Viscous damping coefficient
    
    % Base Intrinsic Tensions (Default uniform values)
    gamma_a_base = 2e-3;       % Apical tension 
    gamma_b_base = -1e-3;      % Basal tension 
    gamma_l_base = -3e-3;      % Lateral tension 
    
    A = 5e-23;                 % Gaussian curvature energy parameter
    K_area = 1e9;              % Area constraint stiffness
    kappa_bend = 1e-12;        % Bending stiffness
    
    vis_interval = 5;          % Visualization update interval
    W_water = 0;               % Work done by water pressure
    W_ecm = 0;                 % Work done by ECM pressure

    %% ========== CASE CONFIGURATION ==========
    
    % --- Initialize local tension vectors to default uniform values ---
    gamma_a_local = ones(n, 1) * gamma_a_base;
    gamma_b_local = ones(n, 1) * gamma_b_base;
    gamma_l_local = ones(n, 1) * gamma_l_base;
    
    % Define angular sector boundaries (for spatial variation)
    angle_edges = [0, pi/8, pi/4, 3*pi/8, 2*pi/4, 5*pi/8, 3*pi/4, 7*pi/8, pi, 9*pi/8, 5*pi/4, 11*pi/8, 6*pi/4, 13/pi/8, 7*pi/4, 15*pi/8, 2*pi];
    num_regions = length(angle_edges) - 1;
    
    % ======================================================================
    % [Case 1: Vary Basal Tension Gamma_b (Spatially Non-uniform)]
    % ======================================================================
    case_id = 1; 
    videoFileName = 'case1_vary_gamma_b.mp4';
    
    % Define spatial pattern for Gamma_b (16 sectors)
    gamma_b_pattern = [20e-3, 2e-3, 20e-3, 2e-3, 10e-3, 2e-3, 20e-3, 2e-3, 20e-3, 2e-3, 2e-3, 10e-3, 2e-3, 20e-3, 2e-3, 10e-3];
    
    % Apply pattern to gamma_b_local
    % (Note: gamma_a_local and gamma_l_local remain uniform defaults)


    % ======================================================================
    % [Case 2: Vary Apical Tension Gamma_a (Spatially Non-uniform)]
    % ======================================================================
    % case_id = 2; 
    % videoFileName = 'case2_vary_gamma_a.mp4';
    % 
    % % Define spatial pattern for Gamma_a 
    % angle_edges_6 = linspace(0, 2*pi, 7);
    % gamma_a_pattern_6 = [10e-3, 2e-3, 10e-3, 2e-3, 10e-3, 2e-3];
    % 
    % % Placeholder: Actual assignment happens after geometry initialization below
    % is_case_2 = true; 


    % ======================================================================
    % [Case 3: Vary Lateral Tension Gamma_l (Spatially Non-uniform)]
    % ======================================================================
    % case_id = 3; 
    % videoFileName = 'case3_vary_gamma_l.mp4';
    % 
    % % Define spatial pattern for Gamma_l 
    % % Top/Bottom quadrants: High Tension; Left/Right: Low Tension
    % is_case_3 = true;


    %% ========== 2. Initial Geometry ==========
    L = intital(); % Ensure 'intital.m' is in path
    r_a = L.r_a;   % Inner ring vertices
    r_b = L.r_b;   % Outer ring vertices
    A0 = L.A0;     % Reference area per cell
    
    % === Key: Symmetrize Initial Positions ===
    fprintf('=== Symmetrizing Initial Geometry ===\n');
    R_a_mean = mean(sqrt(sum(r_a.^2, 2)));
    R_b_mean = mean(sqrt(sum(r_b.^2, 2)));

    n = length(r_a);
    theta = linspace(0, 2*pi, n+1)'; 
    theta(end) = []; % Remove duplicate endpoint

    % Reset positions to perfect circles with mean radii
    r_a = R_a_mean * [cos(theta), sin(theta)];
    r_b = R_b_mean * [cos(theta), sin(theta)];

    % === Apply Spatial Tension Distribution based on Selected Case ===
    theta_b = atan2(r_b(:,2), r_b(:,1));
    theta_b = mod(theta_b, 2*pi);  % Normalize to [0, 2π)

    if case_id == 1
        % Apply Case 1: Vary Gamma_b
        for i = 1:n
            for reg = 1:num_regions
                if theta_b(i) >= angle_edges(reg) && theta_b(i) < angle_edges(reg+1)
                    gamma_b_local(i) = gamma_b_pattern(reg);
                    break;
                end
            end
        end
        fprintf('Running Case 1: Spatially Varying Gamma_b\n');
        
    elseif case_id == 2
        % Apply Case 2: Vary Gamma_a (6 sectors)
        angle_edges_6 = linspace(0, 2*pi, 7);
        gamma_a_pattern_6 = [10e-3, 2e-3, 10e-3, 2e-3, 10e-3, 2e-3];
        for i = 1:n
            for reg = 1:6
                if theta_b(i) >= angle_edges_6(reg) && theta_b(i) < angle_edges_6(reg+1)
                    gamma_a_local(i) = gamma_a_pattern_6(reg);
                    break;
                end
            end
        end
        fprintf('Running Case 2: Spatially Varying Gamma_a\n');
        
    elseif case_id == 3
        % Apply Case 3: Vary Gamma_l (Quadrants)
        for i = 1:n
            % Top (pi/4 to 3pi/4) and Bottom (5pi/4 to 7pi/4) are High Tension
            if (theta_b(i) >= pi/4 && theta_b(i) < 3*pi/4) || ...
               (theta_b(i) >= 5*pi/4 && theta_b(i) < 7*pi/4)
                gamma_l_local(i) = 20e-3; % High
            else
                gamma_l_local(i) = 2e-3;  % Low
            end
        end
        fprintf('Running Case 3: Spatially Varying Gamma_l\n');
    end

    % Pack parameters for energy calculation (if needed)
    params.gamma_a = gamma_a_base; 
    params.gamma_b = gamma_b_base;
    params.gamma_l = gamma_l_base;
    params.K_area = K_area;
    params.A = A;
    params.kappa_bend = kappa_bend;
    params.A0 = A0;
    % Note: Local variations are used directly in force loops below

    x = [r_a(:), r_b(:)];
    E_prev = compute_total_energy(x, n, params); % Optional initial energy check

    %% ========== 3. Initialize Edge Vectors ==========
    v_a = zeros(n, 2); 
    v_b = zeros(n, 2); 
    v_side = zeros(n, 2);
    for i = 1:n
        next_idx = mod(i, n) + 1;
        v_a(i,:) = r_a(next_idx,:) - r_a(i,:);
        v_b(i,:) = r_b(next_idx,:) - r_b(i,:);
        v_side(i,:) = r_b(i,:) - r_a(i,:); % Side vector: Inner -> Outer
    end

    %% ========== 4. Pre-compute Initial Geometric Quantities ==========
    l_a = zeros(n,1); 
    l_b = zeros(n,1); 
    l_side = zeros(n,1);
    e_a = zeros(n,2); 
    e_b = zeros(n,2); 
    e_side = zeros(n,2);
    for i = 1:n
        l_a(i) = norm(v_a(i,:));
        l_b(i) = norm(v_b(i,:));
        l_side(i) = norm(v_side(i,:));
    end
    e_a = v_a ./ l_a;
    e_b = v_b ./ l_b;
    e_side = v_side ./ l_side;

    %% ========== 5. Visualization Setup ==========
    figure('Position', [100, 100, 1400, 700]); 
    subplot(1,2,1); hold on; grid on;
    h_inner = plot([r_a(:,1); r_a(1,1)], [r_a(:,2); r_a(1,2)], 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
    h_outer = plot([r_b(:,1); r_b(1,1)], [r_b(:,2); r_b(1,2)], 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    
    h_side = zeros(n, 1);
    for i = 1:n
        h_side(i) = plot([r_a(i,1), r_b(i,1)], [r_a(i,2), r_b(i,2)], 'g--', 'LineWidth', 1.2);
    end

    g_force_a = quiver(r_a(:,1), r_a(:,2), zeros(n,1), zeros(n,1), 0, 'y', 'AutoScale', 'off');
    g_force_b = quiver(r_b(:,1), r_b(:,2), zeros(n,1), zeros(n,1), 0, 'g', 'AutoScale', 'off');

    % force vector
    % h_force_a = quiver(r_a(:,1), r_a(:,2), zeros(n,1), zeros(n,1), 0, 'r', 'MaxHeadSize', 0.5, 'LineWidth', 1.5, 'AutoScale', 'off');
    % h_force_b = quiver(r_b(:,1), r_b(:,2), zeros(n,1), zeros(n,1), 0, 'b', 'MaxHeadSize', 0.5, 'LineWidth', 1.5, 'AutoScale', 'off');


    title(['Case ', num2str(case_id), ': Organoid Simulation'], 'FontSize', 16, 'FontWeight', 'bold');
    xlabel('x', 'FontSize', 14); ylabel('y', 'FontSize', 14);
    axis equal; grid on; 
    xlim([-1.5e-4, 1.5e-4]); 
    ylim([-1.5e-4, 1.5e-4]);

    % Energy analysis
    subplot(1,2,2); hold on; grid on;
    h_E = plot(NaN, NaN, 'b-', 'LineWidth', 2);
    h_dE = plot(NaN, NaN, 'r--', 'LineWidth', 1.5);
    xlabel('Step'); ylabel('Energy (J)');
    title('Energy analysis');
    legend('Total Energy');
    ylim auto;

    energy_log = [];
    dE_log = [];

    %% ========== 6. Video Initialization ==========
    if exist(videoFileName, 'file'), delete(videoFileName); end
    video = VideoWriter(videoFileName, 'MPEG-4');
    video.FrameRate = 60;
    open(video);

    %% ========== 7. Main Simulation Loop ==========
    for step = 1:num_steps
        
        % --- 7.1 Update Edge Vectors & Pressures ---
        P_ecm = zeros(n,1);
        P_water = zeros(n,1);
        for i = 1:n
            next_idx = mod(i, n) + 1;
            v_a(i,:) = r_a(next_idx,:) - r_a(i,:);
            v_b(i,:) = r_b(next_idx,:) - r_b(i,:);
            v_side(i,:) = r_b(i,:) - r_a(i,:); 
            
            % Pressure Models
            P_ecm(i) = 500/(77e-6)*norm(r_b(i,:)) - 400;
            P_water(i) = 800*sqrt(66e-6)/sqrt(norm(r_a(i,:))) - 700;
        end
            
        % ====== 7.2 Update Geometric Quantities ======
        for i = 1:n
            l_a(i) = norm(v_a(i,:));
            l_b(i) = norm(v_b(i,:));
            l_side(i) = norm(v_side(i,:));
            e_side(i,:) = v_side(i,:) ./ l_side(i);
            e_a(i,:) = v_a(i,:) ./ l_a(i);
            e_b(i,:) = v_b(i,:) ./ l_b(i);
        end

        %% ====== 7.3 Bending Stiffness Forces ======
        F_bend_a = zeros(n, 2);
        F_bend_b = zeros(n, 2);
        E_bend = 0;
        target_cos = 0; % Target: 90 degrees (cos=0)
        bend_coeff = 2 * kappa_bend;

        for i = 1:n
            next_idx = mod(i, n) + 1;
            r_ai = r_a(i,:); 
            r_ain = r_a(next_idx,:);
            r_bi = r_b(i,:); 
            r_bin = r_b(next_idx,:);
            
            u_inner = r_ain - r_ai; 
            u_outer = r_bin - r_bi;
            w_side_i = r_bi - r_ai; 
            w_side_next = r_bin - r_ain;
            
            len_uin = norm(u_inner); 
            len_uout = norm(u_outer);
            len_win = norm(w_side_i); 
            len_wnext = norm(w_side_next);
            
            if min([len_uin, len_uout, len_win, len_wnext]) < 1e-12
                continue; 
            end
            
            e_uin = u_inner/len_uin; 
            e_uout = u_outer/len_uout;
            e_win = w_side_i/len_win; 
            e_wnext = w_side_next/len_wnext;
            
            % Corner 1: Inner Start
            cos1 = max(-1, min(1, dot(e_win, e_uin)));
            E_bend = E_bend + 0.5*bend_coeff*(cos1)^2;
            if abs(cos1) > 1e-8
                f = -bend_coeff*cos1;
                gw = (e_uin - cos1*e_win)/len_win; 
                gu = (e_win - cos1*e_uin)/len_uin;
                F_bend_b(i,:) = F_bend_b(i,:) + f*gw; 
                F_bend_a(i,:) = F_bend_a(i,:) - f*gw;
                F_bend_a(next_idx,:) = F_bend_a(next_idx,:) + f*gu; 
                F_bend_a(i,:) = F_bend_a(i,:) - f*gu;
            end
            % Corner 2: Outer Start
            cos2 = max(-1, min(1, dot(-e_win, e_uout)));
            E_bend = E_bend + 0.5*bend_coeff*(cos2)^2;
            if abs(cos2) > 1e-8
                f = -bend_coeff*cos2;
                gw = (e_uout - cos2*(-e_win))/len_win; 
                gu = ((-e_win) - cos2*e_uout)/len_uout;
                F_bend_a(i,:) = F_bend_a(i,:) + f*gw; 
                F_bend_b(i,:) = F_bend_b(i,:) - f*gw;
                F_bend_b(next_idx,:) = F_bend_b(next_idx,:) + f*gu; 
                F_bend_b(i,:) = F_bend_b(i,:) - f*gu;
            end
            % Corner 3: Inner End
            cos3 = max(-1, min(1, dot(-e_uin, e_wnext)));
            E_bend = E_bend + 0.5*bend_coeff*(cos3)^2;
            if abs(cos3) > 1e-8
                f = -bend_coeff*cos3;
                gu = (e_wnext - cos3*(-e_uin))/len_uin; 
                gn = ((-e_uin) - cos3*e_wnext)/len_wnext;
                F_bend_a(i,:) = F_bend_a(i,:) + f*gu; 
                F_bend_a(next_idx,:) = F_bend_a(next_idx,:) - f*gu;
                F_bend_b(next_idx,:) = F_bend_b(next_idx,:) + f*gn; 
                F_bend_a(next_idx,:) = F_bend_a(next_idx,:) - f*gn;
            end
            % Corner 4: Outer End
            cos4 = max(-1, min(1, dot(-e_uout, e_wnext)));
            E_bend = E_bend + 0.5*bend_coeff*(cos4)^2;
            if abs(cos4) > 1e-8
                f = -bend_coeff*cos4;
                gu = (e_wnext - cos4*(-e_uout))/len_uout; 
                gn = ((-e_uout) - cos4*e_wnext)/len_wnext;
                F_bend_b(i,:) = F_bend_b(i,:) + f*gu; 
                F_bend_b(next_idx,:) = F_bend_b(next_idx,:) - f*gu;
                F_bend_b(next_idx,:) = F_bend_b(next_idx,:) + f*gn; 
                F_bend_a(next_idx,:) = F_bend_a(next_idx,:) - f*gn;
            end
        end

        %% ====== 7.4 Surface Tension Gradients (Using LOCAL Variables) =====
        grad_U_a = zeros(n, 2); 
        grad_U_b = zeros(n, 2);
        grad_side_a = zeros(n, 2); 
        grad_side_b = zeros(n, 2);
        
        for i = 1:n
            next_idx = mod(i, n) + 1;

            % Apical: Use gamma_a_local(i)
            grad_U_a(i,:) = grad_U_a(i,:) + sqrt(3) * gamma_a_local(i) * l_a(i) * (-e_a(i,:));
            grad_U_a(next_idx,:) = grad_U_a(next_idx,:) + sqrt(3) * gamma_a_local(i) * l_a(i) * e_a(i,:);
            
            % Basal: Use gamma_b_local(i)
            grad_U_b(i,:) = grad_U_b(i,:) + sqrt(3) * gamma_b_local(i) * l_b(i) * (-e_b(i,:));
            grad_U_b(next_idx,:) = grad_U_b(next_idx,:) + sqrt(3) * gamma_b_local(i) * l_b(i) * e_b(i,:);

            % Lateral: Use gamma_l_local(i)
            grad_side_a(i,:) = grad_side_a(i,:) ...
                - e_a(i,:) * gamma_l_local(i)/sqrt(12)*(l_side(i)+l_side(next_idx)) ...
                - e_side(i,:) * gamma_l_local(i)/sqrt(12)*(l_a(i)+l_b(i));
        
            grad_side_a(next_idx,:) = grad_side_a(next_idx,:) ...
                + e_a(i,:) * gamma_l_local(i)/sqrt(12)*(l_side(i)+l_side(next_idx)) ...
                - e_side(next_idx,:) * gamma_l_local(i)/sqrt(12)*(l_a(i)+l_b(i));
    
            grad_side_b(i,:) = grad_side_b(i,:) ...
                - e_b(i,:) * gamma_l_local(i)/sqrt(12)*(l_side(i)+l_side(next_idx)) ...
                + e_side(i,:) * gamma_l_local(i)/sqrt(12)*(l_a(i)+l_b(i));
    
            grad_side_b(next_idx,:) = grad_side_b(next_idx,:) ...
                + e_b(i,:) * gamma_l_local(i)/sqrt(12)*(l_side(i)+l_side(next_idx)) ...
                + e_side(next_idx,:) * gamma_l_local(i)/sqrt(12)*(l_a(i)+l_b(i));
        end
        F_surface_a = -grad_U_a;
        F_surface_b = -grad_U_b;
        
        %% ===== 7.5 Pressure Forces (Unchanged) =====
        F_water = zeros(n, 2); 
        F_ecm = zeros(n, 2);
        for i = 1:n
            prev_idx = mod(i-2, n) + 1;
            % Internal Pressure (Outward)
            v_prev = v_a(prev_idx,:); v_i = v_a(i,:);
            t1 = v_i/norm(v_i); 
            t2 = v_prev/norm(v_prev);
            n1 = [t1(2), -t1(1)]; 
            e1 = n1/norm(n1);
            n2 = [t2(2), -t2(1)]; 
            e2 = n2/norm(n2);
            F_water(i,:) = P_water(i)*sqrt(3)/2*(l_a(prev_idx)^2*e2 + l_a(i)^2*e1);
           
            % External Pressure (Inward)
            v_prev = v_b(prev_idx,:); v_i = v_b(i,:);
            t1 = v_i/norm(v_i); 
            t2 = v_prev/norm(v_prev);
            n1 = [-t1(2), t1(1)]; 
            e1 = n1/norm(n1);
            n2 = [-t2(2), t2(1)]; 
            e2 = n2/norm(n2);
            F_ecm(i,:) = P_ecm(i)*sqrt(3)/2*(l_b(prev_idx)^2*e2 + l_b(i)^2*e1);
        end

        %% ===== 7.6 Area Constraint Forces (Unchanged) =====
        F_area_a = zeros(n, 2); 
        F_area_b = zeros(n, 2);
        for i = 1:n
            next_idx = mod(i, n) + 1;
            P1=r_a(i,:); 
            P2=r_b(i,:); 
            P3=r_b(next_idx,:); 
            P4=r_a(next_idx,:);
            A_curr = quad_area(P1, P2, P3, P4);
            dA = A_curr - A0;
            
            % Analytical gradients of area w.r.t vertices
            dA_dP1 = -0.5*[-(P2(2)-P4(2)), P2(1)-P4(1)];
            dA_dP2 = -0.5*[-(P3(2)-P1(2)), P3(1)-P1(1)];
            dA_dP3 = -0.5*[-(P4(2)-P2(2)), P4(1)-P2(1)];
            dA_dP4 = -0.5*[-(P1(2)-P3(2)), P1(1)-P3(1)];
            
            F_area_a(i,:) = F_area_a(i,:) - K_area*dA*dA_dP1;
            F_area_b(i,:) = F_area_b(i,:) - K_area*dA*dA_dP2;
            F_area_b(next_idx,:) = F_area_b(next_idx,:) - K_area*dA*dA_dP3;
            F_area_a(next_idx,:) = F_area_a(next_idx,:) - K_area*dA*dA_dP4;
        end

        %% --- 7.7 Gaussian Curvature Constraint Forces (Unchanged) ---
        F_gauss_a = zeros(n,2); 
        F_gauss_b = zeros(n,2);
        for i=1:n
            next_idx = mod(i, n) + 1;
            prev_idx = mod(i-2, n) + 1;
            F_gauss_a(i,:) = -(A*2*(-e_a(i,:))/(-l_b(i)*l_a(i)^2) + A*(-e_side(i,:))/(-l_side(i)^2*l_side(next_idx))...
                +A*2*(e_a(prev_idx,:))/(-l_b(prev_idx)*l_a(prev_idx)^2) + A*(-e_side(i,:))/(-l_side(i)^2*l_side(prev_idx)));
            F_gauss_b(i,:) = -(A*2*(-e_b(i,:))/(-l_a(i)*l_b(i)^2) + A*(e_side(i,:))/(-l_side(i)^2*l_side(next_idx))...
                +A*2*(e_b(prev_idx,:))/(-l_a(prev_idx)*l_b(prev_idx)^2) + A*(e_side(i,:))/(-l_side(i)^2*l_side(prev_idx)));
        end

        %% --- 7.8 Total Force & Dynamics ---
        vel_a = zeros(n,2); vel_b = zeros(n,2);
        for i=1:n
            F_tot_a = F_area_a(i,:) + F_surface_a(i,:) - grad_side_a(i,:) + F_water(i,:) + F_bend_a(i,:) + F_gauss_a(i,:);
            F_tot_b = F_area_b(i,:) + F_surface_b(i,:) - grad_side_b(i,:) + F_ecm(i,:) + F_bend_b(i,:) + F_gauss_b(i,:);
            
            vel_a(i,:) = F_tot_a / mu;
            vel_b(i,:) = F_tot_b / mu;
        end
  
        % --- 7.9 Update Positions ---
        for i=1:n
            r_a(i,:) = r_a(i,:) + vel_a(i,:) * dt;
            W_water = W_water + dot(vel_a(i,:),F_water(i,:)) * dt;
            r_b(i,:) = r_b(i,:) + vel_b(i,:) * dt;
            W_ecm = W_ecm +  dot(vel_b(i,:),F_ecm(i,:)) * dt  ;
        end
        x_new= [r_a ,r_b];
        for i=1:1:n
            params.gamma_b_local(i)=gamma_b_local(i);
        end
        E_curr = compute_total_energy(x_new,n,params) + E_bend - W_ecm - W_water;
        dE = E_curr - E_prev ;%
        if step >= 1
            energy_log(step) = E_curr;
            dE_log(step) = dE;
        end
        E_prev = E_curr;
        % --- 7.10 Convergence Check ---
        max_vel_a = max(sqrt(sum(vel_a.^2, 2)));  % Maximum speed of inner ring (m/s)
        max_vel_b = max(sqrt(sum(vel_b.^2, 2)));  % Maximum speed of outer ring (m/s)
    
        vel_tol = 3e-10;  % Speed convergence threshold
       if (max_vel_a < vel_tol) && (max_vel_b < vel_tol)
         % Update the image and save the last frame
        set(h_inner, 'XData', [r_a(:,1); r_a(1,1)], 'YData', [r_a(:,2); r_a(1,2)]);
        set(h_outer, 'XData', [r_b(:,1); r_b(1,1)], 'YData', [r_b(:,2); r_b(1,2)]);
        for i = 1:n
            set(h_side(i), 'XData', [r_a(i,1), r_b(i,1)], 'YData', [r_a(i,2), r_b(i,2)]);
        end
        writeVideo(video, getframe(gcf));   
        break;  
       end  

        % --- 7.11 Visualization Update ---
        if mod(step, vis_interval) == 0 || step == 1
            set(h_inner, 'XData', [r_a(:,1); r_a(1,1)], 'YData', [r_a(:,2); r_a(1,2)]);
            set(h_outer, 'XData', [r_b(:,1); r_b(1,1)], 'YData', [r_b(:,2); r_b(1,2)]);
            for i = 1:n
                set(h_side(i), 'XData', [r_a(i,1), r_b(i,1)], 'YData', [r_a(i,2), r_b(i,2)]);
            end
            
            % Update energy chart
            set(h_E, 'XData', 1:step, 'YData', energy_log(1:step)); 

            % Update force vector
            % g_norm_a = F_bend_a * 0; 
            % g_norm_b = F_bend_b * 0;
            % set(g_force_a, 'XData', r_a(:,1), 'YData', r_a(:,2), 'UData', g_norm_a(:,1), 'VData', g_norm_a(:,2));
            % set(g_force_b, 'XData', r_b(:,1), 'YData', r_b(:,2), 'UData', g_norm_b(:,1), 'VData', g_norm_b(:,2));

            time_str = sprintf('step：%d ', step);
            % 更新文本
            if ~exist('time_text_handle', 'var') || ~isvalid(time_text_handle)
                time_text_handle = text(0.05, 0.95, time_str, ...
                'Units', 'normalized', ...
                'FontSize', 20, ...
                'FontWeight', 'bold', ...
                'Color', 'k', ...
                'VerticalAlignment', 'top');
            else
                set(time_text_handle, 'String', time_str);
            end
            drawnow limitrate;
            writeVideo(video, getframe(gcf));
       end
    end

    close(video);
    disp(['Animation saved: ', videoFileName]);
    disp('Simulation Complete!');