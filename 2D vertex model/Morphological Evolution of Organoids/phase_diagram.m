%% PHASE_DIAGRAM_GAMMA_B_VS_PWATER
% Scan gamma_b (uniform) and P_water, compute steady-state total area
clear; clc;

%% ========== 1. Fixed Parameters ==========
n = 40;                     % number of cells
dt = 0.05;
num_steps = 10000;           % fixed steps for each simulation (ensure steady state)
mu = 0.4;

% Fixed tensions (except gamma_b which will be scanned)
gamma_a_base = 2e-3;        % apical tension
gamma_l_base = -3e-3;       % lateral tension

% Energy parameters
A = 1e-22;
K_area = 3e9;
kappa_bend = 1e-11;
tau_relax = 10000;
E_substrate=1300;

% Reference area from initial geometry
L = intital();
A0 =2.5e-10;                  % per cell

%% ========== 2. Scan Ranges ==========
gamma_b_vals = linspace(-0.001, 0.03, 20);   % basal tension scan (adjust range as needed)
P_water_vals = linspace(100, 800, 20);        % internal pressure scan

n_gamma = length(gamma_b_vals);
n_P = length(P_water_vals);
Area_mat = zeros(n_gamma, n_P);             % store total area

%% ========== 3. Main Scan Loop ==========
for i_g = 1:n_gamma
    gamma_b_val = gamma_b_vals(i_g);
    for i_p = 1:n_P
        P_water_val = P_water_vals(i_p);

        % --- Initialize geometry symmetrically ---
        r_a = L.r_a; r_b = L.r_b;
        R_a_mean = mean(sqrt(sum(r_a.^2, 2)));
        R_b_mean = mean(sqrt(sum(r_b.^2, 2)));
        theta = linspace(0, 2*pi, n+1)'; theta(end) = [];
        r_a = R_a_mean * [cos(theta), sin(theta)];
        r_b = R_b_mean * [cos(theta), sin(theta)];

        % --- Set uniform tensions ---
        gamma_a_local = ones(n,1) * gamma_a_base;
        gamma_b_local = ones(n,1) * gamma_b_val;   % uniform basal tension
        gamma_l_local = ones(n,1) * gamma_l_base;

        % --- Precompute edge vectors ---
        v_a = zeros(n,2); v_b = zeros(n,2); v_side = zeros(n,2);
        for k = 1:n
            next = mod(k,n)+1;
            v_a(k,:) = r_a(next,:) - r_a(k,:);
            v_b(k,:) = r_b(next,:) - r_b(k,:);
            v_side(k,:) = r_b(k,:) - r_a(k,:);
        end

        % --- Time stepping ---
        for step = 1:num_steps
            % Update edge lengths and unit vectors
            l_a = zeros(n,1); l_b = zeros(n,1); l_side = zeros(n,1);
            e_a = zeros(n,2); e_b = zeros(n,2); e_side = zeros(n,2);
            for k = 1:n
                next = mod(k,n)+1;
                v_a(k,:) = r_a(next,:) - r_a(k,:);
                v_b(k,:) = r_b(next,:) - r_b(k,:);
                v_side(k,:) = r_b(k,:) - r_a(k,:);
                l_a(k) = norm(v_a(k,:));
                l_b(k) = norm(v_b(k,:));
                l_side(k) = norm(v_side(k,:));
                e_a(k,:) = v_a(k,:) / l_a(k);
                e_b(k,:) = v_b(k,:) / l_b(k);
                e_side(k,:) = v_side(k,:) / l_side(k);
            end

            % Bending forces (same as original)
            F_bend_a = zeros(n,2); F_bend_b = zeros(n,2);
            bend_coeff = 2 * kappa_bend;
            for k = 1:n
                next = mod(k,n)+1;
                r_ai = r_a(k,:); r_ain = r_a(next,:);
                r_bi = r_b(k,:); r_bin = r_b(next,:);
                u_inner = r_ain - r_ai; u_outer = r_bin - r_bi;
                w_side_i = r_bi - r_ai; w_side_next = r_bin - r_ain;
                len_uin = norm(u_inner); len_uout = norm(u_outer);
                len_win = norm(w_side_i); len_wnext = norm(w_side_next);
                if min([len_uin, len_uout, len_win, len_wnext]) < 1e-12, continue; end
                e_uin = u_inner/len_uin; e_uout = u_outer/len_uout;
                e_win = w_side_i/len_win; e_wnext = w_side_next/len_wnext;

                % Corner 1
                cos1 = max(-1, min(1, dot(e_win, e_uin)));
                if abs(cos1) > 1e-8
                    f = -bend_coeff*cos1;
                    gw = (e_uin - cos1*e_win)/len_win;
                    gu = (e_win - cos1*e_uin)/len_uin;
                    F_bend_b(k,:) = F_bend_b(k,:) + f*gw;
                    F_bend_a(k,:) = F_bend_a(k,:) - f*gw;
                    F_bend_a(next,:) = F_bend_a(next,:) + f*gu;
                    F_bend_a(k,:) = F_bend_a(k,:) - f*gu;
                end
                % Corner 2
                cos2 = max(-1, min(1, dot(-e_win, e_uout)));
                if abs(cos2) > 1e-8
                    f = -bend_coeff*cos2;
                    gw = (e_uout - cos2*(-e_win))/len_win;
                    gu = ((-e_win) - cos2*e_uout)/len_uout;
                    F_bend_a(k,:) = F_bend_a(k,:) + f*gw;
                    F_bend_b(k,:) = F_bend_b(k,:) - f*gw;
                    F_bend_b(next,:) = F_bend_b(next,:) + f*gu;
                    F_bend_b(k,:) = F_bend_b(k,:) - f*gu;
                end
                % Corner 3
                cos3 = max(-1, min(1, dot(-e_uin, e_wnext)));
                if abs(cos3) > 1e-8
                    f = -bend_coeff*cos3;
                    gu = (e_wnext - cos3*(-e_uin))/len_uin;
                    gn = ((-e_uin) - cos3*e_wnext)/len_wnext;
                    F_bend_a(k,:) = F_bend_a(k,:) + f*gu;
                    F_bend_a(next,:) = F_bend_a(next,:) - f*gu;
                    F_bend_b(next,:) = F_bend_b(next,:) + f*gn;
                    F_bend_a(next,:) = F_bend_a(next,:) - f*gn;
                end
                % Corner 4
                cos4 = max(-1, min(1, dot(-e_uout, e_wnext)));
                if abs(cos4) > 1e-8
                    f = -bend_coeff*cos4;
                    gu = (e_wnext - cos4*(-e_uout))/len_uout;
                    gn = ((-e_uout) - cos4*e_wnext)/len_wnext;
                    F_bend_b(k,:) = F_bend_b(k,:) + f*gu;
                    F_bend_b(next,:) = F_bend_b(next,:) - f*gu;
                    F_bend_b(next,:) = F_bend_b(next,:) + f*gn;
                    F_bend_a(next,:) = F_bend_a(next,:) - f*gn;
                end
            end

            % Surface tension gradients
            grad_U_a = zeros(n,2); grad_U_b = zeros(n,2);
            grad_side_a = zeros(n,2); grad_side_b = zeros(n,2);
            for k = 1:n
                next = mod(k,n)+1;
                grad_U_a(k,:) = grad_U_a(k,:) + sqrt(3)*gamma_a_local(k)*l_a(k)*(-e_a(k,:));
                grad_U_a(next,:) = grad_U_a(next,:) + sqrt(3)*gamma_a_local(k)*l_a(k)*e_a(k,:);
                grad_U_b(k,:) = grad_U_b(k,:) + sqrt(3)*gamma_b_local(k)*l_b(k)*(-e_b(k,:));
                grad_U_b(next,:) = grad_U_b(next,:) + sqrt(3)*gamma_b_local(k)*l_b(k)*e_b(k,:);

                grad_side_a(k,:) = grad_side_a(k,:) ...
                    - e_a(k,:)*gamma_l_local(k)/sqrt(12)*(l_side(k)+l_side(next)) ...
                    - e_side(k,:)*gamma_l_local(k)/sqrt(12)*(l_a(k)+l_b(k));
                grad_side_a(next,:) = grad_side_a(next,:) ...
                    + e_a(k,:)*gamma_l_local(k)/sqrt(12)*(l_side(k)+l_side(next)) ...
                    - e_side(next,:)*gamma_l_local(k)/sqrt(12)*(l_a(k)+l_b(k));
                grad_side_b(k,:) = grad_side_b(k,:) ...
                    - e_b(k,:)*gamma_l_local(k)/sqrt(12)*(l_side(k)+l_side(next)) ...
                    + e_side(k,:)*gamma_l_local(k)/sqrt(12)*(l_a(k)+l_b(k));
                grad_side_b(next,:) = grad_side_b(next,:) ...
                    + e_b(k,:)*gamma_l_local(k)/sqrt(12)*(l_side(k)+l_side(next)) ...
                    + e_side(next,:)*gamma_l_local(k)/sqrt(12)*(l_a(k)+l_b(k));
            end
            F_surface_a = -grad_U_a;
            F_surface_b = -grad_U_b;

            % Pressure forces (uniform P_water)
            P_ecm = zeros(n,1);
            P_water = P_water_val * ones(n,1);  % constant internal pressure
            F_water = zeros(n,2); F_ecm = zeros(n,2);
            for k = 1:n
                prev = mod(k-2,n)+1;
                % internal pressure (outward on inner ring)
                v_prev = v_a(prev,:); v_i = v_a(k,:);
                t1 = v_i/norm(v_i); t2 = v_prev/norm(v_prev);
                n1 = [t1(2), -t1(1)]; e1 = n1/norm(n1);
                n2 = [t2(2), -t2(1)]; e2 = n2/norm(n2);
                F_water(k,:) = P_water(k)*sqrt(3)/2*(l_a(prev)^2*e2 + l_a(k)^2*e1);
                % external pressure (inward on outer ring) – using original P_ecm model
                P_ecm(k) = 500/(77e-6)*norm(r_b(k,:)) - 400;
                v_prev = v_b(prev,:); v_i = v_b(k,:);
                t1 = v_i/norm(v_i); t2 = v_prev/norm(v_prev);
                n1 = [-t1(2), t1(1)]; e1 = n1/norm(n1);
                n2 = [-t2(2), t2(1)]; e2 = n2/norm(n2);
                F_ecm(k,:) = P_ecm(k)*sqrt(3)/2*(l_b(prev)^2*e2 + l_b(k)^2*e1);
            end

            % Area constraint forces
            F_area_a = zeros(n,2); F_area_b = zeros(n,2);
            for k = 1:n
                next = mod(k,n)+1;
                P1=r_a(k,:); P2=r_b(k,:); P3=r_b(next,:); P4=r_a(next,:);
                A_curr = quad_area(P1,P2,P3,P4);
                dA = A_curr - A0;
                dA_dP1 = -0.5*[-(P2(2)-P4(2)), P2(1)-P4(1)];
                dA_dP2 = -0.5*[-(P3(2)-P1(2)), P3(1)-P1(1)];
                dA_dP3 = -0.5*[-(P4(2)-P2(2)), P4(1)-P2(1)];
                dA_dP4 = -0.5*[-(P1(2)-P3(2)), P1(1)-P3(1)];
                F_area_a(k,:) = F_area_a(k,:) - K_area*dA*dA_dP1;
                F_area_b(k,:) = F_area_b(k,:) - K_area*dA*dA_dP2;
                F_area_b(next,:) = F_area_b(next,:) - K_area*dA*dA_dP3;
                F_area_a(next,:) = F_area_a(next,:) - K_area*dA*dA_dP4;
            end

            % Gaussian curvature forces (simplified)
            F_gauss_a = zeros(n,2); F_gauss_b = zeros(n,2);
            for k = 1:n
                next = mod(k,n)+1;
                prev = mod(k-2,n)+1;
                F_gauss_a(k,:) = -(A*2*(-e_a(k,:))/(-l_b(k)*l_a(k)^2) + A*(-e_side(k,:))/(-l_side(k)^2*l_side(next))...
                    +A*2*(e_a(prev,:))/(-l_b(prev)*l_a(prev)^2) + A*(-e_side(k,:))/(-l_side(k)^2*l_side(prev)));
                F_gauss_b(k,:) = -(A*2*(-e_b(k,:))/(-l_a(k)*l_b(k)^2) + A*(e_side(k,:))/(-l_side(k)^2*l_side(next))...
                    +A*2*(e_b(prev,:))/(-l_a(prev)*l_b(prev)^2) + A*(e_side(k,:))/(-l_side(k)^2*l_side(prev)));
            end

            % Total forces and update positions
            vel_a = zeros(n,2); vel_b = zeros(n,2);
            for k = 1:n
                F_tot_a = F_area_a(k,:) + F_surface_a(k,:) - grad_side_a(k,:) + F_water(k,:) + F_bend_a(k,:) + F_gauss_a(k,:);
                F_tot_b = F_area_b(k,:) + F_surface_b(k,:) - grad_side_b(k,:) + F_ecm(k,:) + F_bend_b(k,:) + F_gauss_b(k,:);
                vel_a(k,:) = F_tot_a / mu;
                vel_b(k,:) = F_tot_b / mu;
            end
            for k = 1:n
                r_a(k,:) = r_a(k,:) + vel_a(k,:) * dt;
                r_b(k,:) = r_b(k,:) + vel_b(k,:) * dt;
            end
        end

        % --- Compute final total area ---
        x_b = r_b(:, 1);
        y_b = r_b(:, 2);
        next_idx = mod((1:n), n) + 1; % 获取下一个顶点的索引，最后一个顶点指向第一个
        total_area = 0.5 * abs(sum(x_b .* y_b(next_idx) - x_b(next_idx) .* y_b));

        Area_mat(i_g, i_p) = total_area;

        fprintf('gamma_b = %.4f, P_water = %.1f, Organoid Area = %.4e\n', gamma_b_val, P_water_val, total_area);
    end
end
%% ========== 4. Plot Phase Diagram ==========
figure;
imagesc(P_water_vals, gamma_b_vals, Area_mat);
set(gca, 'YDir', 'normal');
xlabel('Internal pressure P_{water}');
ylabel('Basal tension \gamma_b');
title('Phase diagram: Total area');
colorbar;
colormap(jet);
%% ========== 4. Plot Phase Diagram (3D with Grid Lines) ==========

% 1. Generate grid for the original data
[P_water_grid, gamma_b_grid] = meshgrid(P_water_vals, gamma_b_vals);

% 2. Generate a denser interpolation grid (100x100)
P_water_fine = linspace(min(P_water_vals), max(P_water_vals), 100);
gamma_b_fine = linspace(min(gamma_b_vals), max(gamma_b_vals), 100);
[P_water_fine_grid, gamma_b_fine_grid] = meshgrid(P_water_fine, gamma_b_fine);

% 3. Smooth the area data using cubic spline interpolation ('cubic')
Area_mat_smooth = interp2(P_water_grid, gamma_b_grid, Area_mat, P_water_fine_grid, gamma_b_fine_grid, 'cubic');


% ==========================================
% 2D Top-Down View (Fixed version retaining exponent labels)
% ==========================================
figure('Position', [150, 150, 950, 700], 'Name', '2D Phase Diagram');

contourf(P_water_fine_grid, gamma_b_fine_grid, Area_mat_smooth, 30, 'LineColor', 'none');
colormap(turbo);

% 1. Colorbar settings (kept unchanged)
c = colorbar;
c.Label.String = 'Total Organoid Area (m^2)';
c.Label.FontSize = 16;
c.Label.FontWeight = 'normal';
c.Label.FontName = 'Arial';
set(c, 'FontName', 'Arial', 'FontSize', 20);

% 2. Axis label settings
xlabel('Hydrostatic pressure P_{hyd} (Pa)', 'FontSize', 20, 'FontWeight', 'normal', 'FontName', 'Arial');
ylabel('Basal Tension \gamma_b (N/m)', 'FontSize', 20, 'FontWeight', 'normal', 'FontName', 'Arial');

% 3. Title settings (Key: shift title position upwards)
hTitle = title('Phase Diagram: Organoid Area', 'FontSize', 20, 'FontWeight', 'normal', 'FontName', 'Arial');

% 4. Axis tick and font settings
ax = gca;
set(ax, 'FontName', 'Arial', 'FontSize', 20, 'FontWeight', 'normal', 'LineWidth', 1.2);

% ==========================================
% Core Fix: Retain exponent labels + resolve overlapping/occlusion issues
% ==========================================
grid on;
axis tight;

% Fix 1: Retain the exponent label while preventing overlap with the title
ax.YAxis.Exponent = 0;  % Retain exponent label (allows automatic position adjustment)

% Fix 2: Increase the left margin (resolves Y-axis ticks being cut off)
% Get current axis position [left, bottom, width, height]
axPos = ax.Position;
% Increase left margin (from 0.13 to 0.18), while proportionally reducing the width
ax.Position = [0.18, axPos(2), axPos(3) - 0.05, axPos(4)];

% Fix 3: Reserve extra space for the exponent label
% Adjust the Y-axis label position to prevent the exponent label from being pushed into the title area
ax.YLabel.Position(1) = ax.YLabel.Position(1) * 1.2;  % Shift Y-axis label slightly to the right