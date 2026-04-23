function L = intital()

    %% ========== Parameter Settings ==========
    n = 40;                      % Number of trapezoidal units
    R_a = 66e-6;                 % Initial inner ring radius
    R_b = 77e-6;                 % Initial outer ring radius
    dt = 0.1;                    % Time step
    num_steps = 20000;           % Simulation steps
    mu = 0.4;                    % Viscous coefficient 
    gamma_a = 2e-3;              % Inner ring intrinsic tension 
    gamma_b = -1e-3;             % Outer ring intrinsic tension 
    gamma_l = -3e-3;             % Lateral side intrinsic tension 
    V0 = 1.5e-15;
    A0 = 3e-10;
    A  = 5e-23;
    K_area = 1e9;
    kappa_bend = 1e-12;
    
    %% ========== 1. Initial Geometry ==========
    theta = (0:n-1) * (2*pi/n);  % Avoid cumulative error from linspace
    r_a = [R_a * cos(theta); R_a * sin(theta)]';  % Inner ring vertices (n x 2)
    r_b = [R_b * cos(theta); R_b * sin(theta)]';  % Outer ring vertices (n x 2)
    
    %% ========== 2. Initialize Edge Vectors ==========
    v_a = zeros(n, 2); v_b = zeros(n, 2); v_side = zeros(n, 2);
    for i = 1:n
        next_idx = mod(i, n) + 1;
        v_a(i,:) = r_a(next_idx,:) - r_a(i,:);
        v_b(i,:) = r_b(next_idx,:) - r_b(i,:);
        v_side(i,:) = r_b(i,:) - r_a(i,:); % Side vector points from apical to basal
    end

    %% ========== 3. Pre-calculate Initial Quantities (Avoid Redundant Calculation) ==========
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

 % ========== 4. Main Loop ==========
for step = 1:num_steps
        
        %% ========== 4.1 Update Edge Vectors ==========
        for i = 1:n
            next_idx = mod(i, n) + 1;
            v_a(i,:) = r_a(next_idx,:) - r_a(i,:);
            v_b(i,:) = r_b(next_idx,:) - r_b(i,:);
            v_side(i,:) = r_b(i,:) - r_a(i,:); 
            % P_ecm(i)=80;
            % P_water(i)=100;
            P_ecm(i) = 500/(77e-6)*norm(r_b(i,:)) - 400;
            P_water(i) = 800*sqrt(66e-6)/sqrt(norm(r_a(i,:))) - 700;
        end
       
        %% ========== 4.2 Update All Geometric Quantities at Once ==========
        for i = 1:n
            l_a(i) = norm(v_a(i,:));
            l_b(i) = norm(v_b(i,:));
            l_side(i) = norm(v_side(i,:));
            L_side(i) = l_side(i);
            e_side(i,:) = v_side(i,:) ./ l_side(i);
        end
        
      %% ========== Bending Stiffness Force (Efficient Version: Iterate units, automatically cover contributions from i and i-1) ==========
F_bend_a = zeros(n, 2);
F_bend_b = zeros(n, 2);

target_cos = 0;  % Target: Perpendicular (90 degrees)

for i = 1:n
    next_idx = mod(i, n) + 1; % i+1
    
    % Get coordinates of four key points involved in the current unit
    r_ai = r_a(i, :);
    r_an = r_a(next_idx, :);
    r_bi = r_b(i, :);
    % r_bn = r_b(next_idx, :); % Outer ring angle mainly uses r_bi and r_bn, but side shares r_ai
    
    %% --- 1. Inner Ring Angle (Formed by side edge of unit i and forward inner ring edge) ---
    % Vector definitions:
    
    u = r_bi - r_ai; 
    v = r_an - r_ai;
    
    len_u = norm(u);
    len_v = norm(v);
    
    if len_u > 1e-12 && len_v > 1e-12
        e_u = u / len_u;
        e_v = v / len_v;
        
        cos_theta = dot(e_u, e_v);
        cos_theta = max(-1, min(1, cos_theta));
        
        if abs(cos_theta - target_cos) > 1e-6
            % F = -dE/dr = -2*k*(cos - target)*d(cos)/dr
            factor = -2 * kappa_bend * (cos_theta - target_cos);
            
            % Gradients
            grad_u = (e_v - cos_theta * e_u) / len_u; % d(cos)/du
            grad_v = (e_u - cos_theta * e_v) / len_v; % d(cos)/dv
            
            % === Key: Distribute forces to all nodes involved in this unit ===
            % u = r_bi - r_ai  => Apply +factor*grad_u to r_bi, -factor*grad_u to r_ai
            % v = r_an - r_ai  => Apply +factor*grad_v to r_an, -factor*grad_v to r_ai
            
            F_bend_b(i, :)     = F_bend_b(i, :)     + factor * grad_u; % Acts on outer ring node i
            F_bend_a(i, :)     = F_bend_a(i, :)     - factor * grad_u; % Acts on inner ring node i (from side)
            
            F_bend_a(next_idx, :) = F_bend_a(next_idx, :) + factor * grad_v; % Acts on inner ring node i+1
            F_bend_a(i, :)        = F_bend_a(i, :)        - factor * grad_v; % Acts on inner ring node i (from tangential)
        end
    end
    
    %% --- 2. Outer Ring Angle (Formed by side edge of unit i and forward outer ring edge) ---
    % Vector definitions:
    % Side w = r_ai - r_bi (Outer -> Inner, to form angle with outer tangential)
    % Tangential z = r_bn - r_bi (i -> i+1)
    
    r_bn = r_b(next_idx, :);
    w = r_ai - r_bi; 
    z = r_bn - r_bi;
    
    len_w = norm(w);
    len_z = norm(z);
    
    if len_w > 1e-12 && len_z > 1e-12
        e_w = w / len_w;
        e_z = z / len_z;
        
        cos_theta_b = dot(e_w, e_z);
        cos_theta_b = max(-1, min(1, cos_theta_b));
        
        if abs(cos_theta_b - target_cos) > 1e-6
            factor = -2 * kappa_bend * (cos_theta_b - target_cos);
            
            grad_w = (e_z - cos_theta_b * e_w) / len_w;
            grad_z = (e_w - cos_theta_b * e_z) / len_z;
            
            % === Key: Distribute forces to all nodes involved in this unit ===         
            F_bend_a(i, :)     = F_bend_a(i, :)     + factor * grad_w; % Acts on inner ring node i
            F_bend_b(i, :)     = F_bend_b(i, :)     - factor * grad_w; % Acts on outer ring node i (from side)
            
            F_bend_b(next_idx, :) = F_bend_b(next_idx, :) + factor * grad_z; % Acts on outer ring node i+1
            F_bend_b(i, :)        = F_bend_b(i, :)        - factor * grad_z; % Acts on outer ring node i (from tangential)
        end
    end
end

        %% ========== 4.3 Calculate Surface Energy Gradient ==========
        grad_U_a = zeros(n, 2); 
        grad_U_b = zeros(n, 2);
        grad_side_a = zeros(n, 2); 
        grad_side_b = zeros(n, 2);

        for i = 1:n
            next_idx = mod(i, n) + 1;

            % Inner ring surface tension
            grad_U_a(i,:) = grad_U_a(i,:) + sqrt(3) * gamma_a * l_a(i) * (-e_a(i,:));
            grad_U_a(next_idx,:) = grad_U_a(next_idx,:) + sqrt(3) * gamma_a * l_a(i) * e_a(i,:);
            
            % Outer ring surface tension
            grad_U_b(i,:) = grad_U_b(i,:) + sqrt(3) * gamma_b * l_b(i) * (-e_b(i,:));
            grad_U_b(next_idx,:) = grad_U_b(next_idx,:) + sqrt(3) * gamma_b * l_b(i) * e_b(i,:);

             % Lateral membrane surface tension 
         
            % Gradient w.r.t r_a(i) (via l_a(i) and l_side(i))
            grad_side_a(i,:) = grad_side_a(i,:) ...
            - e_a(i,:) * gamma_l / sqrt(12)* (l_side(i) + l_side(next_idx)) ...  % via l_a(i)
            - e_side(i,:) *gamma_l / sqrt(12) * (l_a(i) + l_b(i));             % via l_side(i)
        
            % Gradient w.r.t r_a(next) (via l_a(i))
            grad_side_a(next_idx,:) = grad_side_a(next_idx,:) ...
            + e_a(i,:) * gamma_l / sqrt(12) * (l_side(i) + l_side(next_idx)) ... % via l_a(i)
            - e_side(next_idx,:) *gamma_l / sqrt(12) * (l_a(i) + l_b(i));       % via l_side(next)
    
            % Gradient w.r.t r_b(i) (via l_b(i) and l_side(i))
            grad_side_b(i,:) = grad_side_b(i,:) ...
            - e_b(i,:) * gamma_l / sqrt(12)* (l_side(i) + l_side(next_idx)) ...  % via l_b(i)
            + e_side(i,:) * gamma_l / sqrt(12) * (l_a(i) + l_b(i));            % via l_side(i)
    
            % Gradient w.r.t r_b(next) (via l_b(i))
            grad_side_b(next_idx,:) = grad_side_b(next_idx,:) ...
            + e_b(i,:) * gamma_l / sqrt(12) * (l_side(i) + l_side(next_idx)) ...
            + e_side(next_idx,:) * gamma_l / sqrt(12) * (l_a(i) + l_b(i)); 
        end
       
        F_surface_a = -grad_U_a;
        F_surface_b = -grad_U_b;
        
        %% ========== 4.4 Pressure Calculation ==========
        F_water = zeros(n, 2); 
        F_ecm = zeros(n, 2);
        for i = 1:n
            prev_idx = mod(i-2, n) + 1;  % Index of previous edge
            
            % === (A) Internal pressure (pushes outward) ===
            v_prev = v_a(prev_idx,:);  % Use inner ring edge to calculate normal
            v_i = v_a(i,:);
            t1 = v_i / norm(v_i);
            t2 = v_prev / norm(v_prev);
            n1 = [t1(2), -t1(1)]; 
            e1 = n1 / norm(n1);  % Normal outward
            n2 = [t2(2), -t2(1)]; 
            e2 = n2 / norm(n2);
            F_water(i,:) = (P_water(i) ) *  sqrt(3) / 2 * (l_a(prev_idx)^2 * e2 + l_a(i)^2 * e1);
           
            % === (B) External pressure (pushes inward) ===
            v_prev = v_b(prev_idx,:);  
            v_i = v_b(i,:);
            t1 = v_i / norm(v_i);
            t2 = v_prev / norm(v_prev);
            n1 = [-t1(2), t1(1)]; 
            e1 = n1 / norm(n1);  % Normal inward
            n2 = [-t2(2), t2(1)]; 
            e2 = n2 / norm(n2);
            F_ecm(i,:) = (P_ecm(i) ) *  sqrt(3) / 2 * (l_b(prev_idx)^2 * e2 + l_b(i)^2 * e1);
           
        end
        
          %% ========== 4.6 Area Constraint Force ==========
        % Initialize forces (zero out)
        F_area_a = zeros(n, 2);
        F_area_b = zeros(n, 2);

        for i = 1:n
            next_idx = mod(i, n) + 1;
    
            P1 = r_a(i, :);
            P2 = r_b(i, :);
            P3 = r_b(next_idx, :);
            P4 = r_a(next_idx, :);
    
            % Current quadrilateral area
            A_current(i) = quad_area(P1, P2, P3, P4);
    
            % Area error
            dA = A_current(i) - A0;
    
    
            % Analytical gradients of area w.r.t vertices (based on directed area formula)
            dA_dP1 = -0.5 * [ -(P2(2) - P4(2)),  P2(1) - P4(1) ];
            dA_dP2 = -0.5 * [ -(P3(2) - P1(2)),  P3(1) - P1(1) ];
            dA_dP3 = -0.5 * [ -(P4(2) - P2(2)),  P4(1) - P2(1) ];
            dA_dP4 = -0.5 * [ -(P1(2) - P3(2)),  P1(1) - P3(1) ];
    
            % Note: Area energy is (1/2) K_area (A - A0)^2,
            % So Force = -grad(E) = -K_area (A - A0) grad(A)
            % i.e.: F = -K_area * dA * dA_dP
    
            % Accumulate to corresponding vertices (Key: += not =)
            F_area_a(i, :)     = F_area_a(i, :)     - K_area * dA * dA_dP1;
            F_area_b(i, :)     = F_area_b(i, :)     - K_area * dA * dA_dP2;
            F_area_b(next_idx, :)  = F_area_b(next_idx, :)  - K_area * dA * dA_dP3;
            F_area_a(next_idx, :)  = F_area_a(next_idx, :)  - K_area * dA * dA_dP4;
        end
        
        %====== 4.6 Gaussian Free Energy Constraint =========
        for i=1:1:n
            next_idx = mod(i, n) + 1;
            prev_idx = mod(i-2, n) + 1;
            F_guass_a(i,:) = -(A * 2 * (-e_a(i,:)) / (-l_b(i) * l_a(i) ^ 2) + A * (-e_side(i,:)) / (-l_side(i) ^ 2 * l_side(next_idx))...
            +A * 2 * (e_a(prev_idx,:)) / (-l_b(prev_idx) * l_a(prev_idx) ^ 2) + A * (-e_side(i,:)) / (-l_side(i) ^ 2 * l_side(prev_idx)));
            F_guass_b(i,:) = -(A * 2 * (-e_b(i,:)) / (-l_a(i) * l_b(i) ^ 2) + A * (e_side(i,:)) / (-l_side(i) ^ 2 * l_side(next_idx))...
            +A * 2 * (e_b(prev_idx,:)) / (-l_a(prev_idx) * l_b(prev_idx) ^ 2) + A * (e_side(i,:)) / (-l_side(i) ^ 2 * l_side(prev_idx)));
    
        end

        for i=1:1:n
        F_total_a(i,:) = F_surface_a(i,:) - grad_side_a(i,:)  + F_water(i,:) + F_area_a(i,:) + F_guass_a(i,:) + F_bend_a(i,:);%+F_vol_a(i,:)  
        F_a(i) = norm(F_total_a(i,:));
        F_total_b(i,:) = F_surface_b(i,:) - grad_side_b(i,:)  + F_area_b(i,:) + F_ecm(i,:) + F_guass_b(i,:) + F_bend_b(i,:); %+ F_vol_b(i,:)
        F_b(i) = norm(F_total_b(i,:));
   
        % ========== 4.7 Dynamics ==========
        vel_a(i,:) = F_total_a(i,:) / mu;
        V_a(i) = norm(vel_a);
        vel_b(i,:) = F_total_b(i,:) / mu;
        V_b(i) = norm(vel_b);
        end

   
        % ========== 4.8 Update Positions ==========
        for i=1:1:n
            
            r_a(i,:) = r_a(i,:) + vel_a(i,:) * dt;
    
            r_b(i,:) = r_b(i,:) + vel_b(i,:) * dt;
   
        end

    % ========== 4.9 Check Convergence: Inner and Outer ring velocities below threshold ==========
    max_vel_a = max(sqrt(sum(vel_a.^2, 2)));  % Max inner ring velocity (m/s)
    max_vel_b = max(sqrt(sum(vel_b.^2, 2)));  % Max outer ring velocity (m/s)
    
    vel_tol = 1e-10;  % Velocity convergence threshold (adjustable)
    
    if (max_vel_a < vel_tol) && (max_vel_b < vel_tol)
        % % % %  % Update graphics and save last frame
        % % % % set(h_inner, 'XData', [r_a(:,1); r_a(1,1)], 'YData', [r_a(:,2); r_a(1,2)]);
        % % % % set(h_outer, 'XData', [r_b(:,1); r_b(1,1)], 'YData', [r_b(:,2); r_b(1,2)]);
        % % % % for i = 1:n
        % % % %     set(h_side(i), 'XData', [r_a(i,1), r_b(i,1)], 'YData', [r_a(i,2), r_b(i,2)]);
        % % % % end
        % % % % writeVideo(video, getframe(gcf));   
        break;  % Terminate main loop early
    end  

end

    L.r_a = r_a;
    L.r_b = r_b;
    L.A0 = A_current(1);
   
 end