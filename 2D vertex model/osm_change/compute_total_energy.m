function [E, grad] = compute_total_energy(x, n, params)
% x: [r_a(:); r_b(:)] -> (4n x 1) vector
% n: Current number of units
% params: Structure containing all physical parameters
 
    % Unpack parameters
    gamma_a = params.gamma_a;
    gamma_b = params.gamma_b;
    gamma_l = params.gamma_l;
    K_area = params.K_area;
    A0 = params.A0;
    A = params.A;

    % ===== Initialize energy-related variables ========
    
    A_curr = zeros(n);
    % ===== Initialize gradient variables (Required to prevent crash when nargout > 1) =====
    if nargout > 1
        grad_U_a = zeros(n, 2);
        grad_U_b = zeros(n, 2);
        grad_side_a = zeros(n, 2);
        grad_side_b = zeros(n, 2);
        grad_a_area = zeros(n, 2);
        grad_b_area = zeros(n, 2);
        F_guass_a = zeros(n, 2);
        F_guass_b = zeros(n, 2);
        grad_a = zeros(n, 2);
        grad_b = zeros(n, 2);
        
    end
    % Unpack coordinates
    r_a = reshape(x(1:2*n), n, 2);
    r_b = reshape(x(2*n+1:end), n, 2);

    % Initialize energy components
    E_surface = 0;
    E_osm = 0;
    E_ecm = 0;
    E_area = 0;
    E_guass = 0;
    
    % Initialize vectors and lengths
    for i = 1:n
        next_idx = mod(i, n) + 1;
        prev_idx = mod(i-2, n) + 1;
        % Edge vectors
        v_a(i,:) = r_a(next_idx,:) - r_a(i,:);
        v_b(i,:) = r_b(next_idx,:) - r_b(i,:);
        v_side(i,:) = r_b(i,:) - r_a(i,:);
        
        l_a(i) = norm(v_a(i,:));
        l_b(i) = norm(v_b(i,:));
        l_side(i) = norm(v_side(i,:));
        
        % Unit vectors
        e_a(i,:) = v_a(i,:) / l_a(i);
        e_b(i,:) = v_b(i,:) / l_b(i);
        e_side(i,:) = v_side(i,:) / l_side(i); 
        
    end
    
    %% ============ 1. Surface Energy + Gradient ===============
    for i = 1:n
         next_idx = mod(i, n) + 1;
         prev_idx = mod(i-2, n) + 1;
        % Inner ring tension energy
        E_surface = E_surface + gamma_a * l_a(i) ^ 2 * 0.5 * sqrt(3);
        % Outer ring tension energy
        E_surface = E_surface + gamma_b * l_b(i) ^ 2 * 0.5 * sqrt(3);
        % Lateral membrane tension energy
        E_surface = E_surface + gamma_l * ( l_b(i) + l_a(i) ) * (l_side(i) + l_side(next_idx)) / (2 * sqrt(3));    
        % Gaussian constraint energy 
        E_guass = E_guass + 2 * A /(l_a(i)*l_b(i)) + A / (l_side(i) * l_side(next_idx));
       
    end

    %% ============== 3. Area Constraint Energy ==================
    for i = 1:n
        next_idx = mod(i, n) + 1;
        P1 = r_a(i,:); 
        P2 = r_b(i,:); 
        P3 = r_b(next_idx,:); 
        P4 = r_a(next_idx,:);
        A_curr(i) = quad_area(P1, P2, P3, P4);
        dA = A_curr(i) - A0;
        E_area = E_area + 0.5 * K_area * dA ^ 2;
    end
 
    %% =============== 4. Bending Energy ====================
   E_bend = 0;
   kappa_bend = params.kappa_bend;
    for i = 1:n
    next_idx = mod(i, n) + 1;
    prev_idx = mod(i-2, n) + 1;
   
    r_prev = r_a(prev_idx, :);
    r_curr = r_a(i, :);
    r_next = r_a(next_idx, :);
    
    v1 = r_prev - r_curr;
    v2 = r_next - r_curr;
    cos_theta = dot(v1, v2) / (norm(v1) * norm(v2));
    cos_theta = max(-1, min(1, cos_theta));
    
    if abs(cos_theta) < sqrt(3)/2
        bend_E_mag = kappa_bend * (cos_theta - sqrt(3)/2)^2;
        E_bend = E_bend + bend_E_mag;
    end  


    r_prev = r_b(prev_idx, :);
    r_curr = r_b(i, :);
    r_next = r_b(next_idx, :);
    
    v1 = r_prev - r_curr;
    v2 = r_next - r_curr;
    cos_theta = dot(v1, v2) / (norm(v1) * norm(v2));
    cos_theta = max(-1, min(1, cos_theta));
    
    if abs(cos_theta) < sqrt(3)/2
        bend_E_mag = kappa_bend * (cos_theta - sqrt(3)/2)^2;
        E_bend = E_bend + bend_E_mag;
    end  
      
    end

    %% Total Energy
    E = E_surface + E_guass + E_area + E_bend; %;
    

end