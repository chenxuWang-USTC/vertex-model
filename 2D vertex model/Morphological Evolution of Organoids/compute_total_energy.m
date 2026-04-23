function [E, grad] = compute_total_energy(x, n, params)
% COMPUTE_TOTAL_ENERGY Calculates the total free energy of the organoid system.
%
% Inputs:
%   x      : State vector [r_a(:); r_b(:)], size (4n x 1) or (2n x 2) reshaped.
%   n      : Number of trapezoidal units (cells).
%   params : Structure containing physical parameters and local tension vectors.
%            Required fields: gamma_a, gamma_b, gamma_l, K_area, A0, A, kappa_bend.
%            Optional fields for spatial variation: gamma_a_local, gamma_b_local, gamma_l_local.
%
% Outputs:
%   E      : Total scalar energy (J).
%   grad   : (Optional) Gradient of energy w.r.t coordinates (not fully implemented in this version).

    %% 1. Unpack Parameters & Handle Spatial Variations
    % Check if local tension vectors exist in params; otherwise use global scalars
    if isfield(params, 'gamma_a_local') && length(params.gamma_a_local) == n
        gamma_a_vec = params.gamma_a_local;
    else
        gamma_a_vec = ones(n, 1) * params.gamma_a;
    end
    
    if isfield(params, 'gamma_b_local') && length(params.gamma_b_local) == n
        gamma_b_vec = params.gamma_b_local;
    else
        gamma_b_vec = ones(n, 1) * params.gamma_b;
    end
    
    if isfield(params, 'gamma_l_local') && length(params.gamma_l_local) == n
        gamma_l_vec = params.gamma_l_local;
    else
        gamma_l_vec = ones(n, 1) * params.gamma_l;
    end

    K_area = params.K_area;
    A0     = params.A0;
    A_gauss= params.A;       % Gaussian curvature parameter
    kappa_bend = params.kappa_bend;

    %% 2. Reshape Coordinates and Initialize Geometric Arrays
    % Reshape input vector into coordinate matrices (n x 2)
    r_a = reshape(x(1:2*n), n, 2);
    r_b = reshape(x(2*n+1:end), n, 2);

    % Pre-allocate geometric arrays
    l_a = zeros(n, 1);
    l_b = zeros(n, 1);
    l_side = zeros(n, 1);
    e_a = zeros(n, 2);
    e_b = zeros(n, 2);
    e_side = zeros(n, 2);
    
    % Initialize Energy Components
    E_surface = 0;
    E_gauss   = 0;
    E_area    = 0;
    E_bend    = 0;

    %% 3. Compute Geometry (Edge Lengths and Unit Vectors)
    for i = 1:n
        next_idx = mod(i, n) + 1;
        
        % Edge vectors
        v_a_i    = r_a(next_idx, :) - r_a(i, :);
        v_b_i    = r_b(next_idx, :) - r_b(i, :);
        v_side_i = r_b(i, :) - r_a(i, :);
        
        % Lengths
        l_a(i)    = norm(v_a_i);
        l_b(i)    = norm(v_b_i);
        l_side(i) = norm(v_side_i);
        
        % Unit vectors (with protection against division by zero)
        if l_a(i) > 1e-12,    e_a(i, :)    = v_a_i / l_a(i);    end
        if l_b(i) > 1e-12,    e_b(i, :)    = v_b_i / l_b(i);    end
        if l_side(i) > 1e-12, e_side(i, :) = v_side_i / l_side(i); end
    end

    %% 4. Calculate Surface Energy & Gaussian Curvature Energy
    % Note: Using local tension vectors (gamma_X_vec) for spatial heterogeneity
    for i = 1:n
        next_idx = mod(i, n) + 1;
        
        % --- Apical (Inner) Surface Energy ---
        % Formula: 0.5 * sqrt(3) * gamma * L^2 (Based on hexagonal/trapezoidal approximation)
        E_surface = E_surface + 0.5 * sqrt(3) * gamma_a_vec(i) * (l_a(i)^2);
        
        % --- Basal (Outer) Surface Energy ---
        E_surface = E_surface + 0.5 * sqrt(3) * gamma_b_vec(i) * (l_b(i)^2);
        
        % --- Lateral (Side) Surface Energy ---
        % Coupling adjacent side lengths and ring lengths
        E_surface = E_surface + gamma_l_vec(i) * (l_a(i) + l_b(i)) * (l_side(i) + l_side(next_idx)) / (2 * sqrt(3));
        
        % --- Gaussian Curvature Constraint Energy ---
        % Penalizes deviation from ideal packing geometry
        % Term 1: Ring coupling
        % Term 2: Side coupling
        if (l_a(i)*l_b(i)) > 1e-24 && (l_side(i)*l_side(next_idx)) > 1e-24
             E_gauss = E_gauss + (2 * A_gauss) / (l_a(i) * l_b(i)) ...
                             + A_gauss / (l_side(i) * l_side(next_idx));
        end
    end

    %% 5. Calculate Area Constraint Energy
    % Quadrilateral area constraint for each cell
    for i = 1:n
        next_idx = mod(i, n) + 1;
        
        P1 = r_a(i, :);
        P2 = r_b(i, :);
        P3 = r_b(next_idx, :);
        P4 = r_a(next_idx, :);
        
        % Calculate current area using helper function
        A_curr = quad_area(P1, P2, P3, P4);
        
        % Harmonic potential: 0.5 * K * (A - A0)^2
        dA = A_curr - A0;
        E_area = E_area + 0.5 * K_area * (dA^2);
    end

    %% 6. Calculate Bending Energy
    % Penalizes deviation of internal angles from ideal hexagonal angle (120 degrees / cos = -0.5)
    % Note: Original code used sqrt(3)/2 (~0.866, which is 30 degrees). 
    % Assuming original logic intends to maintain specific vertex angles.
    target_cos = sqrt(3)/2; 
    
    % Helper function to compute bending energy for a ring
    % We iterate through nodes to check angles between consecutive edges
    
    % --- Inner Ring Bending ---
    for i = 1:n
        prev_idx = mod(i-2, n) + 1; % MATLAB indexing adjustment for previous node
        next_idx = mod(i, n) + 1;
        
        r_prev = r_a(prev_idx, :);
        r_curr = r_a(i, :);
        r_next = r_a(next_idx, :);
        
        v1 = r_prev - r_curr;
        v2 = r_next - r_curr;
        
        len1 = norm(v1);
        len2 = norm(v2);
        
        if len1 > 1e-12 && len2 > 1e-12
            cos_theta = dot(v1, v2) / (len1 * len2);
            cos_theta = max(-1, min(1, cos_theta)); % Clamp for numerical stability
            
            % Only penalize if deviation is significant (as per original logic)
            % Note: Standard bending usually penalizes any deviation. 
            % Keeping original "if" condition but ensuring robustness.
            if abs(cos_theta - target_cos) > 1e-6
                 E_bend = E_bend + kappa_bend * (cos_theta - target_cos)^2;
            end
        end
    end
    
    % --- Outer Ring Bending ---
    for i = 1:n
        prev_idx = mod(i-2, n) + 1;
        next_idx = mod(i, n) + 1;
        
        r_prev = r_b(prev_idx, :);
        r_curr = r_b(i, :);
        r_next = r_b(next_idx, :);
        
        v1 = r_prev - r_curr;
        v2 = r_next - r_curr;
        
        len1 = norm(v1);
        len2 = norm(v2);
        
        if len1 > 1e-12 && len2 > 1e-12
            cos_theta = dot(v1, v2) / (len1 * len2);
            cos_theta = max(-1, min(1, cos_theta));
            
            if abs(cos_theta - target_cos) > 1e-6
                 E_bend = E_bend + kappa_bend * (cos_theta - target_cos)^2;
            end
        end
    end

    %% 7. Sum Total Energy
    E = E_surface + E_gauss + E_area + E_bend;

end