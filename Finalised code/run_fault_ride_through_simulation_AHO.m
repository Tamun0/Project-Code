function run_fault_ride_through_simulation_AHO()
    % Main function to run a comparative simulation of fault ride-through
    % performance for Droop, VSM, and a detailed Andronov-Hopf VOC.

    clear; clc; close all;

    % ========================================================================
    % --- Simulation Case Selection ---
    % ========================================================================
    % Control Mode: 'Droop', 'VSM', or 'VOC'
    control_mode = 'Droop'; 
    
    % ========================================================================
    % 1. DEFINE SYSTEM PARAMETERS (Struct 'P')
    % ========================================================================
    P.f_nom = 50;
    P.w_nom = 2 * pi * P.f_nom;
    P.V_grid_rms = 230;
    P.V_grid_peak = P.V_grid_rms * sqrt(2);
    
    % --- Power Setpoint ---
    P.P_ref = 10000; % Operating at 10 kW power
    P.Q_ref = 0;     % Operating at 0 VAR
    
    % --- Grid Impedance (Weak Inductive) ---
    P.R_line = 0.2; 
    P.L_line = 5e-3;
    P.X_line = P.w_nom * P.L_line;
    
    % --- Controller Parameters ---
    % Droop
    P.mp_droop = (0.05 * P.w_nom) / P.P_ref; % 5% droop at rated power
    P.T_filter_droop = 0.05; 
    % VSM (Low-inertia model)
    P.J_vsm = 0.2;  P.D_vsm = 50.0; P.T_filter_vsm = 0.02;
    % VOC (Andronov-Hopf Oscillator)
    P.xi = 15;
    P.kv = P.V_grid_peak;
    P.ki = 3 * P.V_grid_peak / (P.P_ref + 1); % Scale ki with power
    P.C = 0.2679;
    P.phi = pi/2; % Configured for inductive grid
    
    % --- Fault Parameters ---
    P.t_fault_start = 1.0;
    P.fault_duration = 0.15;
    P.t_fault_clear = P.t_fault_start + P.fault_duration;
    P.V_fault_peak = 0.1 * P.V_grid_peak;
    
    % --- Simulation Setup ---
    t_end = 3.0;
    
    disp(['Running FAULT RIDE-THROUGH test for ' control_mode '...']);

    % ========================================================================
    % 2. SET UP AND RUN THE SIMULATION
    % ========================================================================
    options = odeset('RelTol', 1e-4, 'AbsTol', 1e-5);
    
    % The state vector and initial conditions depend on the controller
    if strcmp(control_mode, 'VOC')
        % State vector for VOC: x = [id, iq, v_alpha, v_beta]
        % Find initial steady state
        v_alpha0 = P.V_grid_peak; v_beta0 = 0;
        Z = P.R_line + 1j*P.X_line;
        I_complex = ( (v_alpha0 + 1j*v_beta0) - P.V_grid_peak) / Z;
        id0 = real(I_complex); iq0 = imag(I_complex);
        x0 = [id0; iq0; v_alpha0; v_beta0];
    else % Droop or VSM
        % State vector: x = [delta, w_dev, p_filt]
        Z_line = sqrt(P.R_line^2 + P.X_line^2);
        delta0 = asin(P.P_ref * Z_line / P.V_grid_peak);
        if isnan(delta0), delta0 = 0; end
        x0 = [delta0; 0; P.P_ref];
    end

    % Run simulation in three parts: pre-fault, during-fault, post-fault
    [t1, x1] = ode23t(@(t,x) fault_dynamics(t, x, P, control_mode, 0), [0 P.t_fault_start], x0, options);
    [t2, x2] = ode23t(@(t,x) fault_dynamics(t, x, P, control_mode, 1), [P.t_fault_start P.t_fault_clear], x1(end,:)', options);
    [t3, x3] = ode23t(@(t,x) fault_dynamics(t, x, P, control_mode, 0), [P.t_fault_clear t_end], x2(end,:)', options);
    
    t = [t1; t2(2:end); t3(2:end)];
    x = [x1; x2(2:end,:); x3(2:end,:)];

    % ========================================================================
    % 3. POST-PROCESS AND PLOT RESULTS
    % ========================================================================
    disp('Simulation finished. Processing results...');
    
    if strcmp(control_mode, 'VOC')
        id = x(:,1); iq = x(:,2); v_alpha = x(:,3); v_beta = x(:,4);
        V_grid_vec = (t < P.t_fault_start | t > P.t_fault_clear) * P.V_grid_peak + ...
                     (t >= P.t_fault_start & t <= P.t_fault_clear) * P.V_fault_peak;
        P_elec = V_grid_vec .* id; % Simplified power calculation in dq frame
        
        % Calculate frequency from oscillator states
        dv_alpha_dt = gradient(v_alpha) ./ gradient(t);
        dv_beta_dt = gradient(v_beta) ./ gradient(t);
        w_inst = (v_alpha .* dv_beta_dt - v_beta .* dv_alpha_dt) ./ (v_alpha.^2 + v_beta.^2);
        freq_hz = w_inst / (2*pi);
        
        % Angle is angle of oscillator vector
        angle_deg = rad2deg(unwrap(atan2(v_beta, v_alpha)));
    else % Droop or VSM
        delta_rad = x(:,1); w_dev = x(:,2);
        V_grid_vec = (t < P.t_fault_start | t > P.t_fault_clear) * P.V_grid_peak + ...
                     (t >= P.t_fault_start & t <= P.t_fault_clear) * P.V_fault_peak;
        P_elec = (V_grid_vec / Z_line) .* sin(delta_rad);
        freq_hz = P.f_nom + w_dev / (2*pi);
        angle_deg = rad2deg(unwrap(delta_rad));
    end
    
    figure('Name', [control_mode ' Fault Ride-Through Response']);
    
    subplot(3,1,1); plot(t, freq_hz, 'b', 'LineWidth', 2);
    hold on; line([P.t_fault_start P.t_fault_start], ylim, 'Color', 'r', 'LineStyle', '--');
    line([P.t_fault_clear P.t_fault_clear], ylim, 'Color', 'k', 'LineStyle', '--');
    title('Inverter Frequency'); xlabel('Time (s)'); ylabel('Frequency (Hz)'); grid on;

    subplot(3,1,2); plot(t, P_elec/1000, 'm', 'LineWidth', 2);
    hold on; line([P.t_fault_start P.t_fault_start], ylim, 'Color', 'r', 'LineStyle', '--');
    line([P.t_fault_clear P.t_fault_clear], ylim, 'Color', 'k', 'LineStyle', '--');
    title('Inverter Active Power Output'); xlabel('Time (s)'); ylabel('Power (kW)'); grid on;

    subplot(3,1,3); plot(t, angle_deg, 'g', 'LineWidth', 2);
    hold on; line([P.t_fault_start P.t_fault_start], ylim, 'Color', 'r', 'LineStyle', '--');
    line([P.t_fault_clear P.t_fault_clear], ylim, 'Color', 'k', 'LineStyle', '--');
    title('Inverter Internal Angle'); xlabel('Time (s)'); ylabel('Angle (degrees)'); grid on;
end

% ============================================================================
% THE DYNAMICS FUNCTION FOR FAULT SIMULATION
% ============================================================================
function dxdt = fault_dynamics(t, x, P, control_mode, fault_status)
    % Determine grid voltage based on fault status
    if fault_status == 1, V_grid_dynamic = P.V_fault_peak;
    else, V_grid_dynamic = P.V_grid_peak; end
    
    % --- Controller Dynamics ---
    if strcmp(control_mode, 'VOC')
        % Unpack VOC states: [id, iq, v_alpha, v_beta]
        id = x(1); iq = x(2); v_alpha = x(3); v_beta = x(4);
        
        % --- Interface between stationary (ab) and rotating (dq) frames ---
        % The grid is the reference frame, its angle is w_nom*t
        theta_grid = P.w_nom * t;
        
        % Inverter voltage in the grid's rotating frame
        v_d_inv = v_alpha * cos(theta_grid) + v_beta * sin(theta_grid);
        v_q_inv = -v_alpha * sin(theta_grid) + v_beta * cos(theta_grid);
        
        % --- Line Dynamics (State Equations for id, iq) ---
        did_dt = (1/P.L_line) * (-P.R_line*id + P.w_nom*P.L_line*iq + v_d_inv - V_grid_dynamic);
        diq_dt = (1/P.L_line) * (-P.R_line*iq - P.w_nom*P.L_line*id + v_q_inv - 0); % Vq_grid = 0
        
        % --- VOC Feedback Calculation ---
        % Current in the stationary frame
        i_alpha = id * cos(theta_grid) - iq * sin(theta_grid);
        i_beta = id * sin(theta_grid) + iq * cos(theta_grid);
        
        % Current setpoints from power setpoints
        V_sq_inv = v_alpha^2 + v_beta^2;
        if V_sq_inv < 1e-3, V_sq_inv = 1e-3; end % Avoid division by zero
        i_alpha_star = (2/3) * (v_alpha * P.P_ref + v_beta * P.Q_ref) / V_sq_inv;
        i_beta_star = (2/3) * (v_beta * P.P_ref - v_alpha * P.Q_ref) / V_sq_inv;
        
        % Current error with rotation for grid type
        R_phi = [cos(P.phi) -sin(P.phi); sin(P.phi) cos(P.phi)];
        i_err = [i_alpha - i_alpha_star; i_beta - i_beta_star];
        u = P.ki * R_phi * i_err;
        
        % --- Andronov-Hopf Oscillator Dynamics (State Equations for v_alpha, v_beta) ---
        V_nom_sq = P.kv^2;
        dv_alpha_dt = (P.xi/P.kv^2)*(2*V_nom_sq - V_sq_inv)*v_alpha - P.w_nom*v_beta - (P.kv/P.C)*u(1);
        dv_beta_dt  = (P.xi/P.kv^2)*(2*V_nom_sq - V_sq_inv)*v_beta  + P.w_nom*v_alpha - (P.kv/P.C)*u(2);
        
        dxdt = [did_dt; diq_dt; dv_alpha_dt; dv_beta_dt];
        
    else % Droop or VSM (behavioral model)
        % Unpack states: [delta, w_dev, p_filt]
        delta = x(1); w_dev = x(2); p_filt = x(3);
        
        Z_line = sqrt(P.R_line^2 + P.X_line^2);
        P_elec = (V_grid_dynamic / Z_line) * sin(delta);
        
        if strcmp(control_mode, 'Droop')
            dp_filt_dt = (1/P.T_filter_droop) * (P_elec - p_filt);
            dw_dev_dt = (1/P.J_droop) * (P.P_ref - p_filt - P.D_droop * w_dev);
        else % VSM
            dp_filt_dt = (1/P.T_filter_vsm) * (P_elec - p_filt);
            dw_dev_dt = (1/P.J_vsm) * (P.P_ref - p_filt - P.D_vsm * w_dev);
        end
        ddelta_dt = w_dev;
        dxdt = [ddelta_dt; dw_dev_dt; dp_filt_dt];
    end
end
