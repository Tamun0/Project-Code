function run_fault_ride_through_simulation_uVOC()
    % Main function to run a comparative simulation of fault ride-through
    % performance. This version includes a full Unified Virtual Oscillator
    % Control (uVOC) model with fault-specific dynamics.

    clear; clc; close all;

    % ========================================================================
    % --- Simulation Case Selection ---
    % ========================================================================
    % Control Mode: 'Droop_Pure', 'VSM', 'VOC_AHO', or 'uVOC'
    control_mode = 'uVOC'; 
    
    % ========================================================================
    % 1. DEFINE SYSTEM PARAMETERS (Struct 'P')
    % ========================================================================
    P.f_nom = 50;
    P.w_nom = 2 * pi * P.f_nom;
    P.V_grid_nominal = 1.0; % Per-unit voltage of infinite bus
    P.V_peak = 1.0; % Using per-unit, peak voltage is 1.0
    
    % --- Power Setpoint ---
    P.P_ref = 0.8; 
    P.Q_ref = 0;
    
    % --- Controller Parameters ---
    P.mp_droop = (0.05 * P.w_nom) / P.P_ref; P.T_filter_droop = 0.05; 
    P.J_vsm = 5.0;  P.D_vsm = 50.0; P.T_filter_vsm = 0.02;
    % AHO-VOC and uVOC model parameters
    P.xi = 25; P.kv = P.V_peak; P.ki = 3 * P.V_peak / (P.P_ref + 1);
    P.C_voc = 0.1; P.phi = pi/2;
    P.gamma = 5.0; % Gain for voltage synchronization in uVOC

    % --- Grid Impedance (Weak Inductive) ---
    P.R_line = 0.1; 
    P.X_line = 1.0; 
    
    % --- Fault Parameters ---
    P.t_fault_start = 1.0;
    P.fault_duration = 0.15;
    P.t_fault_clear = P.t_fault_start + P.fault_duration;
    P.V_fault = 0.1;
    
    % --- Simulation Setup ---
    t_end = 3.0;
    
    disp(['Running FAULT RIDE-THROUGH test for ' control_mode '...']);

    % ========================================================================
    % 2. SET UP AND RUN THE SIMULATION
    % ========================================================================
    options = odeset('RelTol', 1e-4, 'AbsTol', 1e-5);
    
    if strcmp(control_mode, 'VOC_AHO') || strcmp(control_mode, 'uVOC')
        % State for AHO/uVOC: x = [id, iq, v_alpha, v_beta]
        v_alpha0 = P.V_peak; v_beta0 = 0;
        Z = P.R_line + 1j*P.X_line;
        I_complex = ( (v_alpha0 + 1j*v_beta0) - P.V_grid_nominal) / Z;
        id0 = real(I_complex); iq0 = imag(I_complex);
        x0 = [id0; iq0; v_alpha0; v_beta0];
    else % Droop or VSM
        Z_line = sqrt(P.R_line^2 + P.X_line^2);
        delta0 = asin(P.P_ref * Z_line / P.V_grid_nominal);
        if isnan(delta0), delta0 = 0; end
        if strcmp(control_mode, 'Droop_Pure')
            x0 = [delta0; P.P_ref];
        else
            x0 = [delta0; 0; P.P_ref];
        end
    end

    [t1, x1] = ode23t(@(t,x) fault_dynamics(t, x, P, control_mode, 0), [0 P.t_fault_start], x0, options);
    [t2, x2] = ode23t(@(t,x) fault_dynamics(t, x, P, control_mode, 1), [P.t_fault_start P.t_fault_clear], x1(end,:)', options);
    [t3, x3] = ode23t(@(t,x) fault_dynamics(t, x, P, control_mode, 0), [P.t_fault_clear t_end], x2(end,:)', options);
    
    t = [t1; t2(2:end); t3(2:end)];
    x = [x1; x2(2:end,:); x3(2:end,:)];

    % ========================================================================
    % 3. POST-PROCESS AND PLOT RESULTS
    % ========================================================================
    disp('Simulation finished. Processing results...');
    
    if strcmp(control_mode, 'VOC_AHO') || strcmp(control_mode, 'uVOC')
        id = x(:,1); iq = x(:,2); v_alpha = x(:,3); v_beta = x(:,4);
        V_grid_vec = (t < P.t_fault_start | t > P.t_fault_clear) * P.V_grid_nominal + ...
                     (t >= P.t_fault_start & t <= P.t_fault_clear) * P.V_fault;
        P_elec = V_grid_vec .* id;
        w_inst = (v_alpha .* gradient(v_beta,t) - v_beta .* gradient(v_alpha,t)) ./ (v_alpha.^2 + v_beta.^2);
        freq_hz = w_inst / (2*pi);
        angle_deg = rad2deg(unwrap(atan2(v_beta, v_alpha)));
    else
        if strcmp(control_mode, 'Droop_Pure')
            delta_rad = x(:,1); p_filt = x(:,2);
            w_dev = -P.mp_droop * (p_filt - P.P_ref);
        else
            delta_rad = x(:,1); w_dev = x(:,2);
        end
        V_grid_vec = (t < P.t_fault_start | t > P.t_fault_clear) * P.V_grid_nominal + ...
                     (t >= P.t_fault_start & t <= P.t_fault_clear) * P.V_fault;
        Z_line = sqrt(P.R_line^2 + P.X_line^2);
        P_elec = (V_grid_vec / Z_line) .* sin(delta_rad);
        freq_hz = P.f_nom + w_dev / (2*pi);
        angle_deg = rad2deg(unwrap(delta_rad));
    end
    
    figure('Name', [control_mode ' Fault Ride-Through Response']);
    
    subplot(3,1,1); plot(t, freq_hz, 'b', 'LineWidth', 2);
    hold on; line([P.t_fault_start P.t_fault_start], ylim, 'Color', 'r', 'LineStyle', '--');
    line([P.t_fault_clear P.t_fault_clear], ylim, 'Color', 'k', 'LineStyle', '--');
    title('Inverter Frequency'); xlabel('Time (s)'); ylabel('Frequency (Hz)'); grid on;

    subplot(3,1,2); plot(t, P_elec, 'm', 'LineWidth', 2);
    hold on; line([P.t_fault_start P.t_fault_start], ylim, 'Color', 'r', 'LineStyle', '--');
    line([P.t_fault_clear P.t_fault_clear], ylim, 'Color', 'k', 'LineStyle', '--');
    title('Inverter Active Power Output'); xlabel('Time (s)'); ylabel('Power (pu)'); grid on;

    subplot(3,1,3); plot(t, angle_deg, 'g', 'LineWidth', 2);
    hold on; line([P.t_fault_start P.t_fault_start], ylim, 'Color', 'r', 'LineStyle', '--');
    line([P.t_fault_clear P.t_fault_clear], ylim, 'Color', 'k', 'LineStyle', '--');
    title('Inverter Internal Angle'); xlabel('Time (s)'); ylabel('Angle (degrees)'); grid on;
end

% ============================================================================
% THE DYNAMICS FUNCTION FOR FAULT SIMULATION
% ============================================================================
function dxdt = fault_dynamics(t, x, P, control_mode, fault_status)
    if fault_status == 1, V_grid_dynamic = P.V_fault;
    else, V_grid_dynamic = P.V_grid_nominal; end
    
    if strcmp(control_mode, 'VOC_AHO') || strcmp(control_mode, 'uVOC')
        id = x(1); iq = x(2); v_alpha = x(3); v_beta = x(4);
        
        theta_grid = P.w_nom * t;
        v_d_inv = v_alpha * cos(theta_grid) + v_beta * sin(theta_grid);
        v_q_inv = -v_alpha * sin(theta_grid) + v_beta * cos(theta_grid);
        
        did_dt = (1/P.X_line) * (-P.X_line*id + P.w_nom*P.X_line*iq + v_d_inv - V_grid_dynamic);
        diq_dt = (1/P.X_line) * (-P.X_line*iq - P.w_nom*P.X_line*id + v_q_inv - 0);
        
        % --- VOC Feedback Calculation ---
        if fault_status == 1 && strcmp(control_mode, 'uVOC')
            % --- uVOC: Voltage Synchronization during fault ---
            v_pcc_hat_d = V_grid_dynamic * cos(theta_grid);
            v_pcc_hat_q = V_grid_dynamic * sin(theta_grid);
            u = P.gamma * ([v_alpha; v_beta] - P.kv * [v_pcc_hat_d; v_pcc_hat_q]);
        else
            % --- Standard AHO: Current Synchronization ---
            i_alpha = id * cos(theta_grid) - iq * sin(theta_grid);
            i_beta = id * sin(theta_grid) + iq * cos(theta_grid);
            V_sq_inv = v_alpha^2 + v_beta^2;
            if V_sq_inv < 1e-3, V_sq_inv = 1e-3; end
            i_alpha_star = (2/3) * (v_alpha * P.P_ref + v_beta * P.Q_ref) / V_sq_inv;
            i_beta_star = (2/3) * (v_beta * P.P_ref - v_alpha * P.Q_ref) / V_sq_inv;
            R_phi = [cos(P.phi) -sin(P.phi); sin(P.phi) cos(P.phi)];
            i_err = [i_alpha - i_alpha_star; i_beta - i_beta_star];
            u = P.ki * R_phi * i_err;
        end
        
        V_sq_inv = v_alpha^2 + v_beta^2;
        V_nom_sq = P.kv^2;
        dv_alpha_dt = (P.xi/P.kv^2)*(2*V_nom_sq - V_sq_inv)*v_alpha - P.w_nom*v_beta - (P.kv/P.C_voc)*u(1);
        dv_beta_dt  = (P.xi/P.kv^2)*(2*V_nom_sq - V_sq_inv)*v_beta  + P.w_nom*v_alpha - (P.kv/P.C_voc)*u(2);
        
        dxdt = [did_dt; diq_dt; dv_alpha_dt; dv_beta_dt];
        
    else % Droop or VSM
        Z_line = sqrt(P.R_line^2 + P.X_line^2);
        if strcmp(control_mode, 'Droop_Pure')
            delta = x(1); p_filt = x(2);
            P_elec = (V_grid_dynamic / Z_line) * sin(delta);
            dp_filt_dt = (1/P.T_filter_droop) * (P_elec - p_filt);
            w_dev = -P.mp_droop * (p_filt - P.P_ref);
            ddelta_dt = w_dev;
            dxdt = [ddelta_dt; dp_filt_dt];
        else % VSM
            delta = x(1); w_dev = x(2); p_filt = x(3);
            P_elec = (V_grid_dynamic / Z_line) * sin(delta);
            dp_filt_dt = (1/P.T_filter_vsm) * (P_elec - p_filt);
            dw_dev_dt = (1/P.J_vsm) * (P.P_ref - p_filt - P.D_vsm * w_dev);
            ddelta_dt = w_dev;
            dxdt = [ddelta_dt; dw_dev_dt; dp_filt_dt];
        end
    end
end
