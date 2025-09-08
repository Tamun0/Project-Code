function run_harmonic_simulation_AHO_VOC()
    % Main function to run a robust simulation comparing the power quality
    % (harmonic performance) of different GFM controllers when
    % supplying a nonlinear load in an islanded microgrid.
    % This version includes a full Andronov-Hopf Oscillator (AHO) VOC model.

    clear; clc; close all;

    % ========================================================================
    % --- Simulation Case Selection ---
    % ========================================================================
    % Control Mode: 'Droop', 'VSM', or 'VOC_AHO'
    control_mode = 'VSM'; 
    
    % ========================================================================
    % 1. DEFINE SYSTEM PARAMETERS (Struct 'P')
    % ========================================================================
    P.f_nom = 50;
    P.w_nom = 2 * pi * P.f_nom;
    P.V_nom_rms = 230;
    P.V_peak = P.V_nom_rms * sqrt(2);
    
    % --- Power Setpoints & Loads ---
    P.P_ref = 10000; % 10 kW load
    P.Q_ref = 0;
    P.R_load = (P.V_nom_rms^2) / P.P_ref; % Linear load
    
    % --- Nonlinear Load Parameters ---
    P.I5_mag = 0.2 * (P.V_peak / P.R_load); % 20% 5th harmonic
    P.I7_mag = 0.14 * (P.V_peak / P.R_load);% 14% 7th harmonic
    
    % --- Controller Parameters ---
    % Inner loop gains for Droop & VSM
    P.kp_v = 0.5; P.ki_v = 150; P.kp_i = 10; P.ki_i = 200;
    % VSM-based model
    P.J_vsm = 0.2; P.D_vsm = 5.0; 
    P.T_filter = 0.05; 
    % Pure Droop model
    P.mp_droop = (0.05 * P.w_nom) / P.P_ref;
    % AHO-VOC model
    P.xi = 25; P.kv = P.V_peak; P.ki = 3 * P.V_peak / (P.P_ref + 1);
    P.C_voc = 0.1; P.phi = pi/2;

    % --- Circuit Parameters ---
    P.Lf = 1e-3; P.Rf = 0.1; P.Cf = 100e-6;
    
    % --- Simulation Setup ---
    P.t_harmonic_start = 0.2;
    t_span = [0 0.4];
    
    disp(['Running ROBUST HARMONIC PERFORMANCE test for ' control_mode '...']);

    % ========================================================================
    % 2. SET UP AND RUN THE SIMULATION
    % ========================================================================
    if strcmp(control_mode, 'VOC_AHO')
        % State for AHO-VOC: x = [iL_a, iL_b, iL_c, vC_a, vC_b, vC_c, v_alpha, v_beta]
        x0 = zeros(8, 1);
        x0(7) = P.V_peak; % Start with initial voltage
    elseif strcmp(control_mode, 'Droop')
        % State for Pure Droop: x = [iL_a..c, vC_a..c, delta, p_filt, ...
        %      err_v_d_int, err_v_q_int, err_i_d_int, err_i_q_int]
        x0 = zeros(12, 1);
    else
        % State for VSM-based models: x = [iL_a..c, vC_a..c, delta, w_dev, p_filt, ...
        %      err_v_d_int, err_v_q_int, err_i_d_int, err_i_q_int]
        x0 = zeros(13, 1);
    end
    
    options = odeset('RelTol', 1e-4, 'AbsTol', 1e-5);
    [t, x] = ode15s(@(t,x) harmonic_dynamics_abc(t, x, P, control_mode), t_span, x0, options);

    % ========================================================================
    % 3. POST-PROCESS AND PLOT RESULTS
    % ========================================================================
    disp('Simulation finished. Processing results...');
    
    vC_a = x(:,4);
    
    fs = 1 / (t(2)-t(1)); 
    idx_start = find(t >= 0.3, 1);
    thd_val_percent = 0;
    if ~isempty(idx_start)
        thd_db = thd(vC_a(idx_start:end), fs, 7);
        thd_val_percent = 100 * 10^(thd_db / 20);
    end
    
    figure('Name', [control_mode ' Harmonic Performance']);
    plot(t, vC_a, 'b', 'LineWidth', 1.5);
    hold on;
    line([P.t_harmonic_start P.t_harmonic_start], ylim, 'Color', 'r', 'LineStyle', '--');
    title_str = sprintf('%s Response to Nonlinear Load\nCalculated THD = %.2f%%', strrep(control_mode, '_', ' '), thd_val_percent);
    title(title_str);
    xlabel('Time (s)'); ylabel('Phase A Voltage (V)');
    legend('V_a', 'Nonlinear Load Connected');
    grid on;
end

% ============================================================================
% THE DYNAMICS FUNCTION (in abc frame)
% ============================================================================
function dxdt = harmonic_dynamics_abc(t, x, P, control_mode)
    % --- Unpack state vector based on control mode ---
    iL_abc = x(1:3); vC_abc = x(4:6);
    
    if strcmp(control_mode, 'VOC_AHO')
        v_alpha = x(7); v_beta = x(8);
        w_inst = P.w_nom; % AHO has its own internal frequency
        theta = atan2(v_beta, v_alpha);
    elseif strcmp(control_mode, 'Droop')
        delta = x(7); p_filt = x(8);
        err_v_d_int = x(9); err_v_q_int = x(10);
        err_i_d_int = x(11); err_i_q_int = x(12);
        w_dev = -P.mp_droop * (p_filt - P.P_ref);
        w_inst = P.w_nom + w_dev;
        theta = mod(w_inst * t + delta, 2*pi);
    else % VSM-based Droop or VSM
        delta = x(7); w_dev = x(8); p_filt = x(9);
        err_v_d_int = x(10); err_v_q_int = x(11);
        err_i_d_int = x(12); err_i_q_int = x(13);
        w_inst = P.w_nom + w_dev;
        theta = mod(w_inst * t + delta, 2*pi);
    end

    % --- Controller Logic ---
    if strcmp(control_mode, 'VOC_AHO')
        % --- AHO-VOC Dynamics ---
        % Convert measured current to alpha-beta frame
        T_clarke = (2/3) * [1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2];
        iL_ab = T_clarke * iL_abc;
        i_alpha = iL_ab(1); i_beta = iL_ab(2);
        
        % Current setpoints from power setpoints
        V_sq_inv = v_alpha^2 + v_beta^2;
        if V_sq_inv < 1e-3, V_sq_inv = 1e-3; end
        i_alpha_star = (2/3) * (v_alpha * P.P_ref + v_beta * P.Q_ref) / V_sq_inv;
        i_beta_star = (2/3) * (v_beta * P.P_ref - v_alpha * P.Q_ref) / V_sq_inv;
        
        % Current error with rotation
        R_phi = [cos(P.phi) -sin(P.phi); sin(P.phi) cos(P.phi)];
        i_err = [i_alpha - i_alpha_star; i_beta - i_beta_star];
        u = P.ki * R_phi * i_err;
        
        % Andronov-Hopf Oscillator Equations
        V_nom_sq = P.kv^2;
        dv_alpha_dt = (P.xi/P.kv^2)*(2*V_nom_sq - V_sq_inv)*v_alpha - P.w_nom*v_beta - (P.kv/P.C_voc)*u(1);
        dv_beta_dt  = (P.xi/P.kv^2)*(2*V_nom_sq - V_sq_inv)*v_beta  + P.w_nom*v_alpha - (P.kv/P.C_voc)*u(2);
        
        % The oscillator output directly becomes the inverter voltage reference
        v_inv_ref_ab = [v_alpha; v_beta];
        T_inv_clarke = [1 0; -1/2 sqrt(3)/2; -1/2 -sqrt(3)/2];
        v_inv_ref_abc = T_inv_clarke * v_inv_ref_ab;
        
        % Set derivatives for unused states to zero
        ddelta_dt = 0; dw_dev_dt = 0; dp_filt_dt = 0;
        derr_v_d_int_dt=0; derr_v_q_int_dt=0; derr_i_d_int_dt=0; derr_i_q_int_dt=0;

    else % Droop or VSM dynamics
        T_park = (2/3) * [cos(theta) cos(theta - 2*pi/3) cos(theta + 2*pi/3);
                          -sin(theta) -sin(theta - 2*pi/3) -sin(theta + 2*pi/3)];
        vC_dq = T_park * vC_abc;
        iL_dq = T_park * iL_abc;
        P_elec = vC_dq(1)*iL_dq(1) + vC_dq(2)*iL_dq(2);
        
        % --- Outer GFM Loop ---
        if strcmp(control_mode, 'Droop')
            dp_filt_dt = (1/P.T_filter) * (P_elec - p_filt);
            dw_dev_dt = []; % Not a state
        else % VSM
            dp_filt_dt = (1/P.T_filter) * (P_elec - p_filt);
            dw_dev_dt = (1/P.J_vsm) * (P.P_ref - p_filt - P.D_vsm * w_dev);
        end
        
        % --- Inner Cascaded Voltage and Current Control ---
        v_d_ref = P.V_peak; v_q_ref = 0;
        err_v_d = v_d_ref - vC_dq(1); err_v_q = v_q_ref - vC_dq(2);
        derr_v_d_int_dt = err_v_d; derr_v_q_int_dt = err_v_q;
        id_ref = P.kp_v*err_v_d + P.ki_v*err_v_d_int;
        iq_ref = P.kp_v*err_v_q + P.ki_v*err_v_q_int;
        err_i_d = id_ref - iL_dq(1); err_i_q = iq_ref - iL_dq(2);
        derr_i_d_int_dt = err_i_d; derr_i_q_int_dt = err_i_q;
        v_inv_ref_d = P.kp_i*err_i_d + P.ki_i*err_i_d_int;
        v_inv_ref_q = P.kp_i*err_i_q + P.ki_i*err_i_q_int;
        v_inv_ref_dq = [v_inv_ref_d; v_inv_ref_q];
        
        ddelta_dt = w_dev;
        T_inv_park = [cos(theta) -sin(theta); cos(theta-2*pi/3) -sin(theta-2*pi/3); cos(theta+2*pi/3) -sin(theta+2*pi/3)];
        v_inv_ref_abc = T_inv_park * v_inv_ref_dq;
    end

    % --- Load Current Calculation ---
    i_linear_load_abc = vC_abc / P.R_load;
    i_harmonic_abc = [0;0;0];
    if t >= P.t_harmonic_start
        i_harmonic_abc = [ P.I5_mag * sin(5 * w_inst * t); P.I5_mag * sin(5 * (w_inst * t - 2*pi/3)); P.I5_mag * sin(5 * (w_inst * t + 2*pi/3)) ] + ...
                         [ P.I7_mag * sin(7 * w_inst * t); P.I7_mag * sin(7 * (w_inst * t - 2*pi/3)); P.I7_mag * sin(7 * (w_inst * t + 2*pi/3)) ];
    end
    i_total_load_abc = i_linear_load_abc + i_harmonic_abc;
    
    % --- Circuit Differential Equations (in abc frame) ---
    diL_abc_dt = (1/P.Lf) * (v_inv_ref_abc - vC_abc - iL_abc * P.Rf);
    dvC_abc_dt = (1/P.Cf) * (iL_abc - i_total_load_abc);
    
    % --- Assemble state derivative vector ---
    if strcmp(control_mode, 'VOC_AHO')
        dxdt = [diL_abc_dt; dvC_abc_dt; dv_alpha_dt; dv_beta_dt];
    elseif strcmp(control_mode, 'Droop')
        dxdt = [diL_abc_dt; dvC_abc_dt; ddelta_dt; dp_filt_dt; ...
                derr_v_d_int_dt; derr_v_q_int_dt; derr_i_d_int_dt; derr_i_q_int_dt];
    else
        dxdt = [diL_abc_dt; dvC_abc_dt; ddelta_dt; dw_dev_dt; dp_filt_dt; ...
                derr_v_d_int_dt; derr_v_q_int_dt; derr_i_d_int_dt; derr_i_q_int_dt];
    end
end
