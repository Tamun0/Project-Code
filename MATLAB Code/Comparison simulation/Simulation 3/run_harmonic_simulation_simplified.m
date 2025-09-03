function run_harmonic_simulation_simplified()
    % Main function to run a simulation comparing the power quality
    % of simplified Droop, VSM, and VOC models without inner control loops.

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
    P.J_droop = 0.2; P.D_droop = 5.0; P.T_filter_droop = 0.05; 
    P.J_vsm = 5.0;  P.D_vsm = 1.5; P.T_filter_vsm = 0.02;
    P.T_voc = 0.02;

    % --- Circuit Parameters ---
    P.Lf = 1e-3;
    P.Rf = 0.1;
    P.Cf = 100e-6;
    
    % --- Simulation Setup ---
    P.t_harmonic_start = 0.2;
    t_span = [0 0.4];
    
    disp(['Running SIMPLIFIED HARMONIC PERFORMANCE test for ' control_mode '...']);

    % ========================================================================
    % 2. SET UP AND RUN THE SIMULATION
    % ========================================================================
    % x = [iL_a, iL_b, iL_c, vC_a, vC_b, vC_c, delta, w_dev, p_filt]
    x0 = zeros(9, 1);
    x0(4:6) = [P.V_peak; P.V_peak*cos(-2*pi/3); P.V_peak*cos(2*pi/3)]; % Start with stable voltage
    
    options = odeset('RelTol', 1e-4, 'AbsTol', 1e-5);
    [t, x] = ode15s(@(t,x) harmonic_dynamics_abc_simplified(t, x, P, control_mode), t_span, x0, options);

    % ========================================================================
    % 3. POST-PROCESS AND PLOT RESULTS
    % ========================================================================
    disp('Simulation finished. Processing results...');
    
    vC_a = x(:,4);
    
    fs = 1 / (t(2)-t(1)); 
    idx_start = find(t >= 0.3, 1);
    thd_val_percent = 0;
    if ~isempty(idx_start) && length(vC_a(idx_start:end)) > 1
        thd_db = thd(vC_a(idx_start:end), fs, 7);
        thd_val_percent = 100 * 10^(thd_db / 20);
    end
    
    figure('Name', [control_mode ' Harmonic Performance (Simplified)']);
    
    plot(t, vC_a, 'b', 'LineWidth', 1.5);
    hold on;
    line([P.t_harmonic_start P.t_harmonic_start], ylim, 'Color', 'r', 'LineStyle', '--');
    title_str = sprintf('%s Response to Nonlinear Load\nCalculated THD = %.2f%%', control_mode, thd_val_percent);
    title(title_str);
    xlabel('Time (s)'); ylabel('Phase A Voltage (V)');
    legend('V_a', 'Nonlinear Load Connected');
    grid on;
end

% ============================================================================
% THE DYNAMICS FUNCTION (Simplified Models)
% ============================================================================
function dxdt = harmonic_dynamics_abc_simplified(t, x, P, control_mode)
    % Unpack state vector
    iL_abc = x(1:3);
    vC_abc = x(4:6);
    delta = x(7);
    w_dev = x(8);
    p_filt = x(9);

    % --- Controller Logic (operates in dq frame) ---
    w_inst = P.w_nom + w_dev;
    theta = mod(w_inst * t + delta, 2*pi);
    
    T_park = (2/3) * [cos(theta) cos(theta - 2*pi/3) cos(theta + 2*pi/3);
                      -sin(theta) -sin(theta - 2*pi/3) -sin(theta + 2*pi/3)];
    vC_dq = T_park * vC_abc;
    iL_dq = T_park * iL_abc;
    
    P_elec = vC_dq(1)*iL_dq(1) + vC_dq(2)*iL_dq(2);
    
    % --- Outer GFM Loop ---
    P_ref = P.P_ref;
    switch control_mode
        case {'Droop', 'VSM'}
            if strcmp(control_mode, 'Droop')
                dp_filt_dt = (1/P.T_filter_droop) * (P_elec - p_filt);
                dw_dev_dt = (1/P.J_droop) * (P_ref - p_filt - P.D_droop * w_dev);
            else % VSM
                dp_filt_dt = (1/P.T_filter_vsm) * (P_elec - p_filt);
                dw_dev_dt = (1/P.J_vsm) * (P_ref - p_filt - P.D_vsm * w_dev);
            end
            % FIX: No inner loops. GFM loop directly sets voltage reference.
            v_inv_ref_dq = [P.V_peak; 0];

        case 'VOC'
            dp_filt_dt = (1/P.T_voc) * (P_ref - p_filt);
            target_delta = asin(p_filt * P.Lf * w_inst / P.V_peak^2);
            if isnan(target_delta), target_delta = delta; end
            dw_dev_dt = 20 * (target_delta - delta) - 2 * w_dev; 
            v_inv_ref_dq = [P.V_peak; 0];
    end
    ddelta_dt = w_dev;

    T_inv_park = [cos(theta) -sin(theta); cos(theta-2*pi/3) -sin(theta-2*pi/3); cos(theta+2*pi/3) -sin(theta+2*pi/3)];
    v_inv_ref_abc = T_inv_park * v_inv_ref_dq;

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
    
    dxdt = [diL_abc_dt; dvC_abc_dt; ddelta_dt; dw_dev_dt; dp_filt_dt];
end
