function run_voc_vs_droop_final()
    % Main function to compare the dynamic response of an AHO-based VOC
    % and a classic Droop controller, based on the models from the provided paper.

    clear; clc; close all;

    % --- Simulation Case Selection ---
    control_mode = 'Droop'; % CHANGE THIS to 'VOC' or 'Droop'

    % ========================================================================
    % 1. DEFINE SYSTEM PARAMETERS (Struct 'P')
    % Based on "Comparison of Droop Control..." paper
    % ========================================================================
    P.f_nom = 60; 
    P.w_nom = 2 * pi * P.f_nom;
    P.V_nom_rms = 120; 
    P.V_peak = P.V_nom_rms * sqrt(2);
    
    % --- Line/Filter Parameters ---
    P.Lf = 1.5e-3; 
    P.Rf = 0.8;    
    
    % --- Power Setpoints ---
    P.P_initial = 0;    
    P.P_final = 500;    
    P.Q_ref = 0;        
    
    % --- AHO-VOC Controller Parameters (from paper Table II) ---
    P.kv = 120;
    P.ki = 0.24;
    P.C_voc = 0.2679; % This is a control parameter, not a physical capacitance
    P.xi = 15;
    
    % --- Droop Controller Parameters (from paper Table II) ---
    P.mp = 2.6e-3; 
    P.mq = 5.0e-3;
    P.wc = 2 * pi * 30; % 30 Hz cutoff frequency for power filter

    % --- Simulation Setup ---
    P.t_step = 0.2; 
    t_span = [0 0.4];
    
    disp(['Running comparison simulation for ' control_mode ' control...']);

    % ========================================================================
    % 2. SET UP AND RUN THE SIMULATION
    % ========================================================================
    if strcmp(control_mode, 'VOC')
        % x_voc = [v_alpha, v_beta, i_d, i_q]
        x0 = [P.V_peak; 0; 0; 0];
        [t, x] = ode23t(@(t,x) voc_dynamics(t, x, P), t_span, x0);
        V_rms = sqrt(x(:,1).^2 + x(:,2).^2) / sqrt(2);
        P_out = (3/2) * (x(:,1).*x(:,3) + x(:,2).*x(:,4));
    else % Droop
        % x_droop = [delta, V_peak, P_filt, Q_filt, id, iq]
        x0 = [0; P.V_peak; P.P_initial; P.Q_ref; 0; 0];
        [t, x] = ode23t(@(t,x) droop_dynamics(t, x, P), t_span, x0);
        V_rms = x(:,2) / sqrt(2);
        P_out = x(:,3); % Plotting the filtered power
    end

    % ========================================================================
    % 3. PLOT RESULTS
    % ========================================================================
    disp('Simulation finished. Plotting results...');
    
    figure('Name', ['Dynamic Response of ' control_mode ' Controller']);
    
    subplot(2,1,1);
    plot(t, V_rms, 'b', 'LineWidth', 2);
    hold on;
    line([P.t_step P.t_step], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    title('Inverter Terminal Voltage Magnitude');
    xlabel('Time (s)'); ylabel('Voltage (V_{RMS})');
    grid on;

    subplot(2,1,2);
    plot(t, P_out, 'm', 'LineWidth', 2);
    hold on;
    line([P.t_step P.t_step], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    title('Inverter Active Power Output');
    xlabel('Time (s)'); ylabel('Power (W)');
    grid on;
end

% ============================================================================
% DYNAMICS FUNCTIONS (Based on the paper)
% ============================================================================

function dxdt = voc_dynamics(t, x, P)
    % Implements the AHO-VOC dynamics from Eq. (7) & (2) of the paper
    v_alpha = x(1); v_beta = x(2);
    id = x(3); iq = x(4); 
    
    if t < P.t_step, P_star = P.P_initial; else, P_star = P.P_final; end
    
    V_sq = v_alpha^2 + v_beta^2;
    V_abs = sqrt(V_sq);
    
    i_alpha_star = (2/3) * (v_alpha * P_star - v_beta * P.Q_ref) / (V_sq + 1e-6);
    i_beta_star  = (2/3) * (v_beta * P_star + v_alpha * P.Q_ref) / (V_sq + 1e-6);
    
    theta = atan2(v_beta, v_alpha);
    i_alpha = id * cos(theta) - iq * sin(theta);
    i_beta  = id * sin(theta) + iq * cos(theta);

    % VOC Dynamics Eq. (7)
    dv_alpha_dt = (P.xi/P.kv^2)*(P.V_peak^2-V_sq)*v_alpha - P.w_nom*v_beta - (P.kv/P.C_voc)*P.ki*(i_alpha-i_alpha_star);
    dv_beta_dt  = P.w_nom*v_alpha + (P.xi/P.kv^2)*(P.V_peak^2-V_sq)*v_beta - (P.kv/P.C_voc)*P.ki*(i_beta-i_beta_star);

    % Line/Filter Dynamics Eq. (2)
    e_d = P.V_peak * cos(theta); 
    e_q = -P.V_peak * sin(theta);
    w_i = P.w_nom - (P.kv*P.ki / (3*P.C_voc*V_sq + 1e-6)) * ((3/2)*(v_alpha*i_alpha + v_beta*i_beta) - P_star);
    
    did_dt = -P.Rf/P.Lf * id + w_i*iq + (1/P.Lf)*(V_abs - e_d);
    diq_dt = -w_i*id - P.Rf/P.Lf * iq - (1/P.Lf)*e_q;
    
    dxdt = [dv_alpha_dt; dv_beta_dt; did_dt; diq_dt];
end

function dxdt = droop_dynamics(t, x, P)
    % Implements the Droop control dynamics from Eq. (11) & (12) & (2)
    delta = x(1); V_peak = x(2);
    P_filt = x(3); Q_filt = x(4);
    id = x(5); iq = x(6);
    
    if t < P.t_step, P_star = P.P_initial; else, P_star = P.P_final; end
    
    % Instantaneous Power
    P_inst = (3/2) * V_peak * id;
    Q_inst = -(3/2) * V_peak * iq;
    
    % Power Filter Dynamics Eq. (11)
    dP_filt_dt = P.wc * (P_inst - P_filt);
    dQ_filt_dt = P.wc * (Q_inst - Q_filt);
    
    % Droop Dynamics Eq. (12) - Using a more stable formulation
    V_ref = P.V_peak - P.mq * (Q_filt - P.Q_ref);
    dV_dt = 10 * P.wc * (V_ref - V_peak); % V tracks its reference
    
    w_i = P.w_nom - P.mp * (P_filt - P_star);
    ddelta_dt = w_i - P.w_nom; 
    
    % Line/Filter Dynamics Eq. (2)
    e_d = P.V_peak * cos(delta);
    e_q = -P.V_peak * sin(delta);
    
    did_dt = -P.Rf/P.Lf * id + w_i*iq + (1/P.Lf)*(V_peak - e_d);
    diq_dt = -w_i*id - P.Rf/P.Lf * iq - (1/P.Lf)*e_q;
    
    dxdt = [ddelta_dt; dV_dt; dP_filt_dt; dQ_filt_dt; did_dt; diq_dt];
end
