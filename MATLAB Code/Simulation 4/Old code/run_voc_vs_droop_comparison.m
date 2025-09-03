function run_voc_vs_droop_comparison()
    % Main function to compare the dynamic response of an AHO-based VOC
    % and a classic Droop controller, based on the models from the provided paper.

    clear; clc; close all;

    % --- Simulation Case Selection ---
    control_mode = 'Droop'; % CHANGE THIS to 'VOC' or 'Droop'

    % ========================================================================
    % 1. DEFINE SYSTEM PARAMETERS (Struct 'P')
    % Based on "Comparison of Droop Control..." paper
    % ========================================================================
    P.f_nom = 60; % Using 60 Hz as is common in US-based papers
    P.w_nom = 2 * pi * P.f_nom;
    P.V_nom = 120; % Nominal phase voltage (RMS)
    P.V_peak = P.V_nom * sqrt(2);
    
    % --- Line/Filter Parameters ---
    P.Lf = 1.5e-3; % Filter inductance
    P.Rf = 0.8;    % Filter resistance
    
    % --- Power Setpoints ---
    P.P_initial = 0;    % Start at 0 W
    P.P_final = 500;    % Step to 500 W
    P.Q_ref = 0;        % Reactive power reference
    
    % --- AHO-VOC Controller Parameters (from paper) ---
    P.kv = 120;
    P.ki = 0.24;
    P.C = 0.2679e-6; % Scaled capacitance
    P.xi = 15;
    
    % --- Droop Controller Parameters (from paper) ---
    P.mp = 2.6e-3 * P.w_nom; % Adjusted for rad/s
    P.mq = 5.0e-3;
    P.wc = 2 * pi * 30; % 30 Hz cutoff frequency for power filter

    % --- Simulation Setup ---
    P.t_step = 0.2; % Time of the power reference step
    t_span = [0 0.4];
    
    disp(['Running comparison simulation for ' control_mode ' control...']);

    % ========================================================================
    % 2. SET UP AND RUN THE SIMULATION
    % ========================================================================
    if strcmp(control_mode, 'VOC')
        % x_voc = [v_alpha, v_beta, i_d, i_q]
        x0 = [P.V_peak; 0; 0; 0];
        [t, x] = ode23t(@(t,x) voc_dynamics(t, x, P), t_span, x0);
        V = sqrt(x(:,1).^2 + x(:,2).^2) / sqrt(2); % RMS Voltage
        P_out = (3/2) * (x(:,1).*x(:,3) + x(:,2).*x(:,4));
    else % Droop
        % x_droop = [delta, V, P_filt, Q_filt, id, iq]
        delta0 = 0;
        V0 = P.V_peak / sqrt(2);
        x0 = [delta0; V0; P.P_initial; P.Q_ref; 0; 0];
        [t, x] = ode23t(@(t,x) droop_dynamics(t, x, P), t_span, x0);
        V = x(:,2);
        P_out = x(:,3);
    end

    % ========================================================================
    % 3. PLOT RESULTS
    % ========================================================================
    disp('Simulation finished. Plotting results...');
    
    figure('Name', ['Dynamic Response of ' control_mode ' Controller']);
    
    % Plot Terminal Voltage
    subplot(2,1,1);
    plot(t, V, 'b', 'LineWidth', 2);
    hold on;
    line([P.t_step P.t_step], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    title('Inverter Terminal Voltage Magnitude');
    xlabel('Time (s)'); ylabel('Voltage (V_{RMS})');
    grid on;

    % Plot Power Output
    subplot(2,1,2);
    plot(t, P_out, 'm', 'LineWidth', 2);
    hold on;
    line([P.t_step P.t_step], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    title('Inverter Active Power Output');
    xlabel('Time (s)'); ylabel('Power (W)');
    grid on;
end

% ============================================================================
% DYNAMICS FUNCTIONS
% ============================================================================

function dxdt = voc_dynamics(t, x, P)
    % Implements the AHO-VOC dynamics from Eq. (7) of the paper
    v_alpha = x(1); v_beta = x(2);
    id = x(3); iq = x(4); % Currents are in the inverter's dq frame
    
    % Determine power reference
    if t < P.t_step
        P_star = P.P_initial;
    else
        P_star = P.P_final;
    end
    
    V_sq = v_alpha^2 + v_beta^2;
    V_abs = sqrt(V_sq);
    
    % Calculate current references from power setpoints (simplified)
    i_alpha_star = (2/3) * (v_alpha * P_star - v_beta * P.Q_ref) / V_sq;
    i_beta_star  = (2/3) * (v_beta * P_star + v_alpha * P.Q_ref) / V_sq;
    
    % Convert measured currents to alpha-beta frame
    theta = atan2(v_beta, v_alpha);
    i_alpha = id * cos(theta) - iq * sin(theta);
    i_beta  = id * sin(theta) + iq * cos(theta);

    % VOC Dynamics from Eq. (7)
    dv_alpha_dt = (P.xi/P.kv^2)*(P.V_peak^2-V_sq)*v_alpha - P.w_nom*v_beta - (P.kv/P.C)*P.ki*(i_alpha-i_alpha_star);
    dv_beta_dt  = P.w_nom*v_alpha + (P.xi/P.kv^2)*(P.V_peak^2-V_sq)*v_beta - (P.kv/P.C)*P.ki*(i_beta-i_beta_star);

    % Line/Filter Dynamics from Eq. (2)
    e_d = P.V_peak * cos(theta); % Assume grid angle is 0
    e_q = -P.V_peak * sin(theta);
    w_i = P.w_nom - (P.kv*P.ki / (3*P.C*V_sq)) * ((3/2)*(v_alpha*i_alpha + v_beta*i_beta) - P_star);
    
    did_dt = -P.Rf/P.Lf * id + w_i*iq + (1/P.Lf)*(V_abs - e_d);
    diq_dt = -w_i*id - P.Rf/P.Lf * iq + (1/P.Lf)*(0 - e_q);
    
    dxdt = [dv_alpha_dt; dv_beta_dt; did_dt; diq_dt];
end

function dxdt = droop_dynamics(t, x, P)
    % Implements the Droop control dynamics from Eq. (11) & (12)
    delta = x(1); V_rms = x(2);
    P_filt = x(3); Q_filt = x(4);
    id = x(5); iq = x(6);
    
    % Determine power reference
    if t < P.t_step
        P_star = P.P_initial;
    else
        P_star = P.P_final;
    end
    
    % Calculate instantaneous power
    V_peak = V_rms * sqrt(2);
    P_inst = (3/2) * V_peak * id;
    Q_inst = -(3/2) * V_peak * iq;
    
    % Power Filter Dynamics Eq. (11)
    dP_filt_dt = -P.wc * (P_filt - P_inst);
    dQ_filt_dt = -P.wc * (Q_filt - Q_inst);
    
    % Droop Dynamics Eq. (12)
    dV_dt = -P.mq * dQ_filt_dt;
    w_i = P.w_nom - P.mp * (P_filt - P_star);
    ddelta_dt = w_i - P.w_nom; % Assuming grid frequency is w_nom
    
    % Line/Filter Dynamics Eq. (2)
    e_d = P.V_peak * cos(delta);
    e_q = -P.V_peak * sin(delta);
    
    did_dt = -P.Rf/P.Lf * id + w_i*iq + (1/P.Lf)*(V_peak - e_d);
    diq_dt = -w_i*id - P.Rf/P.Lf * iq + (1/P.Lf)*(0 - e_q);
    
    dxdt = [ddelta_dt; dV_dt; dP_filt_dt; dQ_filt_dt; did_dt; diq_dt];
end
