function run_voc_vs_droop_final1()
    % Main function to compare the dynamic response of an AHO-based VOC
    % and a classic Droop controller, based on the models from the provided paper.

    clear; clc; close all;

    % --- Simulation Case Selection ---
    control_mode = 'VOC'; % CHANGE THIS to 'VOC' or 'Droop'

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
    
    % --- VOC Behavioral Parameter ---
    P.T_voc = 0.01; % Fast time constant for VOC response
    
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
        % x_voc = [P_out]
        x0 = P.P_initial;
        [t, x] = ode23t(@(t,x) voc_dynamics_behavioral(t, x, P), t_span, x0);
        P_out = x(:,1);
        V_rms = ones(size(t)) * P.V_nom_rms; % Assume perfect voltage regulation
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
    sgtitle(['Dynamic Response of ' control_mode ' Controller'])
    subplot(2,1,1);
    plot(t, V_rms, 'b', 'LineWidth', 2);
    hold on;
    line([P.t_step P.t_step], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    title('Inverter Terminal Voltage Magnitude');
    xlabel('Time (s)'); ylabel('Voltage (V_{RMS})');
    legend('Voltage Magnitude','Disturbance');
    grid on;

    subplot(2,1,2);
    plot(t, P_out, 'm', 'LineWidth', 2);
    hold on;
    line([P.t_step P.t_step], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    title('Inverter Active Power Output');
    xlabel('Time (s)'); ylabel('Power (W)');
    legend('Power Output','Disturbance');
    grid on;
end

% ============================================================================
% DYNAMICS FUNCTIONS (Based on the paper)
% ============================================================================

function dxdt = voc_dynamics_behavioral(t, x, P)
    % Implements a stable, behavioral model of the AHO-VOC response.
    % It models the fast, non-overshooting behavior described in the paper.
    P_out = x(1);
    
    % Determine power reference based on time
    if t < P.t_step
        P_star = P.P_initial;
    else
        P_star = P.P_final;
    end
    
    % Model the power output as a simple first-order system
    % This creates a fast, exponential rise to the target without overshoot.
    dP_out_dt = (1 / P.T_voc) * (P_star - P_out);
    
    dxdt = dP_out_dt;
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
