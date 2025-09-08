function run_inverter_comparison_v7_1()
    % Main function to run a comparative simulation of Droop, VSM, and VOC
    % controllers in either a strong or weak grid environment, focusing on
    % the three-phase voltage and current waveforms.

    clear; clc; close all;

    % --- Simulation Case Selection ---
    grid_type = 'weak'; % Can be 'strong' or 'weak'

    % ========================================================================
    % 1. DEFINE SYSTEM AND CONTROLLER PARAMETERS
    % ========================================================================
    P.f_nom = 50;
    P.w_nom = 2 * pi * P.f_nom;
    P.V_nom_rms = 230;
    P.V_peak = P.V_nom_rms * sqrt(2);
    
    P.P_final = 10000; % Final active power setpoint (W)
    P.Q_final = 0;     % Final reactive power setpoint (VAR)
    
    % --- Controller Parameters ---
    % Droop
    P.mp_droop = (0.05 * P.w_nom) / P.P_final;
    P.mq_droop = (0.05 * P.V_peak) / (abs(P.Q_final) + 1e-6);
    P.T_filter_droop = 0.05;
    % VSM
    P.J_vsm = 0.5;
    P.D_vsm = 50.0;
    P.T_filter_vsm = 0.02;
    % VOC (Behavioral Model)
    P.T_voc = 0.02;
    
    % --- FIX: Add Virtual Inductance for Numerical Damping ---
    P.L_virtual = 1e-3; % 1 mH of virtual inductance to soften the model

    % --- Grid Impedance ---
    if strcmp(grid_type, 'strong')
        P.R_grid = 0.1; P.L_grid = 0.1e-3;
        disp('Running STRONG GRID comparison...');
    else % weak
        P.R_grid = 2.0; P.L_grid = 5.0e-3;
        disp('Running WEAK GRID comparison...');
    end
    
    P.L_filt = 2e-3; P.R_filt = 0.1;
    P.t_ramp_end = 0.3;
    t_span = [0 0.5];
    
    % ========================================================================
    % 2. RUN SIMULATION FOR EACH CONTROLLER
    % ========================================================================
    controllers = {'Droop', 'VSM', 'VOC'};
    results = cell(1, length(controllers));

    for i = 1:length(controllers)
        control_mode = controllers{i};
        fprintf('Simulating %s controller...\n', control_mode);
        
        if strcmp(control_mode, 'VOC')
            x0 = [0; 0; 0];
        else
            x0 = zeros(6, 1);
        end
        
        options = odeset('RelTol', 1e-4, 'AbsTol', 1e-5);
        [t, x] = ode23t(@(t,x) system_dynamics_dq_comparative(t, x, P, control_mode), t_span, x0, options);
        
        results{i}.t = t;
        results{i}.x = x;
        results{i}.name = control_mode;
    end

    % ========================================================================
    % 3. POST-PROCESS AND PLOT SEPARATE RESULTS FOR EACH CONTROLLER
    % ========================================================================
    disp('Simulations finished. Processing results...');
    
    for i = 1:length(results)
        % Create a new figure for each controller
        figure('Name', [results{i}.name ' Waveforms for ' grid_type ' Grid'], 'Position', [100 + (i-1)*450, 100, 1200, 500]);
        
        % Post-process to get waveform data
        [v_abc, i_abc] = post_process_waveforms(results{i}.t, results{i}.x, P, results{i}.name);
        
        % --- Subfigure 1: Terminal Voltage (ABC Waveforms) Plot ---
        subplot(1,2,1);
        plot(results{i}.t, v_abc, 'LineWidth', 1.5);
        title({[results{i}.name ' Controller']; ['Three-Phase Terminal Voltage (' grid_type ' Grid)']});
        xlabel('Time (s)'); ylabel('Voltage (V)');
        legend('Va', 'Vb', 'Vc', 'Location', 'northeast');
        grid on; ylim([-400 400]);

        % --- Subfigure 2: Output Current (ABC Waveforms) Plot ---
        subplot(1,2,2);
        plot(results{i}.t, i_abc, 'LineWidth', 1.5);
        title({[results{i}.name ' Controller']; 'Three-Phase Output Current'});
        xlabel('Time (s)'); ylabel('Current (A)');
        legend('Ia', 'Ib', 'Ic', 'Location', 'northeast');
        grid on;
    end
end

% ============================================================================
% POST-PROCESSING FUNCTION FOR WAVEFORMS
% ============================================================================
function [v_abc_terminal, i_abc] = post_process_waveforms(t, x, P, mode)
    v_abc_terminal = zeros(length(t), 3);
    i_abc = zeros(length(t), 3);

    if strcmp(mode, 'VOC')
        % For VOC behavioral model, we need to approximate the waveforms
        p_filt = x(:,3);
        % P = 3/2*V_peak*I_peak for balanced system in this frame
        I_peak_approx = p_filt / (1.5 * P.V_peak); 
        % Approximate sag using total impedance magnitude
        Z_total = sqrt((P.R_filt + P.R_grid)^2 + (P.w_nom * (P.L_filt + P.L_grid + P.L_virtual))^2);
        V_peak_approx = P.V_peak - I_peak_approx .* Z_total; 

        for k = 1:length(t)
            grid_angle = P.w_nom * t(k);
            v_abc_terminal(k,:) = V_peak_approx(k) * [cos(grid_angle), cos(grid_angle-2*pi/3), cos(grid_angle+2*pi/3)];
            i_abc(k,:) = I_peak_approx(k) * [cos(grid_angle), cos(grid_angle-2*pi/3), cos(grid_angle+2*pi/3)];
        end
    else % Droop and VSM
        id = x(:,1); iq = x(:,2);
        for k = 1:length(t)
           grid_angle = P.w_nom * t(k);
           T_inv_park = [cos(grid_angle), -sin(grid_angle);
                         cos(grid_angle - 2*pi/3), -sin(grid_angle - 2*pi/3);
                         cos(grid_angle + 2*pi/3), -sin(grid_angle + 2*pi/3)];
           
           % Calculate terminal voltage based on current and grid impedance
           v_d_terminal = P.V_peak - (P.R_grid * id(k) - P.w_nom * P.L_grid * iq(k));
           v_q_terminal = 0 - (P.R_grid * iq(k) + P.w_nom * P.L_grid * id(k));

           v_abc_terminal(k,:) = (T_inv_park * [v_d_terminal; v_q_terminal])';
           i_abc(k,:) = (T_inv_park * [id(k); iq(k)])';
        end
    end
end

% ============================================================================
% THE UNIFIED SYSTEM DYNAMICS FUNCTION
% ============================================================================
function dxdt = system_dynamics_dq_comparative(t, x, P, mode)
    % --- Power Ramp ---
    if t < P.t_ramp_end
        P_ref = (t / P.t_ramp_end) * P.P_final;
        Q_ref = (t / P.t_ramp_end) * P.Q_final;
    else
        P_ref = P.P_final;
        Q_ref = P.Q_final;
    end

    if strcmp(mode, 'VOC')
        delta = x(1); w_dev = x(2); p_filt = x(3);
        dp_filt_dt = (1/P.T_voc) * (P_ref - p_filt);
        V_grid = P.V_peak;
        Z_line = sqrt((P.R_filt + P.R_grid)^2 + (P.w_nom * (P.L_filt + P.L_grid))^2);
        asin_arg = p_filt * Z_line / V_grid^2;
        asin_arg = max(min(asin_arg, 1), -1);
        target_delta = asin(asin_arg);
        if isnan(target_delta), target_delta = delta; end
        dw_dev_dt = 20 * (target_delta - delta) - 2 * w_dev;
        ddelta_dt = w_dev;
        dxdt = [ddelta_dt; dw_dev_dt; dp_filt_dt];
        return;
    end
    
    id = x(1); iq = x(2); delta = x(3);
    w_dev = x(4); p_filt = x(5); q_filt = x(6);
    w_gfm = P.w_nom + w_dev;
    
    switch mode
        case 'Droop'
            w_ref = P.w_nom - P.mp_droop * (P_ref - p_filt);
            V_mag_ref = P.V_peak - P.mq_droop * (Q_ref - q_filt);
            dw_dev_dt = 10 * ((w_ref - P.w_nom) - w_dev);
        case 'VSM'
            dw_dev_dt = (1/P.J_vsm) * (P_ref - p_filt - P.D_vsm * w_dev);
            V_mag_ref = P.V_peak;
    end
    
    v_d_inv_ref = V_mag_ref; v_q_inv_ref = 0;
    v_d_inv = v_d_inv_ref * cos(delta) - v_q_inv_ref * sin(delta);
    v_q_inv = v_d_inv_ref * sin(delta) + v_q_inv_ref * cos(delta);
    v_d_grid = P.V_peak; v_q_grid = 0;
    
    % FIX: Include virtual inductance in total inductance for numerical stability
    L_total = P.L_filt + P.L_grid + P.L_virtual;
    R_total = P.R_filt + P.R_grid;
    
    di_d_dt = (1/L_total) * (-R_total * id + w_gfm * L_total * iq + v_d_inv - v_d_grid);
    di_q_dt = (1/L_total) * (-R_total * iq - w_gfm * L_total * id + v_q_inv - v_q_grid);
    
    ddelta_dt = w_dev;
    
    % Correctly calculate instantaneous power based on terminal voltage
    v_d_term = P.V_peak - (P.R_grid * id - P.w_nom * P.L_grid * iq);
    v_q_term = 0 - (P.R_grid * iq + P.w_nom * P.L_grid * id);
    
    P_inst = (3/2) * (v_d_term * id + v_q_term * iq);
    Q_inst = (3/2) * (v_q_term * id - v_d_term * iq);
    
    T_filter = P.T_filter_droop;
    if strcmp(mode, 'VSM'), T_filter = P.T_filter_vsm; end
    dp_filt_dt = (1/T_filter) * (P_inst - p_filt);
    dq_filt_dt = (1/T_filter) * (Q_inst - q_filt);
    
    dxdt = [di_d_dt; di_q_dt; ddelta_dt; dw_dev_dt; dp_filt_dt; dq_filt_dt];
end
