function run_robustness_simulation2()
   % Main function to run a comparative simulation testing the robustness
    % of Droop, VSM, and VOC controllers under different grid conditions.

    clear; clc; close all;

    % ========================================================================
    % --- Simulation Case Selection ---
    % ========================================================================
    % Control Mode: 'Droop', 'VSM', or 'VOC'
    control_mode = 'Droop'; 
    
    % Grid Type: 'StrongInductive', 'WeakInductive', or 'WeakResistive'
    grid_type = 'WeakResistive';

    % ========================================================================
    % 1. DEFINE SYSTEM PARAMETERS (Struct 'P')
    % ========================================================================
    P.f_nom = 50;
    P.w_nom = 2 * pi * P.f_nom;
    P.V_grid = 1.0; % Per-unit voltage of infinite bus
    
    % --- Power Setpoints ---
    P.P_initial = 0.5; % Start at 0.5 pu power
    P.P_step = 0.3;    % Step up by 0.3 pu
    
    % --- Controller Parameters ---
    % Droop
    P.mp_droop = (0.05 * P.w_nom) / 1.0; 
    P.T_filter_droop = 0.1; % Filter time constant
    % VSM
    P.J_vsm = 0.2;  
    P.D_vsm = 1.5; % Damping coefficient
    P.T_filter_vsm = 0.02;
    % VOC
    P.T_voc = 0.02; % Fast time constant for first-order response

    % --- Grid Impedance Parameters ---
    switch grid_type
        case 'StrongInductive'
            P.R_line = 0.05; P.X_line = 0.5;
        case 'WeakInductive'
            P.R_line = 0.1; P.X_line = 1.0; 
        case 'WeakResistive'
            P.R_line = 0.8; P.X_line = 0.8; 
    end
    
    % --- Simulation Setup ---
    P.t_disturbance = 1.0;
    t_span = [0 5.0];
    
    disp(['Running ROBUSTNESS test for ' control_mode ' in a ' grid_type ' grid...']);

    % ========================================================================
    % 2. SET UP AND RUN THE SIMULATION
    % ========================================================================
    % x = [delta, w_dev, p_filt]
    % delta: internal angle, w_dev: freq deviation, p_filt: filtered power
    
    Z_line = sqrt(P.R_line^2 + P.X_line^2);
    delta0 = asin(P.P_initial * Z_line / P.V_grid);
    if isnan(delta0), delta0 = 0; end
    x0 = [delta0; 0; P.P_initial];

    [t, x] = ode23t(@(t,x) unified_dynamics(t, x, P, control_mode), t_span, x0);

    % ========================================================================
    % 3. POST-PROCESS AND PLOT RESULTS
    % ========================================================================
    disp('Simulation finished. Processing results...');
    
    delta_rad = x(:,1);
    
    % Calculate electrical power and terminal voltage
    P_elec_pu = (P.V_grid / Z_line) * sin(delta_rad);
    I_approx = P_elec_pu ./ P.V_grid; 
    V_term_pu = sqrt((P.V_grid - I_approx.*P.R_line).^2 + (I_approx.*P.X_line).^2);
    
    figure('Name', [control_mode ' Response in ' grid_type ' Grid']);
    
    subplot(3,1,1);
    plot(t, V_term_pu, 'b', 'LineWidth', 2);
    hold on;
    line([P.t_disturbance P.t_disturbance], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    title('Inverter Terminal Voltage');
    xlabel('Time (s)'); ylabel('Voltage (pu)'); grid on;

    subplot(3,1,2);
    plot(t, P_elec_pu, 'm', 'LineWidth', 2);
    hold on;
    line([P.t_disturbance P.t_disturbance], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    title('Inverter Active Power Output');
    xlabel('Time (s)'); ylabel('Power (pu)'); grid on;

    subplot(3,1,3);
    plot(t, I_approx, 'g', 'LineWidth', 2);
    hold on;
    line([P.t_disturbance P.t_disturbance], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    title('Inverter Output Current Magnitude');
    xlabel('Time (s)'); ylabel('Current (pu)'); grid on;
end

% ============================================================================
% THE UNIFIED BEHAVIORAL DYNAMICS FUNCTION
% ============================================================================
function dxdt = unified_dynamics(t, x, P, control_mode)
    % Unpack state vector
    delta = x(1);
    w_dev = x(2);
    p_filt = x(3);

    % Determine the power reference based on time
    if t < P.t_disturbance
        P_ref = P.P_initial;
    else
        P_ref = P.P_initial + P.P_step;
    end
    
    % Electrical power based on current power-angle relationship
    Z_line = sqrt(P.R_line^2 + P.X_line^2);
    P_elec = (P.V_grid / Z_line) * sin(delta);
    
    % --- Controller Dynamics ---
    switch control_mode
        case 'Droop'
            dp_filt_dt = (1/P.T_filter_droop) * (P_elec - p_filt);
            target_w_dev = -P.mp_droop * (p_filt - P_ref);
            dw_dev_dt = 10 * (target_w_dev - w_dev);
            
        case 'VSM'
            dp_filt_dt = (1/P.T_filter_vsm) * (P_elec - p_filt);
            dw_dev_dt = (1/P.J_vsm) * (P_ref - p_filt - P.D_vsm * w_dev);
            
        case 'VOC'
            % FIX: Use the stable, first-order behavioral model for VOC
            % This models its fast, non-overshooting response directly.
            dp_filt_dt = (1/P.T_voc) * (P_ref - p_filt);
            
            % The angle must follow the power to maintain the power-angle relationship
            target_delta = asin(p_filt * Z_line / P.V_grid);
            if isnan(target_delta), target_delta = delta; end
            
            % The frequency deviation is the error that drives the angle to its target
            dw_dev_dt = 20 * (target_delta - delta) - 2 * w_dev; % Damped response to target angle
    end
    
    ddelta_dt = w_dev;

    dxdt = [ddelta_dt; dw_dev_dt; dp_filt_dt];
end
