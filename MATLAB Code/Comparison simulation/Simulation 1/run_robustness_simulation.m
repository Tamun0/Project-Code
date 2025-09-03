function run_robustness_simulation()
    % Main function to run a comparative simulation testing the robustness
    % of Droop, VSM, and VOC controllers under different grid conditions.
    % This version uses a more physically accurate and robust unified model.

    clear; clc; close all;

    % ========================================================================
    % --- Simulation Case Selection ---
    % ========================================================================
    % Control Mode: 'Droop', 'VSM', or 'VOC'
    control_mode = 'Droop'; 
    
    % Grid Type: 'StrongInductive', 'WeakInductive', or 'WeakResistive'
    grid_type = 'WeakInductive';

    % ========================================================================
    % 1. DEFINE SYSTEM PARAMETERS (Struct 'P')
    % ========================================================================
    P.f_nom = 50;
    P.w_nom = 2 * pi * P.f_nom;
    P.V_grid = 1.0; % Per-unit voltage of infinite bus
    P.V_terminal_ref = 1.0; % VSM terminal voltage reference
    
    % --- Power Setpoints ---
    P.P_initial = 0.5; % Start at 0.5 pu power
    P.P_step = 0.3;    % Step up by 0.3 pu
    
    % --- Controller Parameters ---
    % Droop
    P.mp_droop = 0.9 / 1.0; % 90% droop
    P.T_filter_droop = 0.1; 
    % VSM
    P.J_vsm = 2.0;  % Inertia constant (H) in seconds
    P.D_vsm = 25.0; % Damping coefficient
    % VOC
    P.T_voc = 0.02;

    % --- Grid Impedance Parameters ---
    switch grid_type
        case 'StrongInductive'
            P.R_line = 0.02; P.X_line = 0.2;
        case 'WeakInductive'
            P.R_line = 0.05; P.X_line = 0.5; 
        case 'WeakResistive'
            P.R_line = 0.1; P.X_line = 0.1; 
    end
    
    % --- Simulation Setup ---
    P.t_disturbance = 1.0;
    t_span = [0 5.0];
    
    disp(['Running ROBUSTNESS test for ' control_mode ' in a ' grid_type ' grid...']);

    % ========================================================================
    % 2. SET UP AND RUN THE SIMULATION
    % ========================================================================
    % x = [delta, w_dev]
    % delta: internal angle, w_dev: freq deviation
    
    Z_line = sqrt(P.R_line^2 + P.X_line^2);
    phi_line = atan(P.X_line / P.R_line);
    delta0 = asin(P.P_initial * Z_line / (P.V_terminal_ref * P.V_grid)) + phi_line;
    if isnan(delta0), delta0 = 0; end
    x0 = [delta0; 0];

    [t, x] = ode23t(@(t,x) unified_dynamics(t, x, P, control_mode), t_span, x0);

    % ========================================================================
    % 3. POST-PROCESS AND PLOT RESULTS
    % ========================================================================
    disp('Simulation finished. Processing results...');
    
    delta_rad = x(:,1);
    
    % Calculate electrical power, terminal voltage, and current
    P_elec_pu = (P.V_terminal_ref * P.V_grid / Z_line) * cos(phi_line - delta_rad) - (P.V_grid^2 / Z_line)*cos(phi_line);
    I_approx = P_elec_pu ./ P.V_grid; 
    V_term_pu = sqrt((P.V_grid*cos(-delta_rad) + I_approx.*P.R_line).^2 + (P.V_grid*sin(-delta_rad) + I_approx.*P.X_line).^2);

    
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

    % Determine the power reference based on time
    if t < P.t_disturbance
        P_ref = P.P_initial;
    else
        P_ref = P.P_initial + P.P_step;
    end
    
    % More accurate electrical power calculation for R-X lines
    Z_line = sqrt(P.R_line^2 + P.X_line^2);
    phi_line = atan(P.X_line / P.R_line);
    P_elec = (P.V_terminal_ref * P.V_grid / Z_line) * cos(phi_line - delta) - (P.V_grid^2 / Z_line)*cos(phi_line);
    
    % --- Controller Dynamics ---
    switch control_mode
        case 'Droop'
            % Droop control tries to match power P_ref based on freq deviation
            P_mech = P_elec - (w_dev / P.mp_droop);
            dw_dev_dt = 10 * (P_ref - P_mech); % Simplified first-order response
            
        case 'VSM'
            % VSM uses the swing equation directly
            P_mech = P_ref;
            dw_dev_dt = (1 / (2*P.J_vsm)) * (P_mech - P_elec - P.D_vsm * w_dev);
            
        case 'VOC'
            % VOC is modeled as a fast power-tracking response
            P_mech = P_ref;
            % Angle is driven to match the power reference
            target_delta = phi_line - acos( (P_ref + (P.V_grid^2 / Z_line)*cos(phi_line)) * Z_line / (P.V_terminal_ref*P.V_grid));
            if ~isreal(target_delta), target_delta = delta; end
            dw_dev_dt = 50 * (target_delta - delta) - 5 * w_dev; % Fast damped response to target angle
    end
    
    ddelta_dt = w_dev * P.w_nom; % Convert pu freq dev to rad/s for angle

    dxdt = [ddelta_dt; dw_dev_dt];
end
