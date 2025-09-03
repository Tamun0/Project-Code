function run_fault_ride_through_simulation()
    % Main function to run a comparative simulation of fault ride-through
    % performance for Droop, VSM, and VOC controllers.

    clear; clc; close all;

    % ========================================================================
    % --- Simulation Case Selection ---
    % ========================================================================
    % Control Mode: 'Droop', 'VSM', or 'VOC'
    control_mode = 'VOC'; 
    
    % Grid Type: Using a Weak Inductive grid for a challenging but stable baseline
    grid_type = 'WeakInductive';

    % ========================================================================
    % 1. DEFINE SYSTEM PARAMETERS (Struct 'P')
    % ========================================================================
    P.f_nom = 50;
    P.w_nom = 2 * pi * P.f_nom;
    P.V_grid_nominal = 1.0; % Nominal per-unit voltage of infinite bus
    
    % --- Power Setpoint ---
    P.P_ref = 0.8; % Operating at 0.8 pu power
    
    % --- Controller Parameters (from robust model) ---
    P.J_droop = 0.2; P.D_droop = 5.0; P.T_filter_droop = 0.05; 
    P.J_vsm = 5.0;  P.D_vsm = 1.5; P.T_filter_vsm = 0.02;
    P.T_voc = 0.02;

    % --- Grid Impedance ---
    P.R_line = 0.1; P.X_line = 1.0; 
    
    % --- Fault Parameters ---
    P.t_fault_start = 1.0;       % Fault starts at 1.0s
    P.fault_duration = 0.15;     % Fault lasts 150ms
    P.t_fault_clear = P.t_fault_start + P.fault_duration;
    P.V_fault = 0.1;             % Voltage sags to 0.1 pu during fault
    
    % --- Simulation Setup ---
    t_end = 5.0;
    
    disp(['Running FAULT RIDE-THROUGH test for ' control_mode '...']);

    % ========================================================================
    % 2. SET UP AND RUN THE SIMULATION (IN THREE PARTS)
    % ========================================================================
    % x = [delta, w_dev, p_filt]
    Z_line = sqrt(P.R_line^2 + P.X_line^2);
    delta0 = asin(P.P_ref * Z_line / P.V_grid_nominal);
    if isnan(delta0), delta0 = 0; end
    x0 = [delta0; 0; P.P_ref];

    options = odeset('RelTol', 1e-4, 'AbsTol', 1e-5);
    
    % Part 1: Pre-fault
    t_span1 = [0 P.t_fault_start];
    [t1, x1] = ode23t(@(t,x) fault_dynamics(t, x, P, control_mode, 0), t_span1, x0, options);
    
    % Part 2: During-fault
    x0_part2 = x1(end, :)';
    t_span2 = [P.t_fault_start P.t_fault_clear];
    [t2, x2] = ode23t(@(t,x) fault_dynamics(t, x, P, control_mode, 1), t_span2, x0_part2, options);
    
    % Part 3: Post-fault (Recovery)
    x0_part3 = x2(end, :)';
    t_span3 = [P.t_fault_clear t_end];
    [t3, x3] = ode23t(@(t,x) fault_dynamics(t, x, P, control_mode, 0), t_span3, x0_part3, options);
    
    % Combine results
    t = [t1; t2(2:end); t3(2:end)];
    x = [x1; x2(2:end,:); x3(2:end,:)];

    % ========================================================================
    % 3. POST-PROCESS AND PLOT RESULTS
    % ========================================================================
    disp('Simulation finished. Processing results...');
    
    delta_rad = x(:,1);
    w_dev_pu = x(:,2) / P.w_nom; % Convert rad/s deviation to pu
    
    % Recalculate power for plotting, accounting for changing grid voltage
    V_grid_vec = (t < P.t_fault_start | t > P.t_fault_clear) * P.V_grid_nominal + ...
                 (t >= P.t_fault_start & t <= P.t_fault_clear) * P.V_fault;
    P_elec_pu = (V_grid_vec / Z_line) .* sin(delta_rad);
    
    figure('Name', [control_mode ' Fault Ride-Through Response']);
    
    subplot(3,1,1);
    plot(t, P.f_nom * (1 + w_dev_pu), 'b', 'LineWidth', 2);
    hold on;
    line([P.t_fault_start P.t_fault_start], ylim, 'Color', 'r', 'LineStyle', '--');
    line([P.t_fault_clear P.t_fault_clear], ylim, 'Color', 'k', 'LineStyle', '--');
    title('Inverter Frequency');
    xlabel('Time (s)'); ylabel('Frequency (Hz)'); grid on;

    subplot(3,1,2);
    plot(t, P_elec_pu, 'm', 'LineWidth', 2);
    hold on;
    line([P.t_fault_start P.t_fault_start], ylim, 'Color', 'r', 'LineStyle', '--');
    line([P.t_fault_clear P.t_fault_clear], ylim, 'Color', 'k', 'LineStyle', '--');
    title('Inverter Active Power Output');
    xlabel('Time (s)'); ylabel('Power (pu)'); grid on;

    subplot(3,1,3);
    plot(t, rad2deg(unwrap(delta_rad)), 'g', 'LineWidth', 2);
    hold on;
    line([P.t_fault_start P.t_fault_start], ylim, 'Color', 'r', 'LineStyle', '--');
    line([P.t_fault_clear P.t_fault_clear], ylim, 'Color', 'k', 'LineStyle', '--');
    title('Inverter Internal Angle (Rotor Angle)');
    xlabel('Time (s)'); ylabel('Angle (degrees)'); grid on;
end

% ============================================================================
% THE DYNAMICS FUNCTION FOR FAULT SIMULATION
% ============================================================================
function dxdt = fault_dynamics(t, x, P, control_mode, fault_status)
    % Unpack state vector
    delta = x(1);
    w_dev = x(2);
    p_filt = x(3);

    % Determine grid voltage based on fault status
    if fault_status == 1
        V_grid_dynamic = P.V_fault;
    else
        V_grid_dynamic = P.V_grid_nominal;
    end
    
    % Electrical power based on current power-angle relationship
    Z_line = sqrt(P.R_line^2 + P.X_line^2);
    P_elec = (V_grid_dynamic / Z_line) * sin(delta);
    
    % --- Controller Dynamics ---
    switch control_mode
        case 'Droop'
            dp_filt_dt = (1/P.T_filter_droop) * (P_elec - p_filt);
            dw_dev_dt = (1/P.J_droop) * (P.P_ref - p_filt - P.D_droop * w_dev);
            
        case 'VSM'
            dp_filt_dt = (1/P.T_filter_vsm) * (P_elec - p_filt);
            dw_dev_dt = (1/P.J_vsm) * (P.P_ref - p_filt - P.D_vsm * w_dev);
            
        case 'VOC'
            dp_filt_dt = (1/P.T_voc) * (P.P_ref - p_filt);
            target_delta = asin(p_filt * Z_line / V_grid_dynamic);
            if isnan(target_delta), target_delta = delta; end
            dw_dev_dt = 20 * (target_delta - delta) - 2 * w_dev; 
    end
    
    ddelta_dt = w_dev;

    dxdt = [ddelta_dt; dw_dev_dt; dp_filt_dt];
end
