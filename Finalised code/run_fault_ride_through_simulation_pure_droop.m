function run_fault_ride_through_simulation_pure_droop()
    % Main function to run a comparative simulation of fault ride-through
    % performance. This version includes a "pure" Droop model alongside
    % the VSM and behavioral VOC models for a complete comparison.

    clear; clc; close all;

    % ========================================================================
    % --- Simulation Case Selection ---
    % ========================================================================
    % Control Mode: 'Droop', 'VSM', or 'VOC'
    control_mode = 'VOC'; 
    
    % ========================================================================
    % 1. DEFINE SYSTEM PARAMETERS (Struct 'P')
    % ========================================================================
    P.f_nom = 50;
    P.w_nom = 2 * pi * P.f_nom;
    P.V_grid_nominal = 1.0; % Nominal per-unit voltage of infinite bus
    
    % --- Power Setpoint ---
    P.P_ref = 0.8; % Operating at 0.8 pu power
    
    % --- Controller Parameters ---
    % VSM - for comparison
    P.J_vsm = 0.2; P.D_vsm = 50.0; 
    % Pure Droop
    P.mp_droop = (0.05 * P.w_nom) / P.P_ref; % 5% droop at rated power
    P.T_filter_droop = 0.05; 
    % VSM (High-inertia model)
    %P.J_vsm = 5.0;  P.D_vsm = 50.0; P.T_filter_vsm = 0.02;
    % VOC (Behavioral model time constant)
    P.T_voc = 0.02;

    % --- Grid Impedance (Weak Inductive) ---
    P.R_line = 0.1; 
    P.X_line = 1.0; 
    
    % --- Fault Parameters ---
    P.t_fault_start = 1.0;
    P.fault_duration = 0.15;
    P.t_fault_clear = P.t_fault_start + P.fault_duration;
    P.V_fault = 0.1;
    
    % --- Simulation Setup ---
    t_end = 5.0;
    
    disp(['Running FAULT RIDE-THROUGH test for ' control_mode '...']);

    % ========================================================================
    % 2. SET UP AND RUN THE SIMULATION
    % ========================================================================
    Z_line = sqrt(P.R_line^2 + P.X_line^2);
    delta0 = asin(P.P_ref * Z_line / P.V_grid_nominal);
    if isnan(delta0), delta0 = 0; end
    
    % Initial conditions depend on the model
    if strcmp(control_mode, 'Droop')
        % State for Pure Droop: x = [delta, p_filt]
        x0 = [delta0; P.P_ref];
    else
        % State for VSM-based models: x = [delta, w_dev, p_filt]
        x0 = [delta0; 0; P.P_ref];
    end

    options = odeset('RelTol', 1e-4, 'AbsTol', 1e-5);
    
    % Run in three parts: pre-fault, during-fault, post-fault
    [t1, x1] = ode23t(@(t,x) fault_dynamics(t, x, P, control_mode, 0), [0 P.t_fault_start], x0, options);
    [t2, x2] = ode23t(@(t,x) fault_dynamics(t, x, P, control_mode, 1), [P.t_fault_start P.t_fault_clear], x1(end,:)', options);
    [t3, x3] = ode23t(@(t,x) fault_dynamics(t, x, P, control_mode, 0), [P.t_fault_clear t_end], x2(end,:)', options);
    
    t = [t1; t2(2:end); t3(2:end)];
    x = [x1; x2(2:end,:); x3(2:end,:)];

    % ========================================================================
    % 3. POST-PROCESS AND PLOT RESULTS
    % ========================================================================
    disp('Simulation finished. Processing results...');
    
    if strcmp(control_mode, 'Droop')
        delta_rad = x(:,1);
        p_filt = x(:,2);
        w_dev = -P.mp_droop * (p_filt - P.P_ref); % Frequency is an algebraic result
    else
        delta_rad = x(:,1);
        w_dev = x(:,2);
    end
    
    V_grid_vec = (t < P.t_fault_start | t > P.t_fault_clear) * P.V_grid_nominal + ...
                 (t >= P.t_fault_start & t <= P.t_fault_clear) * P.V_fault;
    P_elec_pu = (V_grid_vec / Z_line) .* sin(delta_rad);
    
    figure('Name', [control_mode ' Fault Ride-Through Response']);
    sgtitle([control_mode ' Fault Ride-Through Response']);
    subplot(3,1,1); plot(t, P.f_nom + w_dev / (2*pi), 'b', 'LineWidth', 2);
    hold on; line([P.t_fault_start P.t_fault_start], ylim, 'Color', 'r', 'LineStyle', '--');
    line([P.t_fault_clear P.t_fault_clear], ylim, 'Color', 'k', 'LineStyle', '--');
    title('Inverter Frequency'); xlabel('Time (s)'); ylabel('Frequency (Hz)'); grid on;

    subplot(3,1,2); plot(t, P_elec_pu, 'm', 'LineWidth', 2);
    hold on; line([P.t_fault_start P.t_fault_start], ylim, 'Color', 'r', 'LineStyle', '--');
    line([P.t_fault_clear P.t_fault_clear], ylim, 'Color', 'k', 'LineStyle', '--');
    title('Inverter Active Power Output'); xlabel('Time (s)'); ylabel('Power (pu)'); grid on;

    subplot(3,1,3); plot(t, rad2deg(unwrap(delta_rad)), 'g', 'LineWidth', 2);
    hold on; line([P.t_fault_start P.t_fault_start], ylim, 'Color', 'r', 'LineStyle', '--');
    line([P.t_fault_clear P.t_fault_clear], ylim, 'Color', 'k', 'LineStyle', '--');
    title('Inverter Internal Angle (Rotor Angle)'); xlabel('Time (s)'); ylabel('Angle (degrees)'); grid on;
end

% ============================================================================
% THE DYNAMICS FUNCTION FOR FAULT SIMULATION
% ============================================================================
function dxdt = fault_dynamics(t, x, P, control_mode, fault_status)
    % Determine grid voltage based on fault status
    if fault_status == 1, V_grid_dynamic = P.V_fault;
    else, V_grid_dynamic = P.V_grid_nominal; end
    
    Z_line = sqrt(P.R_line^2 + P.X_line^2);
    
    % --- Controller Dynamics ---
    switch control_mode
        case 'Droop'
            % Unpack state vector: x = [delta, p_filt]
            delta = x(1);
            p_filt = x(2);
            
            P_elec = (V_grid_dynamic / Z_line) * sin(delta);
            
            % Power filter dynamics
            dp_filt_dt = (1/P.T_filter_droop) * (P_elec - p_filt);
            
            % Frequency is calculated algebraically from the droop law
            w_dev = -P.mp_droop * (p_filt - P.P_ref);
            
            % Rate of change of angle is the frequency deviation
            ddelta_dt = w_dev;
            
            dxdt = [ddelta_dt; dp_filt_dt];

        case 'VSM'
            delta = x(1); w_dev = x(2); p_filt = x(3);
            P_elec = (V_grid_dynamic / Z_line) * sin(delta);
            dp_filt_dt = (1/P.T_filter_droop) * (P_elec - p_filt);
            dw_dev_dt = (1/P.J_vsm) * (P.P_ref - p_filt - P.D_vsm * w_dev);
            ddelta_dt = w_dev;
            dxdt = [ddelta_dt; dw_dev_dt; dp_filt_dt];
            
        % case 'VSM'
        %     delta = x(1); w_dev = x(2); p_filt = x(3);
        %     P_elec = (V_grid_dynamic / Z_line) * sin(delta);
        %     dp_filt_dt = (1/P.T_filter_vsm) * (P_elec - p_filt);
        %     dw_dev_dt = (1/P.J_vsm) * (P.P_ref - p_filt - P.D_vsm * w_dev);
        %     ddelta_dt = w_dev;
        %     dxdt = [ddelta_dt; dw_dev_dt; dp_filt_dt];
            
        case 'VOC'
            delta = x(1); w_dev = x(2); p_filt = x(3);
            P_elec = (V_grid_dynamic / Z_line) * sin(delta);
            dp_filt_dt = (1/P.T_voc) * (P.P_ref - p_filt);
            asin_arg = p_filt * Z_line / V_grid_dynamic;
            if asin_arg > 1.0, asin_arg = 1.0;
            elseif asin_arg < -1.0, asin_arg = -1.0; end
            target_delta = asin(asin_arg);
            if isnan(target_delta), target_delta = delta; end
            dw_dev_dt = 20 * (target_delta - delta) - 2 * w_dev;
            ddelta_dt = w_dev;
            dxdt = [ddelta_dt; dw_dev_dt; dp_filt_dt];
    end
end
