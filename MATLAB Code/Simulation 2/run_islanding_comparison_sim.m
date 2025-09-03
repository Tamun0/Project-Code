function run_islanding_comparison_sim()
    % Main function to run a fast, behavioral simulation comparing how GFL,
    % Droop, VSM, and VOC controllers respond to an islanding event.

    clear; clc; close all;

    % --- Simulation Case Selection ---
    % Control Mode: 'GFL', 'Droop', 'VSM', or 'VOC'
    control_mode = 'VSM'; 

    % ========================================================================
    % 1. DEFINE SYSTEM PARAMETERS (Struct 'P')
    % ========================================================================
    P.f_nom = 50;      % Nominal grid frequency (Hz)
    P.w_nom = 2 * pi * P.f_nom;
    
    P.P_load = 10000;  % The local load that remains after islanding (W)
    P.P_inv_ref = 10000; % The inverter's power setpoint before islanding
    
    % --- Controller Parameters ---
    % Droop
    P.mp_droop = (0.05 * P.w_nom) / P.P_load; % 5% frequency droop at full load
    % VSM
    P.J_vsm = 0.5; % Virtual Inertia
    P.D_vsm = 50.0; % Damping
    
    % --- Simulation Event Times ---
    P.t_island = 0.4;   % Time at which the grid is disconnected
    t_span = [0 1.5];
    
    disp(['Running ISLANDING COMPARISON simulation for ' control_mode '...']);

    % ========================================================================
    % 2. SET UP AND RUN THE SIMULATION
    % ========================================================================
    % x = [frequency_deviation, inverter_power_output]
    x0 = [0; P.P_inv_ref]; % Start synchronized, at reference power

    [t, x] = ode45(@(t,x) behavioral_dynamics_islanding(t, x, P, control_mode), t_span, x0);

    % ========================================================================
    % 3. POST-PROCESS AND PLOT RESULTS
    % ========================================================================
    disp('Simulation finished. Processing results for plotting...');
    
    frequency_hz = P.f_nom + x(:,1) / (2*pi);
    inverter_power_output_kw = x(:,2) / 1000;
    
    figure('Name', ['Islanding Response for ' control_mode ' control']);
    
    % Plot Frequency
    subplot(2,1,1);
    plot(t, frequency_hz, 'b', 'LineWidth', 2);
    hold on;
    line([P.t_island P.t_island], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    title(['Local Grid Frequency (' control_mode ')']);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    legend('Frequency', 'Islanding Event');
    grid on;

    % Plot Power Output
    subplot(2,1,2);
    plot(t, inverter_power_output_kw, 'b', 'LineWidth', 2);
    hold on;
    line([P.t_island P.t_island], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    title(['Inverter Active Power Output (' control_mode ')']);
    xlabel('Time (s)');
    ylabel('Power (kW)');
    legend('Inverter Power', 'Islanding Event');
    grid on;
end

% ============================================================================
% THE BEHAVIORAL DYNAMICS FUNCTION
% ============================================================================
function dxdt = behavioral_dynamics_islanding(t, x, P, control_mode)
    % Unpack state vector
    w_dev = x(1); % Frequency deviation from nominal (rad/s)
    P_out = x(2); % Current power output of the inverter (W)

    % Check if the grid is connected
    if t < P.t_island
        % --- Grid-Connected Behavior ---
        % The grid is strong and holds the frequency at nominal.
        dw_dev_dt = -100 * w_dev; % Drive frequency deviation to zero
        dP_out_dt = 10 * (P.P_inv_ref - P_out); % Ensure power is at setpoint
    else
        % --- Islanded Behavior ---
        power_mismatch = P_out - P.P_load;
        
        switch control_mode
            case 'GFL'
                % GFL inverter has lost its grid reference and fails.
                dw_dev_dt = -1000; % Frequency collapses
                dP_out_dt = -10 * P_out; % Power collapses
                
            case 'Droop'
                % Frequency is determined by the droop characteristic.
                target_w_dev = -P.mp_droop * (P_out - P.P_inv_ref);
                dw_dev_dt = 10 * (target_w_dev - w_dev);
                % Power output changes to match the load
                dP_out_dt = 10 * (P.P_load - P_out);
                
            case 'VSM'
                % Frequency changes based on the swing equation.
                dw_dev_dt = (1/P.J_vsm) * (-power_mismatch - P.D_vsm * w_dev);
                % Power output changes to match the load
                dP_out_dt = 10 * (P.P_load - P_out);
                
            case 'VOC'
                % VOC is modeled to quickly stabilize frequency and power.
                dw_dev_dt = -50 * w_dev; % Drive frequency deviation to zero quickly
                dP_out_dt = 50 * (P.P_load - P_out); % Match load very fast
        end
    end
    
    dxdt = [dw_dev_dt; dP_out_dt];
end
