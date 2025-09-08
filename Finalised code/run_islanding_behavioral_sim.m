function run_islanding_behavioral_sim()
    % Main function to run a fast, behavioral simulation of GFM vs GFL islanding.
    % This script uses a simplified model focusing on controller dynamics to
    % run very quickly while demonstrating the core concept.

    clear; clc; close all;

    % --- Simulation Case Selection ---
    control_mode = 'GFM'; % CHANGE THIS to 'GFM' or 'GFL'

    % ========================================================================
    % 1. DEFINE SYSTEM PARAMETERS (Struct 'P')
    % ========================================================================
    P.f_nom = 50;      % Nominal grid frequency (Hz)
    P.w_nom = 2 * pi * P.f_nom;
    
    P.P_load = 10000;  % The local load connected to the inverter (W)
    P.P_inv_ref = 10000; % The inverter's power setpoint
    
    % --- GFM Droop Gain ---
    P.mp = (0.05 * P.w_nom) / P.P_load; % 5% frequency droop at full load
    
    % --- Simulation Event Times ---
    P.t_island = 0.4;   % Time at which the grid is disconnected
    t_span = [0 0.8];
    
    disp(['Running ISLANDING simulation for ' control_mode '...']);

    % ========================================================================
    % 2. SET UP AND RUN THE SIMULATION
    % ========================================================================
    % x = [frequency_deviation]
    x0 = 0; % Start synchronized with the grid (zero frequency deviation)

    % Use a standard solver, as the system is no longer stiff
    [t, x] = ode45(@(t,x) behavioral_dynamics(t, x, P, control_mode), t_span, x0);

    % ========================================================================
    % 3. POST-PROCESS AND PLOT RESULTS
    % ========================================================================
    disp('Simulation finished. Processing results for plotting...');
    
    frequency_hz = P.f_nom + x(:,1) / (2*pi);
    inverter_power_output = zeros(size(t));
    
    for k = 1:length(t)
        if strcmp(control_mode, 'GFM')
            % In GFM mode, power is determined by the load after islanding
            if t(k) < P.t_island
                inverter_power_output(k) = P.P_inv_ref;
            else
                inverter_power_output(k) = P.P_load;
            end
        else % GFL mode
            % In GFL mode, it injects reference power until islanded, then fails
            if t(k) < P.t_island
                inverter_power_output(k) = P.P_inv_ref;
            else
                inverter_power_output(k) = 0; % Output collapses
            end
        end
    end

    figure('Name', ['Behavioral Islanding Response for ' control_mode ' control']);
    sgtitle(['Islanding Response for ' control_mode ' control']);
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
    plot(t, inverter_power_output / 1000, 'b', 'LineWidth', 2);
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
function dxdt = behavioral_dynamics(t, x, P, control_mode)
    % Unpack state vector
    w_dev = x(1); % Frequency deviation from nominal (w - w_nom)

    % Check if the grid is connected
    if t < P.t_island
        % --- Grid-Connected Behavior ---
        % The grid is strong and holds the frequency at nominal.
        % The deviation is driven towards zero.
        dw_dev_dt = -100 * w_dev; % A simple restoring force to keep it at 0
    else
        % --- Islanded Behavior ---
        if strcmp(control_mode, 'GFM')
            % GFM inverter is now supplying the load. Its frequency is
            % determined by the power mismatch and its droop gain.
            % Power mismatch = P_inverter - P_load
            % From droop: w_dev = -mp * (P_inverter - P_ref)
            % Let P_inverter = P_load
            power_mismatch = P.P_inv_ref - P.P_load;
            
            % The frequency will change until the droop equation is satisfied
            target_w_dev = -P.mp * power_mismatch;
            dw_dev_dt = 10 * (target_w_dev - w_dev); % Drive towards new steady state
            
        else % 'GFL'
            % The GFL inverter has lost its grid reference. Its frequency
            % becomes undefined, and it stops producing power.
            % We model this as the frequency collapsing.
            dw_dev_dt = -1000; % Rapidly drive frequency to a failed state
        end
    end
    
    dxdt = dw_dev_dt;
end
