function run_vsm_simulation_1()
    % Main function to run a simulation highlighting the inertial and
    % damping properties of a Virtual Synchronous Machine (VSM).
    % A VSM (Inverter 1) and a simple droop GFM (Inverter 2) supply a
    % common load. A disturbance is applied to see how the VSM stabilizes the system.

    clear; clc; close all;

    % ========================================================================
    % 1. DEFINE SYSTEM PARAMETERS (Struct 'P')
    % ========================================================================
    P.f_nom = 50;
    P.w_nom = 2 * pi * P.f_nom;
    
    % --- Power & Load Parameters ---
    P.P1_ref = 10000; % VSM power setpoint
    P.P2_ref = 10000; % Droop GFM power setpoint
    P.P_load_initial = P.P1_ref + P.P2_ref; % Total initial load
    P.P_load_step = -5000; % Disturbance: a sudden 5kW load reduction
    
    % --- VSM Control Parameters (Inverter 1) ---
    % From paper: J is inertia, D is damping
    P.J_vsm = 0.5;  % Virtual Inertia Constant (s)
    P.D_vsm = 50.0; % Damping Coefficient (Increased for stability)
    
    % --- Simple Droop GFM Parameters (Inverter 2) ---
    P.mp2 = (0.05 * P.w_nom) / P.P2_ref; % 5% droop
    P.T_droop = 0.1; % Time constant for droop response
    
    % --- Simulation Event Times ---
    P.t_disturbance = 0.8; % Time of the load step disturbance
    t_span = [0 2.0];
    
    disp('Running VSM INERTIAL RESPONSE simulation...');

    % ========================================================================
    % 2. SET UP AND RUN THE SIMULATION
    % ========================================================================
    % x = [w_dev, P1_out, P2_out]
    % (system frequency deviation, and power output for each inverter)
    x0 = [0; P.P1_ref; P.P2_ref]; 

    [t, x] = ode23t(@(t,x) vsm_dynamics(t, x, P), t_span, x0);

    % ========================================================================
    % 3. POST-PROCESS AND PLOT RESULTS
    % ========================================================================
    disp('Simulation finished. Processing results for plotting...');
    
    frequency_hz = P.f_nom + x(:,1) / (2*pi); 
    inv1_power_kw = x(:,2) / 1000;
    inv2_power_kw = x(:,3) / 1000;
    
    figure('Name', 'VSM vs. Droop GFM Response to Disturbance');
    
    % Plot Frequency
    subplot(2,1,1);
    plot(t, frequency_hz, 'b', 'LineWidth', 2);
    hold on;
    line([P.t_disturbance P.t_disturbance], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    title('System Frequency');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    legend('Frequency', 'Load Disturbance');
    grid on;

    % Plot Power Outputs
    subplot(2,1,2);
    plot(t, inv1_power_kw, 'm', 'LineWidth', 2);
    hold on;
    plot(t, inv2_power_kw, 'c', 'LineWidth', 2);
    line([P.t_disturbance P.t_disturbance], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    title('Inverter Active Power Outputs');
    xlabel('Time (s)');
    ylabel('Power (kW)');
    legend('Inverter 1 (VSM)', 'Inverter 2 (Droop)', 'Load Disturbance');
    grid on;
end

% ============================================================================
% THE BEHAVIORAL DYNAMICS FUNCTION FOR VSM
% ============================================================================
function dxdt = vsm_dynamics(t, x, P)
    % Unpack state vector
    w_dev = x(1);
    P1_out = x(2); % VSM Power Output
    P2_out = x(3); % Droop GFM Power Output

    % Determine the total load on the system based on time
    if t < P.t_disturbance
        P_total_load = P.P_load_initial;
    else
        P_total_load = P.P_load_initial + P.P_load_step;
    end
    
    % Total power mismatch in the islanded system
    power_mismatch = (P1_out + P2_out) - P_total_load;
    
    % --- System Frequency Dynamics ---
    % The rate of change of frequency is based on the total power mismatch
    % acting on the VSM's virtual inertia. The simple droop inverter has no inertia.
    % J*dw/dt = -power_mismatch
    dw_dev_dt = (1/P.J_vsm) * (-power_mismatch);
    
    % --- Inverter 1: VSM Power Response ---
    % Its power output tries to match its internal reference, which includes
    % the damping term D*w_dev.
    % This is the swing equation: P_out = P_ref - D*w_dev
    target_power_1 = P.P1_ref - P.D_vsm * w_dev;
    dP1_out_dt = 25 * (target_power_1 - P1_out); % First-order response toward target

    % --- Inverter 2: Simple Droop GFM Power Response ---
    % Its power output only responds to the frequency deviation based on
    % its droop gain, with a first-order delay.
    target_power_2 = P.P2_ref - w_dev / P.mp2;
    dP2_out_dt = (1/P.T_droop) * (target_power_2 - P2_out);

    dxdt = [dw_dev_dt; dP1_out_dt; dP2_out_dt];
end
