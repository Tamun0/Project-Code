function run_vsm_smib_detailed_simulation()
    % Main function to run a detailed and robust simulation of a VSM
    % connected to an infinite bus (SMIB). This model is based on the
    % detailed VSM structure presented in the provided papers to accurately
    % highlight the inertial and damping characteristics.

    clear; clc; close all;

    % ========================================================================
    % 1. DEFINE SYSTEM PARAMETERS (Struct 'P')
    % ========================================================================
    P.f_nom = 50;      % Nominal grid frequency (Hz)
    P.w_nom = 2 * pi * P.f_nom;
    
    P.V_terminal = 1.0; % Per-unit voltage at inverter terminal
    P.V_grid = 1.0;     % Per-unit voltage of infinite bus
    
    % --- Transmission Line ---
    P.X_line = 0.5; % Per-unit impedance of the line connecting VSM to grid
    
    % --- Power Setpoints ---
    P.P_initial = 0.5; % Start at 0.5 pu power
    P.P_step = 0.3;    % Step up by 0.3 pu
    
    % --- VSM Control Parameters (based on papers) ---
    P.Ta = 5.0;  % Virtual mechanical time constant (J or 2H) in seconds
    P.kd = 20.0; % Damping coefficient
    P.kw = 20.0; % Droop gain
    
    % --- Power Measurement Filter ---
    P.T_filter = 0.02; % Time constant for the power measurement filter

    % --- Simulation Event Times ---
    P.t_disturbance = 1.0; % Time of the power reference step
    t_span = [0 10.0];     % Simulate for 10 seconds to see oscillations
    
    disp('Running VSM (SMIB) DETAILED SIMULATION...');

    % ========================================================================
    % 2. SET UP AND RUN THE SIMULATION
    % ========================================================================
    % x = [delta, w_vsm, p_meas]
    % delta: the VSM's internal rotor angle (rad)
    % w_vsm: the VSM's internal frequency (pu)
    % p_meas: the filtered measurement of electrical power (pu)
    
    delta0 = asin(P.P_initial * P.X_line / (P.V_terminal * P.V_grid));
    x0 = [delta0; 1.0; P.P_initial]; % Start at steady state

    [t, x] = ode23t(@(t,x) vsm_detailed_dynamics(t, x, P), t_span, x0);

    % ========================================================================
    % 3. POST-PROCESS AND PLOT RESULTS
    % ========================================================================
    disp('Simulation finished. Processing results for plotting...');
    
    delta_rad = x(:,1);
    w_vsm_pu = x(:,2);
    p_meas_pu = x(:,3);
    
    % Calculate other quantities for plotting
    frequency_hz = w_vsm_pu * P.f_nom;
    electrical_power_pu = (P.V_terminal * P.V_grid / P.X_line) * sin(delta_rad);
    
    figure('Name', 'VSM Response to Power Step Command');
    
    % Plot Frequency
    sgtitle('VSM Response to Power Step Command')
    subplot(3,1,1);
    plot(t, frequency_hz, 'b', 'LineWidth', 2);
    hold on;
    line([P.t_disturbance P.t_disturbance], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    title('VSM Frequency');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    grid on;

    % Plot Power Output
    subplot(3,1,2);
    plot(t, electrical_power_pu, 'm', 'LineWidth', 2);
    hold on;
    plot(t, p_meas_pu, 'k:', 'LineWidth', 1.5);
    line([P.t_disturbance P.t_disturbance], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    title('VSM Electrical Power Output');
    xlabel('Time (s)');
    ylabel('Power (pu)');
    legend('P_{elec}', 'P_{meas}', 'Disturbance');
    grid on;

    % Plot Rotor Angle
    subplot(3,1,3);
    plot(t, rad2deg(delta_rad), 'g', 'LineWidth', 2);
    hold on;
    line([P.t_disturbance P.t_disturbance], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    title('VSM Internal Angle (Rotor Angle)');
    xlabel('Time (s)');
    ylabel('Angle (degrees)');
    grid on;
end

% ============================================================================
% THE DYNAMICS FUNCTION: Detailed VSM Swing Equation
% ============================================================================
function dxdt = vsm_detailed_dynamics(t, x, P)
    % Unpack state vector
    delta = x(1);
    w_vsm = x(2);
    p_meas = x(3);

    % Determine the mechanical power reference based on time
    if t < P.t_disturbance
        p_ref = P.P_initial;
    else
        p_ref = P.P_initial + P.P_step;
    end
    
    % Calculate the instantaneous electrical power based on the current angle
    p_elec = (P.V_terminal * P.V_grid / P.X_line) * sin(delta);
    
    % --- Low-Pass Filter for Power Measurement ---
    % This is a standard part of the VSM model
    dp_meas_dt = (1 / P.T_filter) * (p_elec - p_meas);

    % --- VSM Swing Equation (from papers) ---
    % Ta * d(w_vsm)/dt = p_ref - p_meas - Damping_Term
    % The damping term has two parts: kd * w_dev and kw * w_dev
    w_dev = w_vsm - 1.0; % Frequency deviation in pu
    damping_term = (P.kd + P.kw) * w_dev;
    
    dw_vsm_dt = (1 / P.Ta) * (p_ref - p_meas - damping_term);
    
    % The rate of change of the angle is the angular frequency
    % Need to convert w_vsm (pu) to rad/s deviation for integration
    ddelta_dt = w_dev * P.w_nom;

    dxdt = [ddelta_dt; dw_vsm_dt; dp_meas_dt];
end
