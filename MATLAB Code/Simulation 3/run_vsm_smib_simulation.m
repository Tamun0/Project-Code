function run_vsm_smib_simulation()
    % Main function to run a simulation of a Virtual Synchronous Machine (VSM)
    % connected to an infinite bus (Single Machine Infinite Bus - SMIB).
    % This model is based on the swing equation described in the provided papers
    % to highlight the VSM's inertial and damping characteristics.

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
    
    % --- VSM Control Parameters ---
    % These are the core VSM parameters from the swing equation
    P.J = 5.0;  % Virtual Inertia Constant (J or 2H in some notations)
    P.D = 0.8;  % Damping Coefficient
    
    % --- Simulation Event Times ---
    P.t_disturbance = 1.0; % Time of the power reference step
    t_span = [0 10.0];     % Simulate for 10 seconds to see oscillations
    
    disp('Running VSM (SMIB) INERTIAL RESPONSE simulation...');

    % ========================================================================
    % 2. SET UP AND RUN THE SIMULATION
    % ========================================================================
    % x = [delta, w_dev]
    % delta: the VSM's internal rotor angle (rad)
    % w_dev: the VSM's frequency deviation from nominal (rad/s)
    
    % Calculate initial steady-state angle for P_initial
    delta0 = asin(P.P_initial * P.X_line / (P.V_terminal * P.V_grid));
    x0 = [delta0; 0]; % Start at steady state: correct angle, zero freq deviation

    [t, x] = ode45(@(t,x) vsm_swing_equation(t, x, P), t_span, x0);

    % ========================================================================
    % 3. POST-PROCESS AND PLOT RESULTS
    % ========================================================================
    disp('Simulation finished. Processing results for plotting...');
    
    delta_rad = x(:,1);
    w_dev_rad_s = x(:,2);
    
    % Calculate other quantities for plotting
    frequency_hz = P.f_nom + w_dev_rad_s / (2*pi);
    electrical_power_pu = (P.V_terminal * P.V_grid / P.X_line) * sin(delta_rad);
    
    figure('Name', 'VSM Response to Power Step Command');
    
    % Plot Frequency
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
    line([P.t_disturbance P.t_disturbance], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    title('VSM Electrical Power Output');
    xlabel('Time (s)');
    ylabel('Power (pu)');
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
% THE DYNAMICS FUNCTION: VSM SWING EQUATION
% ============================================================================
function dxdt = vsm_swing_equation(t, x, P)
    % Unpack state vector
    delta = x(1);
    w_dev = x(2);

    % Determine the mechanical power reference based on time
    if t < P.t_disturbance
        P_mech = P.P_initial;
    else
        P_mech = P.P_initial + P.P_step;
    end
    
    % Calculate the electrical power based on the current angle
    % This is the standard power-angle curve equation
    P_elec = (P.V_terminal * P.V_grid / P.X_line) * sin(delta);
    
    % --- VSM Swing Equation ---
    % This is the core of the VSM, based on the models in the papers.
    % J * d(w_dev)/dt = P_mech - P_elec - D * w_dev
    
    dw_dev_dt = (1/P.J) * (P_mech - P_elec - P.D * w_dev);
    
    % The rate of change of the angle is the frequency deviation
    ddelta_dt = w_dev;

    dxdt = [ddelta_dt; dw_dev_dt];
end
