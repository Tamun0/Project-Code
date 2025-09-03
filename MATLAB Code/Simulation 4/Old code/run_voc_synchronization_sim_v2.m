function run_voc_synchronization_sim_v2()
    % Main function to run a simulation of two dVOC inverters autonomously
    % synchronizing in an islanded microgrid. This model is based on the
    % dVOC equations from the paper by Dörfler & Groß.

    clear; clc; close all;

    % ========================================================================
    % 1. DEFINE SYSTEM PARAMETERS (Struct 'P')
    % ========================================================================
    P.f_nom = 50;
    P.w_nom = 2 * pi * P.f_nom;
    P.V_peak = 325; % Peak phase voltage (approx. 230V RMS)
    P.S_base = 10000; % Base power for per-unit calculations
    
    % --- dVOC Controller Parameters (from paper) ---
    P.eta = 50;  % Synchronization gain (increased for faster sync)
    P.alpha = 5; % Voltage magnitude control gain
    P.kappa = -pi/2; % Angle for network impedance (pi/2 for inductive, 0 for resistive)
    
    % --- Power Setpoints (per unit) ---
    P.p_ref1 = 0.5; % Inverter 1 active power (5 kW)
    P.q_ref1 = 0;   % Inverter 1 reactive power
    P.p_ref2 = 0.3; % Inverter 2 active power (3 kW)
    P.q_ref2 = 0;   % Inverter 2 reactive power
    
    % --- Circuit Parameters ---
    total_power = P.S_base * (P.p_ref1 + P.p_ref2);
    P.R_load = (P.V_peak^2 / 2) / total_power; % Resistive load
    P.X_line = 0.1; % Small line reactance between inverters
    
    % --- Simulation Setup ---
    t_span = [0 1.5];
    
    disp('Running dVOC SELF-SYNCHRONIZATION simulation...');

    % ========================================================================
    % 2. SET UP AND RUN THE SIMULATION
    % ========================================================================
    % x = [v1_alpha, v1_beta, v2_alpha, v2_beta]
    % The state is the internal voltage vector of each oscillator.
    
    % Initial conditions: Start them out of sync.
    v1_0 = [P.V_peak; 0]; % Inverter 1 starts at angle 0.
    angle_offset = pi/2;  % Start Inverter 2 at 90 degrees.
    v2_0 = [P.V_peak * cos(angle_offset); P.V_peak * sin(angle_offset)];
    
    x0 = [v1_0; v2_0];

    [t, x] = ode23t(@(t,x) voc_dynamics(t, x, P), t_span, x0);

    % ========================================================================
    % 3. POST-PROCESS AND PLOT RESULTS
    % ========================================================================
    disp('Simulation finished. Processing results for plotting...');
    
    % FIX: Correctly extract state variables from the output matrix 'x'
    v1_alpha = x(:,1); v1_beta = x(:,2);
    v2_alpha = x(:,3); v2_beta = x(:,4);
    
    % Calculate phase angles and power for plotting
    angle1 = atan2(v1_beta, v1_alpha);
    angle2 = atan2(v2_beta, v2_alpha);
    angle_diff_deg = rad2deg(unwrap(angle1 - angle2)); 
    
    % Re-calculate currents for plotting using complex numbers for clarity
    v1_c = v1_alpha + 1i * v1_beta;
    v2_c = v2_alpha + 1i * v2_beta;
    i12_c = (v1_c - v2_c) / (1i * P.X_line);
    i_load1_c = v1_c / (2 * P.R_load);
    i_load2_c = v2_c / (2 * P.R_load);
    i1_c = i12_c + i_load1_c;
    i2_c = -i12_c + i_load2_c;
    
    % Calculate power output for each inverter
    P1_out_kw = (real(v1_c) .* real(i1_c) + imag(v1_c) .* imag(i1_c)) / 1000;
    P2_out_kw = (real(v2_c) .* real(i2_c) + imag(v2_c) .* imag(i2_c)) / 1000;

    figure('Name', 'dVOC Self-Synchronization');
    
    % Plot Angle Difference
    subplot(2,1,1);
    plot(t, angle_diff_deg, 'b', 'LineWidth', 2);
    title('Phase Angle Difference Between Inverters');
    xlabel('Time (s)');
    ylabel('Angle Difference (degrees)');
    grid on;

    % Plot Power Outputs
    subplot(2,1,2);
    plot(t, P1_out_kw, 'm', 'LineWidth', 2);
    hold on;
    plot(t, P2_out_kw, 'c', 'LineWidth', 2);
    title('Inverter Active Power Outputs');
    xlabel('Time (s)');
    ylabel('Power (kW)');
    legend('Inverter 1 (5kW ref)', 'Inverter 2 (3kW ref)');
    grid on;
end

% ============================================================================
% THE DYNAMICS FUNCTION: dVOC Equation
% ============================================================================
function dxdt = voc_dynamics(t, x, P)
    % Unpack state vector
    v1 = x(1:2); % Voltage vector for Inverter 1
    v2 = x(3:4); % Voltage vector for Inverter 2

    % --- Calculate current based on voltage difference across line impedance ---
    % Convert to complex numbers for easier calculation
    v1_c = v1(1) + 1i * v1(2);
    v2_c = v2(1) + 1i * v2(2);
    
    % Current flowing from inverter 1 to inverter 2
    i12_c = (v1_c - v2_c) / (1i * P.X_line);
    
    % Assume load is split between the two nodes
    i_load1_c = v1_c / (2 * P.R_load);
    i_load2_c = v2_c / (2 * P.R_load);
    
    % Total current from each inverter
    i1_c = i12_c + i_load1_c;
    i2_c = -i12_c + i_load2_c;
    
    % Convert back to vectors
    i1 = [real(i1_c); imag(i1_c)];
    i2 = [real(i2_c); imag(i2_c)];

    % --- dVOC Equation for each inverter (from paper Eq. 24) ---
    dv1_dt = calc_dvdt(v1, i1, P.p_ref1, P.q_ref1, P);
    dv2_dt = calc_dvdt(v2, i2, P.p_ref2, P.q_ref2, P);

    dxdt = [dv1_dt; dv2_dt];
end

function dvdt = calc_dvdt(v, i, p_ref, q_ref, P)
    % This helper function implements the dVOC equation for a single inverter
    
    J = [0 -1; 1 0]; % Rotation matrix
    
    % --- Harmonic Oscillator Term ---
    term1 = P.w_nom * J * v;
    
    % --- Synchronization Term ---
    % K matrix contains power setpoints (in per-unit)
    K = (1/P.V_peak^2) * (P.S_base * [p_ref, -q_ref; q_ref, p_ref]);
    R_kappa = [cos(P.kappa), -sin(P.kappa); sin(P.kappa), cos(P.kappa)];
    term2 = P.eta * (K*v - R_kappa*i);
    
    % --- Magnitude Control Term ---
    Phi = (P.V_peak^2 - (v(1)^2 + v(2)^2)) / P.V_peak^2;
    term3 = P.eta * P.alpha * Phi * v;
    
    dvdt = term1 + term2 + term3;
end
