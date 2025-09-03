function run_voc_synchronization_sim()
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
    
    % --- dVOC Controller Parameters (from paper) ---
    P.eta = 20;  % Synchronization gain
    P.alpha = 5; % Voltage magnitude control gain
    P.kappa = 0; % Angle for network impedance (0 for resistive grid)
    
    % --- Power Setpoints (per unit, assuming S_base = 10kVA) ---
    P.p_ref1 = 0.5; % Inverter 1 active power
    P.q_ref1 = 0;   % Inverter 1 reactive power
    P.p_ref2 = 0.5; % Inverter 2 active power
    P.q_ref2 = 0;   % Inverter 2 reactive power
    
    % --- Circuit Parameters ---
    P.R_load = (P.V_peak^2 / 2) / (10000 * (P.p_ref1 + P.p_ref2)); % Resistive load
    
    % --- Simulation Setup ---
    t_span = [0 1.5];
    
    disp('Running dVOC SELF-SYNCHRONIZATION simulation...');

    % ========================================================================
    % 2. SET UP AND RUN THE SIMULATION
    % ========================================================================
    % x = [v1_alpha, v1_beta, v2_alpha, v2_beta]
    % The state is the internal voltage vector of each oscillator.
    
    % Initial conditions: Start them out of sync.
    % Inverter 1 starts at angle 0.
    v1_0 = [P.V_peak; 0]; 
    % Inverter 2 starts at a significant angle offset (e.g., 60 degrees).
    angle_offset = pi/3;
    v2_0 = [P.V_peak * cos(angle_offset); P.V_peak * sin(angle_offset)];
    
    x0 = [v1_0; v2_0];

    [t, x] = ode23t(@(t,x) voc_dynamics(t, x, P), t_span, x0);

    % ========================================================================
    % 3. POST-PROCESS AND PLOT RESULTS
    % ========================================================================
    disp('Simulation finished. Processing results for plotting...');
    
    v1_alpha = x(:,1); v1_beta = x(:,2);
    v2_alpha = x(:,3); v2_beta = x(:,4);
    
    % Calculate phase angles and power for plotting
    angle1 = atan2(v1_beta, v1_alpha);
    angle2 = atan2(v2_beta, v2_alpha);
    angle_diff_deg = rad2deg(unwrap(angle1 - angle2)); % unwrap handles 2pi jumps
    
    i_load_alpha = (v1_alpha + v2_alpha) / P.R_load;
    i_load_beta  = (v1_beta + v2_beta) / P.R_load;
    
    % Assume inverters share current equally for power calculation
    P1_out = (v1_alpha .* i_load_alpha/2 + v1_beta .* i_load_beta/2) / 1000; % in kW
    P2_out = (v2_alpha .* i_load_alpha/2 + v2_beta .* i_load_beta/2) / 1000; % in kW

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
    plot(t, P1_out, 'm', 'LineWidth', 2);
    hold on;
    plot(t, P2_out, 'c', 'LineWidth', 2);
    title('Inverter Active Power Outputs');
    xlabel('Time (s)');
    ylabel('Power (kW)');
    legend('Inverter 1', 'Inverter 2');
    grid on;
end

% ============================================================================
% THE DYNAMICS FUNCTION: dVOC Equation
% ============================================================================
function dxdt = voc_dynamics(t, x, P)
    % Unpack state vector
    v1 = x(1:2); % Voltage vector for Inverter 1
    v2 = x(3:4); % Voltage vector for Inverter 2

    % --- Calculate total load current based on combined voltage ---
    v_total = v1 + v2;
    i_load = v_total / P.R_load;
    
    % For this behavioral model, assume each inverter provides half the current
    i1 = i_load / 2;
    i2 = i_load / 2;

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
    % K matrix contains power setpoints
    K = (1/P.V_peak^2) * [p_ref, -q_ref; q_ref, p_ref];
    R_kappa = [cos(P.kappa), -sin(P.kappa); sin(P.kappa), cos(P.kappa)];
    term2 = P.eta * (K*v - R_kappa*i);
    
    % --- Magnitude Control Term ---
    Phi = (P.V_peak^2 - (v(1)^2 + v(2)^2)) / P.V_peak^2;
    term3 = P.eta * P.alpha * Phi * v;
    
    dvdt = term1 + term2 + term3;
end