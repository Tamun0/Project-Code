function run_inverter_simulation_v7()
    % Main function to run a final, robust power system simulation.
    % This version uses a slower power ramp and a small virtual inductance
    % to ensure numerical stability for both strong and weak grid scenarios.

    clear; clc; close all;

    % --- Simulation Case Selection ---
    grid_type = 'strong'; % Can be 'strong' or 'weak'

    % ========================================================================
    % 1. DEFINE SYSTEM PARAMETERS (Struct 'P')
    % ========================================================================
    P.f_nom = 50; %Nominal Frequency
    P.w_nom = 2 * pi * P.f_nom;% Nominal Angular Frequency
    P.V_nom = 230; %Nominal Voltage
    P.V_peak = P.V_nom * sqrt(2); %Peak Voltage
    
    P.P_final = 10000; % Final active power setpoint (W)
    P.Q_final = 0;     % Final reactive power setpoint (VAR)
    
    P.mp = 0.05 / P.P_final; 
    P.mq = 0.05 / (abs(P.Q_final) + 1e-6);
    P.Cf = 1 / (2 * pi * 5);
    
    P.L_filt = 2e-3;
    P.R_filt = 0.1;
    
    % FIX 1: ADD VIRTUAL INDUCTANCE for numerical damping
    P.L_virtual = 10e-3; % A small virtual inductance (1 mH)
    
    if strcmp(grid_type, 'strong')
        P.R_grid = 0.1;
        P.L_grid = 0.1e-3;
        disp('Running STRONG GRID simulation...');
    else
        P.R_grid = 2.0;
        P.L_grid = 5.0e-3;
        disp('Running WEAK GRID simulation...');
    end
    
    % FIX 2: SLOW THE RAMP TIME FURTHER
    P.t_ramp_end = 0.3; % Ramp power over 0.3 seconds

    % ========================================================================
    % 2. SET UP THE SIMULATION
    % ========================================================================
    t_span = [0 0.5];

    x0 = [0; 0; 0; 0; 0];                     

    options = odeset('RelTol', 1e-4, 'AbsTol', 1e-5, 'MaxStep', 5e-5); 
    [t, x] = ode23t(@(t,x) system_dynamics_dq(t, x, P), t_span, x0, options);

    % ========================================================================
    % 3. POST-PROCESS AND PLOT RESULTS
    % ========================================================================
    disp('Simulation finished. Processing results for plotting...');
    
    id = x(:,1);
    iq = x(:,2);
    
    v_abc_terminal = zeros(length(t), 3);
    i_abc = zeros(length(t), 3);
    
    for k = 1:length(t)
       grid_angle = P.w_nom * t(k);
       T_inv_park = [cos(grid_angle), -sin(grid_angle);
                     cos(grid_angle - 2*pi/3), -sin(grid_angle - 2*pi/3);
                     cos(grid_angle + 2*pi/3), -sin(grid_angle + 2*pi/3)];
       
       v_d_terminal = P.V_peak + P.R_grid * id(k) - P.w_nom * P.L_grid * iq(k);
       v_q_terminal = 0 + P.R_grid * iq(k) + P.w_nom * P.L_grid * id(k);

       v_abc_terminal(k,:) = (T_inv_park * [v_d_terminal; v_q_terminal])';
       i_abc(k,:) = (T_inv_park * [id(k); iq(k)])';
    end
    
    figure('Name', ['Results for ' grid_type ' Grid']);
    
    subplot(2,1,1);
    plot(t, v_abc_terminal);
    title('Three-Phase Terminal Voltage at PCC (V)');
    xlabel('Time (s)'); ylabel('Voltage (V)');
    legend('Va', 'Vb', 'Vc'); grid on; ylim([-400 400]);

    subplot(2,1,2);
    plot(t, i_abc);
    title('Three-Phase Inverter Output Current (A)');
    xlabel('Time (s)'); ylabel('Current (A)');
    legend('Ia', 'Ib', 'Ic'); grid on;

end

% ============================================================================
% THE SYSTEM DYNAMICS FUNCTION (IN DQ FRAME)
% ============================================================================
function dxdt = system_dynamics_dq(t, x, P)
    % Unpack the state vector
    id = x(1);
    iq = x(2);
    theta_gfm = x(3); 
    P_filt = x(4);
    Q_filt = x(5);

    % --- Implement Soft-Start Power Ramp ---
    if t < P.t_ramp_end
        P_ref = (t / P.t_ramp_end) * P.P_final;
        Q_ref = (t / P.t_ramp_end) * P.Q_final;
    else
        P_ref = P.P_final;
        Q_ref = P.Q_final;
    end

    % --- GFM Droop Controller ---
    w_gfm = P.w_nom - (P_ref - P_filt) * P.mp;
    V_mag_ref = P.V_peak - (Q_ref - Q_filt) * P.mq;
    v_d_inv_ref = V_mag_ref;
    v_q_inv_ref = 0;
    
    % --- Grid Voltage and Coordinate Transformation ---
    v_d_grid = P.V_peak;
    v_q_grid = 0;
    v_d_inv = v_d_inv_ref * cos(theta_gfm) - v_q_inv_ref * sin(theta_gfm);
    v_q_inv = v_d_inv_ref * sin(theta_gfm) + v_q_inv_ref * cos(theta_gfm);

    % --- State-Space Equations for the R-L circuit ---
    % FIX: Add the virtual inductance to the total inductance
    L_total = P.L_filt + P.L_grid + P.L_virtual;
    R_total = P.R_filt + P.R_grid;
    
    di_d_dt = (1/L_total) * (-R_total * id + P.w_nom * L_total * iq + v_d_inv - v_d_grid);
    di_q_dt = (1/L_total) * (-R_total * iq - P.w_nom * L_total * id + v_q_inv - v_q_grid);
    
    % --- Controller State Dynamics ---
    dtheta_dt = w_gfm - P.w_nom;
    
    v_d_filt = v_d_inv - (P.R_filt * id - P.w_nom * P.L_filt * iq);
    v_q_filt = v_q_inv - (P.R_filt * iq + P.w_nom * P.L_filt * id);
    
    P_inst = v_d_filt * id + v_q_filt * iq;
    Q_inst = v_q_filt * id - v_d_filt * iq;
    dP_filt_dt = (1/P.Cf) * (P_inst - P_filt);
    dQ_filt_dt = (1/P.Cf) * (Q_inst - Q_filt);
    
    dxdt = [di_d_dt; di_q_dt; dtheta_dt; dP_filt_dt; dQ_filt_dt];
end
