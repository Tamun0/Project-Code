function run_grid_strength_test()
    % Main function to test the robustness of Droop, VSM, and VOC
    % controllers against varying grid strength (impedance).

    clear; clc; close all;

    % ========================================================================
    % 1. DEFINE BASE SYSTEM AND CONTROLLER PARAMETERS
    % ========================================================================
    P.f_nom = 50;
    P.w_nom = 2 * pi * P.f_nom;
    P.V_nom_rms = 230;
    P.V_peak = P.V_nom_rms * sqrt(2);
    
    P.P_rated = 10000; % Inverter rated power (W)
    P.Q_ref = 0;
    
    % --- Controller Parameters ---
    % Droop
    P.mp_droop = (0.05 * P.w_nom) / P.P_rated; % 5% frequency droop at rated power
    P.mq_droop = (0.05 * P.V_peak) / P.P_rated; % 5% voltage droop
    P.T_filter_droop = 0.05;
    % VSM
    P.J_vsm = 5.0;  % High virtual inertia (H=5s equivalent)
    P.D_vsm = 50.0; % Damping
    P.T_filter_vsm = 0.02;
    % AHO-VOC
    P.xi = 25; P.kv = P.V_peak; P.ki = 3 * P.V_peak / (P.P_rated + 1);
    P.C_voc = 0.1;

    % --- Simulation Scenario ---
    P.P_load = 0.8 * P.P_rated; % 80% load step
    t_span = [0 2];
    
    % ========================================================================
    % 2. DEFINE GRID STRENGTH SCENARIOS
    % ========================================================================
    % We will vary the grid impedance from strong (low Z) to weak (high Z)
    grid_impedances = logspace(-2, 0.5, 15); % From 0.01 to ~3 Ohms
    
    controllers = {'Droop', 'VSM', 'VOC'};
    results = struct();

    for c = 1:length(controllers)
        results.(controllers{c}).max_freq_dev_hz = [];
        results.(controllers{c}).min_volt_rms = [];
        
        for z = 1:length(grid_impedances)
            grid_Z = grid_impedances(z);
            
            % Differentiate between an inductive and a resistive grid
            is_resistive = (grid_Z > 1.0); % Simple flag for higher impedance resistive grid
            if is_resistive
                grid_R = grid_Z * 0.707;
                grid_X = grid_Z * 0.707; % X/R = 1
                P.phi = 0; % VOC setting for resistive grid
            else
                grid_R = grid_Z * 0.1;  
                grid_X = grid_Z * 0.99; % High X/R ratio
                P.phi = pi/2; % VOC setting for inductive grid
            end
            
            P.R_grid = grid_R;
            P.L_grid = grid_X / P.w_nom;
            
            control_mode = controllers{c};
            fprintf('Simulating %s with Grid Impedance Z = %.2f Ohm\n', control_mode, grid_Z);
            
            % --- Run Simulation ---
            if strcmp(control_mode, 'Droop'), x0 = [0;0]; else, x0 = [0;0;0]; end
            if strcmp(control_mode, 'VOC'), x0 = [P.V_peak; 0; 0; 0]; end

            options = odeset('RelTol', 1e-4, 'AbsTol', 1e-5);
            [t, x] = ode23t(@(t,x) gfm_dynamics(t, x, P, control_mode), t_span, x0, options);
            
            % --- Analyze and Store Results ---
            if strcmp(control_mode, 'Droop')
                w_dev = -P.mp_droop * (x(:,2) - P.P_load);
            else
                w_dev = x(:,2);
            end
            freq_hz = P.f_nom + w_dev / (2*pi);
            results.(controllers{c}).max_freq_dev_hz(end+1) = max(abs(freq_hz - P.f_nom));
            
            id = x(:,3); iq = x(:,4);
            v_d_term = P.V_peak - (P.R_grid * id - P.w_nom * P.L_grid * iq);
            v_q_term = 0 - (P.R_grid * iq + P.w_nom * P.L_grid * id);
            v_rms = sqrt(v_d_term.^2 + v_q_term.^2) / sqrt(2);
            results.(controllers{c}).min_volt_rms(end+1) = min(v_rms);
        end
    end

    % ========================================================================
    % 3. PLOT SUMMARY RESULTS
    % ========================================================================
    figure('Name', 'GFM Controller Robustness to Grid Strength', 'Position', [100 100 1200 600]);
    
    % --- Frequency Deviation Plot ---
    subplot(1,2,1);
    hold on;
    plot(grid_impedances, results.Droop.max_freq_dev_hz, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
    plot(grid_impedances, results.VSM.max_freq_dev_hz, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
    plot(grid_impedances, results.VOC.max_freq_dev_hz, 'g-^', 'LineWidth', 2, 'MarkerSize', 8);
    grid on; box on;
    title('Frequency Stability vs. Grid Impedance');
    xlabel('Grid Impedance (Ohms)');
    ylabel('Max Frequency Deviation (Hz)');
    legend('Droop', 'VSM', 'VOC', 'Location', 'northwest');
    set(gca, 'FontSize', 12, 'XScale', 'log');

    % --- Voltage Sag Plot ---
    subplot(1,2,2);
    hold on;
    plot(grid_impedances, results.Droop.min_volt_rms, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
    plot(grid_impedances, results.VSM.min_volt_rms, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
    plot(grid_impedances, results.VOC.min_volt_rms, 'g-^', 'LineWidth', 2, 'MarkerSize', 8);
    grid on; box on;
    title('Voltage Sag vs. Grid Impedance');
    xlabel('Grid Impedance (Ohms)');
    ylabel('Minimum Terminal Voltage (V_{RMS})');
    legend('Droop', 'VSM', 'VOC', 'Location', 'southwest');
    set(gca, 'FontSize', 12, 'XScale', 'log');

end

% ============================================================================
% THE UNIFIED DYNAMICS FUNCTION
% ============================================================================
function dxdt = gfm_dynamics(t, x, P, mode)
    P_ref = P.P_load; % Simple step load
    
    if strcmp(mode, 'VOC')
        v_alpha = x(1); v_beta = x(2); id = x(3); iq = x(4);
        i_alpha = id*cos(P.w_nom*t) - iq*sin(P.w_nom*t);
        i_beta = id*sin(P.w_nom*t) + iq*cos(P.w_nom*t);
        i_alpha_star = (2/(3*(v_alpha^2+v_beta^2))) * (v_alpha*P_ref + v_beta*P.Q_ref);
        i_beta_star = (2/(3*(v_alpha^2+v_beta^2))) * (v_beta*P_ref - v_alpha*P.Q_ref);
        i_err = [i_alpha-i_alpha_star; i_beta-i_beta_star];
        R_phi = [cos(P.phi) -sin(P.phi); sin(P.phi) cos(P.phi)];
        u = P.ki * R_phi * i_err;
        V_sq = v_alpha^2 + v_beta^2; V_nom_sq = P.kv^2;
        dv_alpha_dt = (P.xi/P.kv^2)*(2*V_nom_sq-V_sq)*v_alpha - P.w_nom*v_beta - (P.kv/P.C_voc)*u(1);
        dv_beta_dt = (P.xi/P.kv^2)*(2*V_nom_sq-V_sq)*v_beta + P.w_nom*v_alpha - (P.kv/P.C_voc)*u(2);
        
        v_d_inv = v_alpha*cos(P.w_nom*t) + v_beta*sin(P.w_nom*t);
        v_q_inv = -v_alpha*sin(P.w_nom*t) + v_beta*cos(P.w_nom*t);
        
        L_total = P.L_grid; R_total = P.R_grid;
        did_dt = (1/L_total)*(-R_total*id + P.w_nom*L_total*iq + v_d_inv - P.V_peak);
        diq_dt = (1/L_total)*(-R_total*iq - P.w_nom*L_total*id + v_q_inv - 0);
        dxdt = [dv_alpha_dt; dv_beta_dt; did_dt; diq_dt];
        return;
    end
    
    delta = x(1);
    if ~strcmp(mode, 'Droop'), w_dev = x(2); end
    p_filt = x(end);
    
    V_grid = P.V_peak;
    Z_line = sqrt(P.R_grid^2 + (P.w_nom * P.L_grid)^2);
    P_elec = (V_grid^2 / Z_line) * sin(delta);
    
    if strcmp(mode, 'Droop')
        dp_filt_dt = (1/P.T_filter_droop) * (P_elec - p_filt);
        ddelta_dt = -P.mp_droop * (p_filt - P_ref);
        dxdt = [ddelta_dt; dp_filt_dt];
    else % VSM
        dp_filt_dt = (1/P.T_filter_vsm) * (P_elec - p_filt);
        dw_dev_dt = (1/P.J_vsm) * (P_ref - p_filt - P.D_vsm * w_dev);
        ddelta_dt = w_dev;
        dxdt = [ddelta_dt; dw_dev_dt; dp_filt_dt];
    end
end