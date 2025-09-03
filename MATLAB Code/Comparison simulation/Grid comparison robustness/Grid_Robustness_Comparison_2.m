% GFM Controller Comparison: Droop, VSM, and VOC
% This script simulates and compares the transient response of three
% grid-forming (GFM) controllers to a load step in both inductive and
% resistive grid environments. The dynamics for each controller are based
% on the principles discussed in the provided research papers.

% Clear workspace and close all figures
clear; clc; close all;

%% System and Simulation Parameters
P.V_nom = 230;              % Nominal phase voltage (V)
P.f_nom = 50;               % Nominal frequency (Hz)
P.w_nom = 2*pi*P.f_nom;     % Nominal angular frequency (rad/s)
P.S_rated = 10e3;           % Inverter rated power (VA)
P.P_rated = P.S_rated;      % For simplicity, P_rated = S_rated
P.Q_rated = P.S_rated;

% Simulation time
t_start = 0;
t_end = 5;
t_step = 1e-4;
t = t_start:t_step:t_end;
n_steps = length(t);

% Load step event
t_load_step = 1; % Time of load step (s)
P_load_initial = 0.2 * P.P_rated; % Initial load (20% of rated)
P_load_final = 0.8 * P.P_rated;   % Final load (80% of rated)
Q_load = 0; % Assuming purely resistive load for simplicity

% Create load profile
P_load = P_load_initial * ones(1, n_steps);
P_load(t >= t_load_step) = P_load_final;

%% Grid Models
% The grid is modeled as a series impedance with different X/R ratios.
% This distinguishes between Inductive and Resistive characteristics.

% 1. Inductive Grid (High X/R ratio)
Grid(1).R = 0.1; % Low resistance
Grid(1).L = 0.01; % High inductance
Grid(1).X = P.w_nom * Grid(1).L;
Grid(1).name = 'Inductive Grid (X/R >> 1)';

% 2. Resistive Grid (Low X/R ratio)
Grid(2).R = 1.0; % High resistance
Grid(2).L = 0.001; % Low inductance
Grid(2).X = P.w_nom * Grid(2).L;
Grid(2).name = 'Resistive Grid (X/R << 1)';

%% Controller Models and Parameters
Controllers = {};

% --- 1. Droop Controller ---
% Dynamics based on "Comparative Study of Droop Control Methods..."
C_droop.kp_f = (0.5 * 2 * pi) / P.P_rated; % P-f droop gain
C_droop.kp_V = (0.05 * P.V_nom) / P.Q_rated; % Q-V droop gain
C_droop.Tc = 0.1; % Low-pass filter time constant

% Conventional Droop for Inductive Grids (P-f, Q-V)
Controllers{1}.name = 'Droop (Conventional)';
Controllers{1}.type = 'Droop';
Controllers{1}.params = C_droop;
Controllers{1}.params.mode = 'Inductive';

% Inverse Droop for Resistive Grids (P-V, Q-f)
Controllers{2}.name = 'Droop (Inverse)';
Controllers{2}.type = 'Droop';
Controllers{2}.params = C_droop;
Controllers{2}.params.mode = 'Resistive';

% --- 2. Virtual Synchronous Machine (VSM) ---
% Dynamics based on "Design and Analysis of Virtual Synchronous Machines..."
C_vsm.J = 2 * 0.5 / P.w_nom; % Virtual Inertia (H=0.5s)
C_vsm.Dp = P.P_rated / (0.5 * 2 * pi); % Damping factor
C_vsm.Kq = 100; % Reactive power integral gain
C_vsm.Dq = (0.05 * P.V_nom) / P.Q_rated; % Voltage droop gain

Controllers{3}.name = 'VSM';
Controllers{3}.type = 'VSM';
Controllers{3}.params = C_vsm;

% --- 3. Virtual Oscillator Control (VOC) ---
% Dynamics based on "A Grid-compatible Virtual Oscillator Controller..."
C_voc.xi = 15; % Speed constant
C_voc.kv = P.V_nom;
C_voc.ki = 3 * P.V_nom / P.S_rated;
C_voc.C = 0.2679; % Virtual capacitance
C_voc.L = 26.268e-6; % Virtual inductance

% VOC for Inductive Grids (phi = pi/2)
Controllers{4}.name = 'VOC (Inductive)';
Controllers{4}.type = 'VOC';
Controllers{4}.params = C_voc;
Controllers{4}.params.phi = pi/2;

% VOC for Resistive Grids (phi = 0)
Controllers{5}.name = 'VOC (Resistive)';
Controllers{5}.type = 'VOC';
Controllers{5}.params = C_voc;
Controllers{5}.params.phi = 0;


%% Simulation Loop
results = struct();

for g = 1:length(Grid)
    for c = 1:length(Controllers)
        
        fprintf('Simulating: %s on %s\n', Controllers{c}.name, Grid(g).name);
        
        % Initialize states
        w = P.w_nom;
        theta = 0;
        V = P.V_nom;
        
        % Initialize controller-specific states
        if strcmp(Controllers{c}.type, 'Droop')
            P_filt = P_load_initial;
            Q_filt = Q_load;
        elseif strcmp(Controllers{c}.type, 'VSM')
            psi_v = V / w; % Virtual flux
        elseif strcmp(Controllers{c}.type, 'VOC')
            v_alpha = V; v_beta = 0; % Oscillator states
        end
        
        % History arrays for plotting
        w_hist = zeros(1, n_steps);
        V_hist = zeros(1, n_steps);
        
        % Time-domain simulation
        for i = 1:n_steps
            Z_grid = sqrt(Grid(g).R^2 + Grid(g).X^2);
            P_out = (V * P.V_nom / Z_grid) * sin(theta) * (Grid(g).X / Z_grid) + (V * P.V_nom / Z_grid) * cos(theta) * (Grid(g).R / Z_grid) - (V^2 / Z_grid^2) * Grid(g).R;
            Q_out = (V * P.V_nom / Z_grid) * cos(theta) * (Grid(g).X / Z_grid) - (V * P.V_nom / Z_grid) * sin(theta) * (Grid(g).R / Z_grid) - (V^2 / Z_grid^2) * Grid(g).X;
            
            switch Controllers{c}.type
                case 'Droop'
                    p = Controllers{c}.params;
                    dP_filt = (1/p.Tc) * (P_out - P_filt);
                    dQ_filt = (1/p.Tc) * (Q_out - Q_filt);
                    P_filt = P_filt + dP_filt * t_step;
                    Q_filt = Q_filt + dQ_filt * t_step;
                    if strcmp(p.mode, 'Inductive')
                        w = P.w_nom - p.kp_f * P_filt; V = P.V_nom - p.kp_V * Q_filt;
                    else
                        w = P.w_nom - p.kp_f * Q_filt; V = P.V_nom - p.kp_V * P_filt;
                    end
                case 'VSM'
                    p = Controllers{c}.params;
                    Te = P_out / w; Tm = P_load(i) / P.w_nom;
                    dw = (1/p.J) * (Tm - Te - p.Dp * (w - P.w_nom));
                    d_psi_v = p.Kq * (p.Dq*(P.V_nom - V) - Q_out);
                    w = w + dw * t_step;
                    psi_v = psi_v + d_psi_v * t_step;
                    V = psi_v * w;
                case 'VOC'
                    p = Controllers{c}.params;
                    P_star = P_load(i); Q_star = 0;
                    i_alpha_star = (2 / (3 * (v_alpha^2 + v_beta^2))) * (v_alpha * P_star + v_beta * Q_star);
                    i_beta_star = (2 / (3 * (v_alpha^2 + v_beta^2))) * (v_beta * P_star - v_alpha * Q_star);
                    i_alpha = (2 / (3 * (v_alpha^2 + v_beta^2))) * (v_alpha * P_out + v_beta * Q_out);
                    i_beta = (2 / (3 * (v_alpha^2 + v_beta^2))) * (v_beta * P_out - v_alpha * Q_out);
                    R_phi = [cos(p.phi) -sin(p.phi); sin(p.phi) cos(p.phi)];
                    i_err = [i_alpha - i_alpha_star; i_beta - i_beta_star];
                    u = p.ki * R_phi * i_err;
                    V_sq = v_alpha^2 + v_beta^2; V_nom_sq = (p.kv)^2;
                    dv_alpha = (p.xi / p.kv^2) * (2*V_nom_sq - V_sq) * v_alpha - P.w_nom * v_beta - (p.kv / p.C) * u(1);
                    dv_beta  = (p.xi / p.kv^2) * (2*V_nom_sq - V_sq) * v_beta  + P.w_nom * v_alpha - (p.kv / p.C) * u(2);
                    v_alpha = v_alpha + dv_alpha * t_step;
                    v_beta = v_beta + dv_beta * t_step;
                    V = sqrt(v_alpha^2 + v_beta^2);
                    w = P.w_nom - (p.kv*p.ki / (3*p.C*V^2)) * (sin(p.phi)*(P_out-P_star) - cos(p.phi)*(Q_out-Q_star));
            end
            theta = theta + w * t_step;
            w_hist(i) = w; V_hist(i) = V;
        end
        results(g,c).w = w_hist; results(g,c).V = V_hist;
        results(g,c).controller = Controllers{c}; results(g,c).grid = Grid(g);
    end
end

%% Plotting Results - Replicated and Enhanced Structure
colors = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.4660, 0.6740, 0.1880];
line_styles = {'-', '--', ':'};
font_size = 12;

% --- Plot for Inductive Grid ---
fig_inductive = figure('Name', 'Controller Performance in Inductive Grid', 'Position', [100, 100, 800, 900]);
sgtitle('Controller Response to a 60% Load Step in an Inductive Grid', 'FontSize', font_size + 4, 'FontWeight', 'bold');

% Subplot 1: Frequency Response
subplot(2, 1, 1);
hold on;
plot_handles_ind = [];
legend_entries_ind = {};
controllers_to_plot_ind = {'Droop (Conventional)', 'VSM', 'VOC (Inductive)'};
for c = 1:length(Controllers)
    if ismember(Controllers{c}.name, controllers_to_plot_ind)
        idx = find(strcmp(controllers_to_plot_ind, Controllers{c}.name));
        h = plot(t, results(1,c).w / (2*pi), 'LineWidth', 2.5, 'LineStyle', line_styles{idx}, 'Color', colors(idx,:));
        plot_handles_ind = [plot_handles_ind, h];
        legend_entries_ind{end+1} = Controllers{c}.name;
    end
end
grid on; box on;
title('Frequency Response', 'FontSize', font_size + 2);
xlabel('Time (s)', 'FontSize', font_size);
ylabel('Frequency (Hz)', 'FontSize', font_size);
ylim([P.f_nom - 1, P.f_nom + 0.5]);
xline(t_load_step, 'k--', 'Load Step', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
legend(plot_handles_ind, legend_entries_ind, 'Location', 'best', 'FontSize', font_size - 1);
set(gca, 'FontSize', font_size);

% Subplot 2: Voltage Response
subplot(2, 1, 2);
hold on;
for c = 1:length(Controllers)
    if ismember(Controllers{c}.name, controllers_to_plot_ind)
        idx = find(strcmp(controllers_to_plot_ind, Controllers{c}.name));
        plot(t, results(1,c).V, 'LineWidth', 2.5, 'LineStyle', line_styles{idx}, 'Color', colors(idx,:));
    end
end
grid on; box on;
title('Voltage Response', 'FontSize', font_size + 2);
xlabel('Time (s)', 'FontSize', font_size);
ylabel('Voltage (V)', 'FontSize', font_size);
ylim([P.V_nom * 0.9, P.V_nom * 1.05]);
xline(t_load_step, 'k--', 'Load Step', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
legend(legend_entries_ind, 'Location', 'best', 'FontSize', font_size - 1);
set(gca, 'FontSize', font_size);

% --- Plot for Resistive Grid ---
fig_resistive = figure('Name', 'Controller Performance in Resistive Grid', 'Position', [950, 100, 800, 900]);
sgtitle('Controller Response to a 60% Load Step in a Resistive Grid', 'FontSize', font_size + 4, 'FontWeight', 'bold');

% Subplot 1: Frequency Response
subplot(2, 1, 1);
hold on;
plot_handles_res = [];
legend_entries_res = {};
controllers_to_plot_res = {'Droop (Inverse)', 'VOC (Resistive)'};
for c = 1:length(Controllers)
    if ismember(Controllers{c}.name, controllers_to_plot_res)
        idx = find(strcmp(controllers_to_plot_res, Controllers{c}.name));
        h = plot(t, results(2,c).w / (2*pi), 'LineWidth', 2.5, 'LineStyle', line_styles{idx}, 'Color', colors(idx,:));
        plot_handles_res = [plot_handles_res, h];
        legend_entries_res{end+1} = Controllers{c}.name;
    end
end
grid on; box on;
title('Frequency Response', 'FontSize', font_size + 2);
xlabel('Time (s)', 'FontSize', font_size);
ylabel('Frequency (Hz)', 'FontSize', font_size);
ylim([P.f_nom - 1, P.f_nom + 0.5]);
xline(t_load_step, 'k--', 'Load Step', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
legend(plot_handles_res, legend_entries_res, 'Location', 'best', 'FontSize', font_size - 1);
set(gca, 'FontSize', font_size);

% Subplot 2: Voltage Response
subplot(2, 1, 2);
hold on;
for c = 1:length(Controllers)
    if ismember(Controllers{c}.name, controllers_to_plot_res)
        idx = find(strcmp(controllers_to_plot_res, Controllers{c}.name));
        plot(t, results(2,c).V, 'LineWidth', 2.5, 'LineStyle', line_styles{idx}, 'Color', colors(idx,:));
    end
end
grid on; box on;
title('Voltage Response', 'FontSize', font_size + 2);
xlabel('Time (s)', 'FontSize', font_size);
ylabel('Voltage (V)', 'FontSize', font_size);
ylim([P.V_nom * 0.9, P.V_nom * 1.05]);
xline(t_load_step, 'k--', 'Load Step', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
legend(legend_entries_res, 'Location', 'best', 'FontSize', font_size - 1);
set(gca, 'FontSize', font_size);

disp('Simulation complete. Plots generated.');
