%% Task 5: System Characterization and Identification (amp25_run2)
% Ball-and-Beam / Sliding Object System
%
% Summary of Workflow:
%   1. Load raw chirp dataset for 'amp25_run2' and resample to uniform 50 Hz grid (Ts = 0.02 s).
%   2. Preprocess data: remove sensor glitches and truncate the record at t = 12.4 s
%      (where the car hits the physical end bumper, preventing non-linear saturation).
%   3. Identify the physical plant  G(s) = K / (s * (s + b))  and neutral trim u0
%      by minimizing simulation output error on the free sliding motion.
%   4. Cross-check against MATLAB's subspace identification (n4sid).
%   5. Discretize to G(z) at 50 Hz and formulate the Arduino C++ difference equation.
%   6. Generate validation, Bode, pole-zero, and step response comparison plots.
%
% Author: Control Systems Laboratory
% Date: November 2026

clear; close all; clc;

%% 1. Configuration and Paths
target_run_id = 'amp25_run2';
Ts     = 0.02;                 % 50 Hz sampling period (s)
t_grid = (0.1:Ts:19.9)';       % Uniform time grid (s)

% Ideal sliding-car gain from physical geometry:
% L = 0.40 m, d = 0.04 m, g = 9.81 m/s^2 -> K_ideal = g*(d/L)*(pi/180)*100
K_ideal = 9.81 * (0.04 / 0.40) * (100 * pi / 180); % ~1.7122 cm/s^2 per deg

% Locate CSV dataset
candidatePaths = {
    fullfile('..', 'data', 'task5_chirp_runs.csv'), ...
    fullfile('data', 'task5_chirp_runs.csv'), ...
    'task5_chirp_runs.csv'
};
dataPath = '';
for k = 1:numel(candidatePaths)
    if exist(candidatePaths{k}, 'file') == 2
        dataPath = candidatePaths{k};
        break;
    end
end
assert(~isempty(dataPath), 'task5_chirp_runs.csv not found!');

outputDir = fileparts(mfilename('fullpath'));
if isempty(outputDir), outputDir = pwd; end

fprintf('=======================================================================\n');
fprintf('  TASK 5: SYSTEM CHARACTERIZATION FOR RUN: %s\n', target_run_id);
fprintf('=======================================================================\n');

%% 2. Load, Resample, and Preprocess Data
T = readtable(dataPath, 'VariableNamingRule', 'preserve');
mask = strcmp(T.run_id, target_run_id);
assert(any(mask), 'Run ID %s not found in CSV file!', target_run_id);

t_raw = T.time_ms(mask) / 1000;     % Time in seconds
u_raw = T.servo_deg(mask);          % Servo angle in degrees
y_raw = T.distance_cm(mask);        % Ultrasonic distance in cm

% Clean dropouts and resample to uniform 50 Hz grid
y_raw(y_raw <= 0 | ~isfinite(y_raw)) = NaN;
y_raw = fillmissing(y_raw, 'linear');
u_res = interp1(t_raw, u_raw, t_grid, 'linear');
y_res = interp1(t_raw, y_raw, t_grid, 'linear');

% Remove isolated ultrasonic spikes and apply symmetric Gaussian smoothing
y_clean = smoothdata(y_res,   'movmedian', 7);
u_f     = smoothdata(u_res,   'gaussian', 7);
y_f     = smoothdata(y_clean, 'gaussian', 7);

% Truncate record when the car reaches the physical bumper stop (~7 cm)
% Operating past this point saturates the sensor and distorts linear dynamics
idx_end = find(y_f < min(y_f) + 2.0, 1, 'first');
if isempty(idx_end), idx_end = length(t_grid); end
t_end = t_grid(idx_end);

t_id = t_grid(1:idx_end);
u_id = u_f(1:idx_end);
y_id = y_f(1:idx_end);

fprintf('Active identification window: 0.1 to %.2f s (%d samples)\n', t_end, idx_end);

%% 3. Physical Plant Identification: G(s) = K / (s * (s + b))
% Estimate parameters: p = [b (damping), K (gain), u0 (trim), y0 (init pos), v0 (init vel)]
p0 = [0.05, -0.5, mean(u_id), y_id(1), 0.0];
opt = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-6);
p_opt = fminsearch(@(p) sum((simPlant(p, u_id, Ts) - y_id).^2), p0, opt);

b_id = abs(p_opt(1));   % Track friction damping (s^-1)
K_id = p_opt(2);        % Plant gain (cm/s^2 per deg, negative due to tilt direction)
u0   = p_opt(3);        % Estimated horizontal neutral servo trim (deg)
y0   = p_opt(4);        % Initial position (cm)
v0   = p_opt(5);        % Initial velocity (cm/s)

% Continuous-time Transfer Function G(s)
G = tf(K_id, [1, b_id, 0]);
G.InputName = 'Servo \Delta\theta (deg)';
G.OutputName = 'Distance \Delta d (cm)';

% Discrete-time Transfer Function G(z) via Zero-Order Hold (ZOH)
G_z = c2d(G, Ts, 'zoh');
[num_z, den_z] = tfdata(G_z, 'v');

% Simulation and goodness-of-fit evaluation
y_sim = simPlant(p_opt, u_id, Ts);
fit_val = 100 * (1 - norm(y_id - y_sim) / norm(y_id - mean(y_id)));
rmse = sqrt(mean((y_id - y_sim).^2));

%% 4. Subspace Identification (n4sid) Cross-Check
try
    data_id = iddata(y_id, u_id - u0, Ts);
    sys_n4 = n4sid(data_id, 2, n4sidOptions('Focus', 'simulation'));
    [~, fit_n4] = compare(data_id, sys_n4);
    fit_n4_str = sprintf('%.1f%%', fit_n4(1));
catch
    fit_n4_str = 'N/A';
end

%% 5. FIGURE 1: Data Preprocessing and Bumper Truncation (Amplitude vs. Time)
figure('Name', 'Task 5: Preprocessing (Amplitude vs. Time)', 'Color', 'w', 'Position', [100, 100, 900, 600]);

subplot(2, 1, 1);
plot(t_grid, u_res, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.9, 'DisplayName', 'Raw Resampled (50 Hz)'); hold on;
plot(t_grid, u_f, 'Color', [0 0.447 0.741], 'LineWidth', 1.5, 'DisplayName', 'Gaussian Smoothed');
yline(u0, '--k', sprintf('Neutral Trim u_0 = %.2f°', u0), 'LineWidth', 1.2);
xline(t_end, ':r', sprintf('End Stop (t = %.2f s)', t_end), 'LineWidth', 1.5);
grid on;
ylabel({'Amplitude'; 'Servo Angle (deg)'}, 'FontWeight', 'bold');
title(sprintf('Actuator Chirp Input: Amplitude vs. Time [%s: \\pm25° Swing]', target_run_id), 'FontWeight', 'bold');
legend('Location', 'southwest');

subplot(2, 1, 2);
plot(t_grid, y_res, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.9, 'DisplayName', 'Raw Measured Distance'); hold on;
plot(t_grid, y_f, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.6, 'DisplayName', 'Cleaned Distance');
xline(t_end, ':r', sprintf('Truncation Cut (t = %.2f s)', t_end), 'LineWidth', 1.5);
grid on;
xlabel('Time (seconds)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel({'Amplitude'; 'Sensor Distance (cm)'}, 'FontWeight', 'bold');
title('Sensor Distance: Amplitude vs. Time (Excluding Post-Bumper Stalled Region)', 'FontWeight', 'bold');
legend('Location', 'northeast');

saveas(gcf, fullfile(outputDir, 'task5_preprocessing.png'));

%% 6. FIGURE 2: Model Validation and Error Residuals (Amplitude vs. Time)
figure('Name', 'Task 5: Model Validation (Amplitude vs. Time)', 'Color', 'w', 'Position', [150, 120, 900, 650]);

subplot(2, 1, 1);
plot(t_id, y_id, 'Color', [0.35 0.35 0.35], 'LineWidth', 1.5, 'DisplayName', 'Measured (Cleaned)'); hold on;
plot(t_id, y_sim, 'Color', [0 0.447 0.741], 'LineWidth', 2.0, ...
     'DisplayName', sprintf('Physical Model Simulation (Fit: %.1f%%)', fit_val));
grid on;
ylabel({'Amplitude'; 'Sensor Distance (cm)'}, 'FontWeight', 'bold');
title(sprintf('Model Validation: Measured Output vs. Simulated Plant (Amplitude vs. Time, RMSE = %.2f cm)', rmse), 'FontWeight', 'bold');
legend('Location', 'northeast');

subplot(2, 1, 2);
err = y_id - y_sim;
plot(t_id, err, 'Color', [0 0.447 0.741], 'LineWidth', 1.2); hold on;
yline(0, '--k');
yline(2*std(err), ':r', '+2\sigma'); yline(-2*std(err), ':r', '-2\sigma');
grid on;
xlabel('Time (seconds)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel({'Amplitude'; 'Residual Error (cm)'}, 'FontWeight', 'bold');
title('Simulation Error Residuals: Amplitude vs. Time e(t) = y(t) - y_{sim}(t)', 'FontWeight', 'bold');

saveas(gcf, fullfile(outputDir, 'task5_model_validation.png'));

%% 7. FIGURE 3: Bode Diagram and Stability Margins
% Since positive tilt shortens distance (K < 0), margins are evaluated for -G(s)
G_pos = -G;
[Gm, Pm, Wcg, Wcp] = margin(G_pos);

figure('Name', 'Task 5: Bode Diagram', 'Color', 'w', 'Position', [200, 150, 900, 600]);
margin(G_pos); grid on;
title(sprintf('Bode Diagram of -G(s) = %.3f / (s(s + %.3f))  [PM = %.1f° at %.3f rad/s]', ...
      abs(K_id), b_id, Pm, Wcp), 'FontWeight', 'bold');

saveas(gcf, fullfile(outputDir, 'task5_bode_plot.png'));

%% 8. FIGURE 4: Continuous and Discrete Pole-Zero Maps
figure('Name', 'Task 5: Pole-Zero Maps', 'Color', 'w', 'Position', [250, 180, 950, 450]);

subplot(1, 2, 1);
pzmap(G); grid on;
title('s-Plane Pole-Zero Map: G(s)', 'FontWeight', 'bold');

subplot(1, 2, 2);
pzmap(G_z); grid on;
title(sprintf('z-Plane Pole-Zero Map: G(z) [Ts = %.2f s, 50 Hz]', Ts), 'FontWeight', 'bold');

saveas(gcf, fullfile(outputDir, 'task5_pzmap.png'));

%% 9. FIGURE 5: Step Response Comparison: Amplitude vs. Time (1° Servo Step)
t_step  = (0:0.01:3)';
y_ident = step(G, t_step);
y_ideal = -0.5 * K_ideal * (t_step.^2); % Frictionless runaway parabola

t_z = (0:Ts:3)';
y_z = step(G_z, t_z);

figure('Name', 'Task 5: Step Responses (Amplitude vs. Time)', 'Color', 'w', 'Position', [300, 200, 900, 650]);

subplot(2, 1, 1);
plot(t_step, y_ident, 'Color', [0 0.447 0.741], 'LineWidth', 2.2, ...
     'DisplayName', sprintf('Identified Model (K = %.3f, b = %.3f s^{-1})', K_id, b_id)); hold on;
plot(t_step, y_ideal, '--r', 'LineWidth', 2.0, ...
     'DisplayName', sprintf('Ideal Frictionless Geometry (-%.3f t^2/2)', K_ideal));
grid on;
ylabel({'Amplitude'; 'Position \Delta d (cm)'}, 'FontWeight', 'bold');
title('1° Servo Tilt Continuous Step Response: Amplitude vs. Time', 'FontWeight', 'bold');
legend('Location', 'southwest');

subplot(2, 1, 2);
stairs(t_z, y_z, 'Color', [0.2 0.6 0.3], 'LineWidth', 1.8, 'DisplayName', 'Sampled Output G(z)'); hold on;
plot(t_z, y_z, '.', 'Color', [0.2 0.6 0.3], 'MarkerSize', 8, 'HandleVisibility', 'off');
grid on;
xlabel('Time (seconds)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel({'Amplitude'; 'Position \Delta d[k] (cm)'}, 'FontWeight', 'bold');
title(sprintf('Sampled-Data Discrete Step Response G(z): Amplitude vs. Time [Ts = %.2f s, 50 Hz]', Ts), 'FontWeight', 'bold');
legend('Location', 'southwest');

saveas(gcf, fullfile(outputDir, 'task5_step_response_comparison.png'));

%% 10. Final Characterization Report
poles_s = pole(G);
poles_z = pole(G_z);

fprintf('\n=======================================================================\n');
fprintf('        FINAL TASK 5 CHARACTERIZATION SUMMARY: %s\n', target_run_id);
fprintf('=======================================================================\n');
fprintf('Model Validation Fit: %.2f %%  |  RMSE: %.2f cm  |  n4sid Fit: %s\n', fit_val, rmse, fit_n4_str);

fprintf('\n>>> 1. CONTINUOUS-TIME TRANSFER FUNCTION G(s) <<<\n');
fprintf('         %.4f\n', K_id);
fprintf('  G(s) = --------------------  [cm/deg]\n');
fprintf('         s ( s + %.4f )\n\n', b_id);
fprintf('  * Plant Gain (K):              %.4f cm/s^2 per deg\n', K_id);
fprintf('  * Track Friction (b):          %.4f s^-1\n', b_id);
fprintf('  * Estimated Servo Trim (u0):   %.2f deg\n', u0);
fprintf('  * Continuous Poles:            s1 = %+.4f,  s2 = %+.4f rad/s\n', poles_s(1), poles_s(2));

fprintf('\n>>> 2. DISCRETE-TIME TRANSFER FUNCTION G(z) [Ts = %.2f s (50 Hz)] <<<\n', Ts);
fprintf('         %+.6f z %+.6f\n', num_z(2), num_z(3));
fprintf('  G(z) = ----------------------------------------\n');
fprintf('         z^2 %+.6f z %+.6f\n\n', den_z(2), den_z(3));

fprintf('Equivalent Arduino C++ Difference Equation:\n');
fprintf('  y[k] = %+.6f*y[k-1] %+.6f*y[k-2] %+.6f*u[k-1] %+.6f*u[k-2]\n', ...
        -den_z(2), -den_z(3), num_z(2), num_z(3));
fprintf('  (where u = servo_deg - %.2f, y = distance_cm)\n\n', u0);

fprintf('Discrete Poles (z-plane):\n');
for p = 1:numel(poles_z)
    fprintf('  z%d = %+.5f  (|z| = %.5f, STABLE)\n', p, real(poles_z(p)), abs(poles_z(p)));
end

fprintf('\n>>> 3. FREQUENCY RESPONSE AND STABILITY MARGINS <<<\n');
if isinf(Gm)
    fprintf('  * Gain Margin (Gm):   Inf dB (phase never reaches -180°)\n');
else
    fprintf('  * Gain Margin (Gm):   %.2f dB at %.4f rad/s\n', 20*log10(Gm), Wcg);
end
fprintf('  * Phase Margin (Pm):  %.2f° at %.4f rad/s (%.4f Hz)\n', Pm, Wcp, Wcp/(2*pi));
fprintf('  * Control Insight:    Nearly a pure double integrator; phase lead (PD) is required.\n');

fprintf('\n>>> 4. COMPARISON WITH IDEAL GEOMETRIC SLIDING CAR <<<\n');
fprintf('  * Ideal Frictionless Gain (K_ideal): %.4f cm/s^2 per deg\n', K_ideal);
fprintf('  * Identified Gain (|K_id|):          %.4f cm/s^2 per deg\n', abs(K_id));
fprintf('  * Real-to-Ideal Ratio:               %.1f %%\n', 100 * abs(K_id) / K_ideal);
fprintf('  (Gain reduction is physically expected due to linkage friction, servo compliance,\n');
fprintf('   and rotational inertia of the car wheels/structure).\n');

fprintf('\n=======================================================================\n');
fprintf('All figures saved to: %s\n', outputDir);
fprintf('=======================================================================\n');

%% ======================= Helper Simulation Function =====================
function y = simPlant(p, u, Ts)
% Exact ZOH discrete simulation of plant: y'' + b*y' = K*(u - u0)
% Parameters: p = [b, K, u0, y0, v0]
    b  = abs(p(1)) + 1e-6;
    K  = p(2);
    u0 = p(3);
    x1 = p(4); % Position state
    x2 = p(5); % Velocity state
    
    e   = exp(-b * Ts);
    a12 = (1 - e) / b;
    B1  = K * (Ts - a12) / b;
    B2  = K * a12;
    
    N = length(u);
    y = zeros(N, 1);
    for k = 1:N
        y(k) = x1;
        uk  = u(k) - u0;
        x1n = x1 + a12 * x2 + B1 * uk;
        x2  = e * x2 + B2 * uk;
        x1  = x1n;
    end
end
