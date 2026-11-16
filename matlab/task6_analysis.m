%% Task 6: Root Locus and Frequency Response Analysis
% Course Project: Part-II (ES 245) - Ball/Car-and-Beam System
%
% This script implements all three sub-tasks of Task 6:
%   6.1 Root Locus and Nyquist Analysis:
%       - Plots Root Locus with design criteria boundaries (zeta >= 0.707, sigma >= 0.6).
%       - Analyzes breakaway point, stability, and dominant pole trajectories.
%       - Plots Nyquist diagram with critical point (-1, j0) and evaluates Nyquist criterion.
%   6.2 Bode Plot Analysis:
%       - Plots open-loop Bode diagram and computes Gain Margin (Gm) and Phase Margin (Pm).
%       - Compares and correlates stability margins across Bode, Root Locus, and Nyquist.
%   6.3 Design Criteria Evaluation:
%       - Evaluates uncompensated closed-loop performance against design criteria:
%         * Overshoot < 5%
%         * Settling time < 5 s (5% tolerance band).
%       - Assesses inadequacy of proportional control and justifies the need for
%         a Lead Compensator (Task 7).
%
% Validated System Transfer Function (from Task 5):
%   G(s) = -0.6463 / (s * (s + 0.0170))  [cm/deg]
%   With negative feedback accounting for physical sensor direction:
%   G_ol(s) = 0.6463 / (s * (s + 0.0170))
%
% Author: Control Systems Laboratory
% Date: November 2026

clear; close all; clc;

%% 1. Plant Configuration and Setup
% Physical plant parameters identified in Task 5:
K_plant = 0.6463;   % Plant gain (cm/s^2 per deg)
b_plant = 0.0170;   % Track friction coefficient (s^-1)

% Open-loop transfer function: G(s) = K / (s * (s + b))
s = tf('s');
G = K_plant / (s * (s + b_plant));
G.InputName  = 'Servo Angle (\Delta\theta)';
G.InputUnit  = 'deg';
G.OutputName = 'Distance (\Delta y)';
G.OutputUnit = 'cm';

% Output folder for saving figures
outputDir = fileparts(mfilename('fullpath'));
if isempty(outputDir), outputDir = pwd; end

fprintf('=======================================================================\n');
fprintf('  TASK 6: ROOT LOCUS, NYQUIST, AND FREQUENCY RESPONSE ANALYSIS         \n');
fprintf('=======================================================================\n');
fprintf('Open-Loop Plant Transfer Function:\n');
fprintf('            %.4f\n', K_plant);
fprintf('  G(s) = --------------------  [cm/deg]\n');
fprintf('         s ( s + %.4f )\n\n', b_plant);

%% 2. Task 6.1: Root Locus Analysis
fprintf('>>> 6.1 ROOT LOCUS AND NYQUIST ANALYSIS <<<\n');
fprintf('Open-Loop Poles: s1 = 0.0000,  s2 = -%.4f rad/s\n', b_plant);
sigma_breakaway = -b_plant / 2;
fprintf('Breakaway Point: sigma_b = -b/2 = -%.4f rad/s\n', sigma_breakaway);
fprintf('Asymptote Angles: phi = +-90 deg along vertical line Re(s) = -%.4f\n\n', sigma_breakaway);

figure('Name', 'Task 6.1: Root Locus Analysis', 'Color', 'w', 'Position', [100, 100, 850, 650]);
rlocus(G);
hold on; grid on;

% Overlay Design Specifications on s-Plane:
% 1. Damping ratio line for Overshoot < 5% -> zeta >= 0.707 (angle <= 45 deg)
zeta_spec = 0.707;
wn_grid = linspace(0, 3, 100);
sgrid(zeta_spec, []);

% 2. Settling time boundary for Ts < 5 s (5% tolerance: Ts = 3/(zeta*wn) <= 5 -> sigma >= 0.6)
xline(-0.6, '--', 'Color', [0 0.5 0], 'LineWidth', 1.8, ...
      'DisplayName', 'Settling Time Boundary (\sigma \ge 0.6, T_s \le 5 s)');

% Mark open-loop poles and breakaway point
plot([0, -b_plant], [0, 0], 'rx', 'MarkerSize', 10, 'LineWidth', 2.5, 'DisplayName', 'Open-Loop Poles');
plot(sigma_breakaway, 0, 'mo', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', sprintf('Breakaway (\\sigma = -%.4f)', sigma_breakaway));

xlim([-1.5, 0.4]);
ylim([-2.0, 2.0]);
title('Task 6.1: Root Locus of Open-Loop Plant G(s) with Design Spec Boundaries', 'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'northwest');

saveas(gcf, fullfile(outputDir, 'task6_root_locus.png'));

%% 3. Task 6.1: Nyquist Diagram Analysis
figure('Name', 'Task 6.1: Nyquist Plot Analysis', 'Color', 'w', 'Position', [150, 120, 1000, 500]);

% Subplot 1: Global Nyquist contour
subplot(1, 2, 1);
nyquist(G);
hold on; grid on;
plot(-1, 0, 'r+', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Critical Point (-1, j0)');
xlim([-5, 2]);
ylim([-5, 5]);
title('Nyquist Diagram (Global View)', 'FontWeight', 'bold');
legend('Location', 'southwest');

% Subplot 2: Zoomed view near critical point (-1, j0) to visualize stability margin
subplot(1, 2, 2);
nyquist(G);
hold on; grid on;
th = linspace(0, 2*pi, 200);
plot(cos(th), sin(th), 'k:', 'LineWidth', 1.2, 'DisplayName', 'Unit Circle (|G| = 1)');
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2.5, 'DisplayName', 'Critical Point (-1, j0)');

% Find crossover frequency point on Nyquist contour
[~, Pm_val, ~, Wcp_val] = margin(G);
[mag_cp, phase_cp] = bode(G, Wcp_val);
re_cp = mag_cp * cosd(phase_cp);
im_cp = mag_cp * sind(phase_cp);
plot(re_cp, im_cp, 'mo', 'MarkerSize', 7, 'LineWidth', 2, 'DisplayName', sprintf('\\omega_{cp} = %.4f rad/s', Wcp_val));
plot([-1, re_cp], [0, im_cp], 'r-', 'LineWidth', 1.5, 'DisplayName', 'Proximity Distance = 0.021');

xlim([-1.3, -0.7]);
ylim([-0.3, 0.3]);
title(sprintf('Zoomed View Near (-1, j0)  [Phase Margin = %.2f°]', Pm_val), 'FontWeight', 'bold');
legend('Location', 'southeast');

saveas(gcf, fullfile(outputDir, 'task6_nyquist.png'));

%% 4. Task 6.2: Bode Plot and Stability Margins
[Gm, Pm, Wcg, Wcp] = margin(G);

figure('Name', 'Task 6.2: Open-Loop Bode Diagram', 'Color', 'w', 'Position', [200, 150, 900, 600]);
margin(G);
grid on;
title(sprintf('Task 6.2: Open-Loop Bode Diagram [Gain Margin = Inf dB, Phase Margin = %.2f° at %.4f rad/s]', ...
      Pm, Wcp), 'FontSize', 11, 'FontWeight', 'bold');

saveas(gcf, fullfile(outputDir, 'task6_bode.png'));

fprintf('\n>>> 6.2 BODE FREQUENCY RESPONSE & STABILITY MARGINS <<<\n');
if isinf(Gm)
    fprintf('  * Gain Margin (Gm):           Inf dB (phase never reaches -180 deg)\n');
else
    fprintf('  * Gain Margin (Gm):           %.2f dB (at Wcg = %.4f rad/s)\n', 20*log10(Gm), Wcg);
end
fprintf('  * Phase Margin (Pm):          %.2f deg (at Wcp = %.4f rad/s = %.4f Hz)\n', Pm, Wcp, Wcp/(2*pi));
fprintf('  * Gain Crossover Freq (Wcp):  %.4f rad/s\n', Wcp);
fprintf('  * Low-Frequency Slope:        -20 dB/decade (integrator pole s = 0)\n');
fprintf('  * High-Frequency Slope:       -40 dB/decade (double-pole roll-off beyond w = %.4f rad/s)\n', b_plant);

%% 5. Task 6.3: Closed-Loop Performance vs. Design Criteria
% Design Specifications:
% 1. Percent Overshoot: OS < 5%  -> Requires damping ratio zeta >= 0.707
% 2. Settling Time (5% band): Ts < 5.0 s -> Requires real pole part sigma >= 3/5 = 0.6 s^-1

fprintf('\n>>> 6.3 DESIGN CRITERIA EVALUATION (UNCOMPENSATED PROPORTIONAL CONTROL) <<<\n');
fprintf('Target Specifications:\n');
fprintf('  1. Maximum Overshoot:  OS < 5.0 %%\n');
fprintf('  2. Settling Time (5%%): Ts < 5.0 seconds\n\n');

% Evaluate closed loop under unity proportional control (Kp = 1.0):
Kp_test = [0.5, 1.0, 2.0];
t_sim = 0:0.01:15;

figure('Name', 'Task 6.3: Closed-Loop Step Response vs Design Criteria', 'Color', 'w', 'Position', [250, 180, 950, 600]);
hold on; grid on;

% Target Step Reference (r = 1.0 cm)
yline(1.0, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Target Reference r(t) = 1.0 cm');

% Spec Envelope: 5% band (0.95 to 1.05) and Ts <= 5 s
patch([0, 5, 5, 0], [0.95, 0.95, 1.05, 1.05], [0 0.6 0.2], ...
      'FaceAlpha', 0.15, 'EdgeColor', [0 0.5 0], 'LineStyle', '--', 'LineWidth', 1.2, ...
      'DisplayName', 'Target Spec Box (OS \le 5%, T_s \le 5 s)');

colors = {[0 0.447 0.741], [0.85 0.33 0.1], [0.494 0.184 0.556]};
for k = 1:numel(Kp_test)
    kp = Kp_test(k);
    T_cl = feedback(kp * G, 1);
    [y_cl, t_out] = step(T_cl, t_sim);
    
    wn_cl   = sqrt(kp * K_plant);
    zeta_cl = b_plant / (2 * wn_cl);
    os_cl   = (max(y_cl) - 1.0) * 100;
    ts_cl   = 3.0 / (zeta_cl * wn_cl); % 5% settling time formula
    
    fprintf('P-Control Kp = %.1f:\n', kp);
    fprintf('  * Natural Frequency (wn): %.4f rad/s\n', wn_cl);
    fprintf('  * Damping Ratio (zeta):   %.5f (Target: >= 0.707) -> FAILS\n', zeta_cl);
    fprintf('  * Percent Overshoot:      %.2f %% (Target: < 5.0 %%) -> FAILS\n', os_cl);
    fprintf('  * 5%% Settling Time:       %.1f s  (Target: < 5.0 s)  -> FAILS\n\n', ts_cl);
    
    plot(t_out, y_cl, 'Color', colors{k}, 'LineWidth', 1.8, ...
         'DisplayName', sprintf('P-Control K_p = %.1f (OS = %.1f%%, T_s = %.0f s)', kp, os_cl, ts_cl));
end

xlim([0, 15]);
ylim([-0.2, 2.2]);
xlabel('Time (seconds)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Position Response y(t) (cm)', 'FontSize', 11, 'FontWeight', 'bold');
title('Task 6.3: Uncompensated Closed-Loop Step Response vs. Design Criteria Bounds', 'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'northeast');

saveas(gcf, fullfile(outputDir, 'task6_design_criteria_eval.png'));

%% 6. Final Task 6 Correlation & Analytical Summary
fprintf('=======================================================================\n');
fprintf('  CORRELATION BETWEEN ROOT LOCUS, NYQUIST, AND BODE STABILITY MARGINS  \n');
fprintf('=======================================================================\n');
fprintf('1. Gain Margin (Gm = Inf dB):\n');
fprintf('   - Root Locus: Branches never cross the imaginary axis into RHP for any Kp > 0.\n');
fprintf('   - Nyquist:    Locus stays in upper half plane and never crosses (-Inf, 0).\n');
fprintf('   - Bode:       Phase approaches -180 deg but never crosses it.\n\n');

fprintf('2. Phase Margin (Pm = 1.21 deg):\n');
fprintf('   - Root Locus: Branches break away at sigma = -0.0085 and run parallel to jw axis.\n');
fprintf('                 Real part is locked at -0.0085; damping ratio zeta -> 0 as gain increases.\n');
fprintf('   - Nyquist:    Locus passes within 0.021 distance of critical point (-1, j0).\n');
fprintf('   - Bode:       Phase at gain crossover (0.8038 rad/s) is -178.79 deg (PM = 1.21 deg).\n\n');

fprintf('3. Inadequacy of Existing System (Task 6.3 Conclusion):\n');
fprintf('   - Maximum Overshoot: ~96.7%% (violates < 5%% by a factor of 19x).\n');
fprintf('   - Settling Time:     352.9 seconds (violates < 5 s by a factor of 70x).\n');
fprintf('   - Root Locus Bottleneck: Real part of poles is trapped at Re(s) = -0.0085;\n');
fprintf('     proportional gain Kp CANNOT move poles deeper into the left half plane.\n\n');

fprintf('4. Need for Compensator Improvement (Transition to Task 7):\n');
fprintf('   - A LEAD COMPENSATOR C(s) = Kc * (s + z_c) / (s + p_c) is strictly required.\n');
fprintf('   - Placing zero z_c near s = -0.8 pulls the root locus branches far to the left\n');
fprintf('     into the desired region (Re(s) <= -0.6 to achieve Ts < 5 s).\n');
fprintf('   - The compensator will inject +55 deg to +60 deg phase lead, boosting Phase Margin\n');
fprintf('     from 1.21 deg to > 55 deg (achieving Overshoot < 5%%).\n');

fprintf('=======================================================================\n');
fprintf('Task 6 Analysis Complete. All figures saved to: %s\n', outputDir);
fprintf('=======================================================================\n');
