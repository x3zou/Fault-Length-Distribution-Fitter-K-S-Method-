%TEST_FAULT_KS  Demonstrate fit_fault_ks on a synthetic fault-length dataset.
%
%  The synthetic data mimics a typical field observation:
%   • LOW  end  (x < 1.0) : curve flattens  (slope ~ -0.4) — under-sampling
%   • MIDDLE    (1.0–3.0) : true power-law  y = -2*x + 8   — target zone
%   • HIGH end  (x > 3.0) : curve steepens (slope ~ -5)    — truncation
%
%  Run this script as-is; it will call fit_fault_ks and produce a figure.

clear; clc; close all;
rng(2024);                    % reproducible noise

% ═══════════════════════════════════════════════════════════════════════════════
%  1.  GENERATE SYNTHETIC DATASET
% ═══════════════════════════════════════════════════════════════════════════════

b_true  = 2.0;
c_true  = 8.0;
sigma   = 0.12;               % Gaussian noise standard deviation

% ── Low end: flattened (under-sampled small faults) ───────────────────────────
n_low  = 14;
x_low  = linspace(0.10, 0.95, n_low)';
% Gradually transition: slope starts near 0 and approaches -b_true at x=1
t      = (x_low - 0.10) / (0.95 - 0.10);          % 0→1
slope_low = -0.4;
y_low  = (-b_true*1.0 + c_true) ...                % anchor at x=1
         + slope_low .* (x_low - 1.0) ...          % deviated slope
         + sigma * randn(n_low, 1);

% ── Middle: true power-law ────────────────────────────────────────────────────
n_mid  = 35;
x_mid  = linspace(1.00, 3.00, n_mid)';
y_mid  = -b_true .* x_mid + c_true + sigma * randn(n_mid, 1);

% ── High end: steepened (truncated large faults) ──────────────────────────────
n_high = 10;
x_high = linspace(3.05, 3.80, n_high)';
slope_high = -5.5;
y_high = (-b_true*3.0 + c_true) ...
         + slope_high .* (x_high - 3.0) ...
         + sigma * randn(n_high, 1);

% ── Assemble and sort ─────────────────────────────────────────────────────────
x_data = [x_low;  x_mid;  x_high];
y_data = [y_low;  y_mid;  y_high];
[x_data, ord] = sort(x_data);
y_data = y_data(ord);

fprintf('Dataset: %d points total  (%d low / %d mid / %d high)\n', ...
        length(x_data), n_low, n_mid, n_high);
fprintf('True model: y = -%.1f * x + %.1f\n\n', b_true, c_true);

% ═══════════════════════════════════════════════════════════════════════════════
%  2.  RUN THE AUTOMATIC FITTING
% ═══════════════════════════════════════════════════════════════════════════════

[b_fit, c_fit, x_fit, y_fit, idx_lo, idx_hi] = ...
    fit_fault_ks(x_data, y_data, 0.25, 0.05);

% ═══════════════════════════════════════════════════════════════════════════════
%  3.  PLOT RESULTS
% ═══════════════════════════════════════════════════════════════════════════════

fig = figure('Name', 'Fault-Length Distribution Fit', ...
             'Color', 'w', 'Position', [100 80 900 640]);

% ── Main log-log plot ─────────────────────────────────────────────────────────
ax1 = subplot(2, 2, [1 2]);
hold(ax1, 'on');

% Shade excluded regions
lo_x = x_data([1, idx_lo, idx_lo, 1]);
lo_y = [min(y_data)-1, min(y_data)-1, max(y_data)+1, max(y_data)+1];
fill(ax1, lo_x, lo_y, [1 0.85 0.75], 'EdgeColor', 'none', 'FaceAlpha', 0.55);

hi_x = x_data([idx_hi, end, end, idx_hi]);
fill(ax1, hi_x, lo_y, [0.75 0.85 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.55);

% All data
h_all = scatter(ax1, x_data, y_data, 40, [0.6 0.6 0.6], 'o', ...
                'LineWidth', 0.8, 'DisplayName', 'All data');

% Fitted window
h_fit = scatter(ax1, x_fit, y_fit, 55, [0.1 0.45 0.8], 'o', ...
                'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, ...
                'DisplayName', sprintf('Fitted window (n = %d)', length(x_fit)));

% True model
x_plt  = linspace(min(x_data), max(x_data), 200);
h_true = plot(ax1, x_plt, -b_true.*x_plt + c_true, 'k--', ...
              'LineWidth', 1.5, 'DisplayName', ...
              sprintf('True: y = -%.1f x + %.1f', b_true, c_true));

% Fitted model (extended over full x range for visibility)
h_est  = plot(ax1, x_plt, -b_fit.*x_plt + c_fit, 'r-', ...
              'LineWidth', 2.2, 'DisplayName', ...
              sprintf('Fitted: y = -%.3f x + %.3f', b_fit, c_fit));

% Cutoff lines
xline(ax1, x_data(idx_lo),  '--', 'Color', [0.8 0.3 0], 'LineWidth', 1.4, ...
      'Label', sprintf('Low cut  x = %.2f', x_data(idx_lo)), ...
      'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
xline(ax1, x_data(idx_hi),  '--', 'Color', [0.0 0.3 0.8], 'LineWidth', 1.4, ...
      'Label', sprintf('High cut  x = %.2f', x_data(idx_hi)), ...
      'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

legend(ax1, [h_all, h_fit, h_true, h_est], 'Location', 'northeast', ...
       'FontSize', 9);
xlabel(ax1, 'log_{10}(Fault Length)');
ylabel(ax1, 'log_{10}(Frequency)');
title(ax1, 'Fault-Length Distribution: Automatic Range Selection (K-S Test)');
ylim(ax1, [min(y_data)-0.3, max(y_data)+0.3]);
grid(ax1, 'on');

% Annotation boxes
annotation('textbox', [0.13 0.63 0.13 0.07], ...
    'String', 'Under-sampled\n(flat)', ...
    'FitBoxToText', 'on', 'BackgroundColor', [1 0.85 0.75], ...
    'EdgeColor', 'none', 'FontSize', 8, 'HorizontalAlignment', 'center');
annotation('textbox', [0.73 0.63 0.13 0.07], ...
    'String', 'Truncated\n(steep)', ...
    'FitBoxToText', 'on', 'BackgroundColor', [0.75 0.85 1.0], ...
    'EdgeColor', 'none', 'FontSize', 8, 'HorizontalAlignment', 'center');

% ── Residual plot for fitted window ──────────────────────────────────────────
ax2 = subplot(2, 2, 3);
resid_fit = y_fit - (-b_fit .* x_fit + c_fit);
scatter(ax2, x_fit, resid_fit, 45, [0.1 0.45 0.8], 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.4);
yline(ax2, 0, 'k-', 'LineWidth', 1.2);
yline(ax2,  2*std(resid_fit), 'r:', 'LineWidth', 1.2);
yline(ax2, -2*std(resid_fit), 'r:', 'LineWidth', 1.2, 'Label', '±2σ');
xlabel(ax2, 'log_{10}(Fault Length)');
ylabel(ax2, 'Residual');
title(ax2, 'Residuals (Fitted Window)');
grid(ax2, 'on');

% ── Q-Q plot of residuals ─────────────────────────────────────────────────────
ax3 = subplot(2, 2, 4);
z = (resid_fit - mean(resid_fit)) / std(resid_fit);
qqplot(z);                    % built-in Q-Q plot vs standard normal
title(ax3, 'Q-Q Plot of Standardised Residuals');
xlabel(ax3, 'Standard Normal Quantiles');
ylabel(ax3, 'Sample Quantiles');
grid(ax3, 'on');

% ── Super-title ───────────────────────────────────────────────────────────────
sgtitle(fig, sprintf('b_{true} = %.2f  →  b_{fit} = %.4f     c_{true} = %.2f  →  c_{fit} = %.4f', ...
        b_true, b_fit, c_true, c_fit), 'FontSize', 11, 'FontWeight', 'bold');

% ═══════════════════════════════════════════════════════════════════════════════
%  4.  PRINT THE TEST DATASET
% ═══════════════════════════════════════════════════════════════════════════════

fprintf('=== Synthetic test dataset (x, y) ===\n');
fprintf('  Index   x_data    y_data   Zone\n');
zones = [repmat("LOW", n_low,1); repmat("MID", n_mid,1); repmat("HIGH", n_high,1)];
zones = zones(ord);
for k = 1 : length(x_data)
    marker = '';
    if k == idx_lo, marker = '  ← low  cut'; end
    if k == idx_hi, marker = '  ← high cut'; end
    fprintf('  %3d   %7.4f   %7.4f   %s%s\n', ...
            k, x_data(k), y_data(k), zones(k), marker);
end