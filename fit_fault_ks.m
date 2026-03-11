function [b, c, x_fit, y_fit, idx_low, idx_high] = fit_fault_ks(x, y, min_frac, alpha)
%FIT_FAULT_KS  Automatically determine the best-fit range for a fault-length
%              frequency distribution in log-log space using the K-S test.
%
%  The distribution y = -b*x + c is linear (power-law) over the "valid" range.
%  At the LOW end  the curve FLATTENS  (under-sampling of small faults).
%  At the HIGH end the curve STEEPENS  (truncation of large faults).
%
%  ALGORITHM
%  ----------
%  For every candidate window [i, j] (subject to a minimum length):
%    1. Fit a straight line  y = p1*x + p0  by least squares.
%    2. Standardise the residuals to zero mean / unit variance.
%    3. Apply the one-sample K-S test against N(0,1).
%       A high p-value means the residuals are consistent with normality,
%       i.e. the data in this window genuinely follow a straight line.
%    4. Record (i, j, p-value, K-S statistic, window length).
%  Selection rule:
%    • Among all windows whose K-S p-value >= alpha  (null not rejected),
%      keep the LONGEST one  (most data used).
%    • If no window passes at the chosen alpha, relax and return the
%      window with the single highest p-value, with a warning.
%
%  USAGE
%  -----
%    [b, c, x_fit, y_fit, idx_low, idx_high] = fit_fault_ks(x, y)
%    [b, c, x_fit, y_fit, idx_low, idx_high] = fit_fault_ks(x, y, min_frac, alpha)
%
%  INPUTS
%    x, y      – log10(length) and log10(frequency) vectors (sorted by x,
%                ascending).  May be row or column vectors.
%    min_frac  – minimum fraction of total points a window must contain
%                (default 0.25).  Prevents trivially short windows.
%    alpha     – K-S test significance level (default 0.05).
%
%  OUTPUTS
%    b, c      – fitted slope and intercept:  y = -b*x + c   (b > 0)
%    x_fit     – x values used for the final fit
%    y_fit     – y values used for the final fit
%    idx_low   – first index of the accepted window (in the original vectors)
%    idx_high  – last  index of the accepted window (in the original vectors)
%
%  DEPENDENCIES
%    Requires the Statistics and Machine Learning Toolbox (kstest).
%
%  EXAMPLE
%    run test_fault_ks.m

% ── Input handling ────────────────────────────────────────────────────────────
if nargin < 3 || isempty(min_frac), min_frac = 0.25; end
if nargin < 4 || isempty(alpha),    alpha    = 0.05; end

x = x(:);   % ensure column vectors
y = y(:);
n = length(x);

if n < 6
    error('fit_fault_ks: need at least 6 data points.');
end

min_win = max(5, round(n * min_frac));   % absolute minimum window size

% ── Exhaustive window search ──────────────────────────────────────────────────
% Pre-allocate result table: [i_start, i_end, p_value, ks_D, window_length]
max_windows = round(n*(n-1)/2);
results = zeros(max_windows, 5);
row = 0;

for i = 1 : n - min_win + 1
    for j = i + min_win - 1 : n

        xi = x(i:j);
        yi = y(i:j);

        % ---- linear fit ----
        p       = polyfit(xi, yi, 1);
        resid   = yi - polyval(p, xi);

        % ---- standardise residuals ----
        s = std(resid);
        if s < eps          % perfectly collinear – treat as ideal
            ks_stat = 0;
            pval    = 1;
        else
            z = (resid - mean(resid)) / s;
            [~, pval, ks_stat] = kstest(z);   % test vs N(0,1)
        end

        row = row + 1;
        results(row, :) = [i, j, pval, ks_stat, j - i + 1];
    end
end
results = results(1:row, :);   % trim unused rows

% ── Window selection ──────────────────────────────────────────────────────────
pass_mask = results(:, 3) >= alpha;

if any(pass_mask)
    passing     = results(pass_mask, :);
    [~, best]   = max(passing(:, 5));          % longest passing window
    best_row    = passing(best, :);
else
    [~, best]   = max(results(:, 3));          % best available p-value
    best_row    = results(best, :);
    warning('fit_fault_ks:noPass', ...
        'No window passed the K-S test at alpha = %.3f. '  + ...
        'Returning the window with the highest p-value (%.4f).', ...
        alpha, results(best, 3));
end

idx_low  = best_row(1);
idx_high = best_row(2);

% ── Final fit on selected window ──────────────────────────────────────────────
x_fit = x(idx_low : idx_high);
y_fit = y(idx_low : idx_high);

p_final = polyfit(x_fit, y_fit, 1);
b =  -p_final(1);          % convention: slope stored as positive b
c =   p_final(2);

% ── Diagnostic print ─────────────────────────────────────────────────────────
fprintf('\n--- fit_fault_ks result ---\n');
fprintf('  Accepted window : indices %d – %d  (%d / %d points)\n', ...
        idx_low, idx_high, idx_high - idx_low + 1, n);
fprintf('  x range         : [%.3f, %.3f]\n', x(idx_low), x(idx_high));
fprintf('  Fitted line     : y = -%.4f * x + %.4f\n', b, c);
fprintf('  K-S p-value     : %.4f\n', best_row(3));
fprintf('  K-S D statistic : %.4f\n\n', best_row(4));
end