function [mean_all, std_all, mean_residual_all, std_residual_all, mean_residual_rel_all, std_residual_rel_all, max_residual_all] = calculate_metrics(x_out_tx, y_out_tx, y_out_std_tx, label)

[mean_vals, std_vals] = cellfun(@mean_and_std_rel, y_out_tx, y_out_std_tx);
mean_all = mean(mean_vals);
std_all = mean(std_vals);
[mean_residual, max_residual, std_residual, mean_residual_rel, std_residual_rel] = cellfun(@analyze_deviation_from_line, x_out_tx, y_out_tx);
mean_residual_all = mean(mean_residual);
max_residual_all = max(max_residual);
std_residual_all = mean(std_residual);
mean_residual_rel_all = mean(mean_residual_rel);
std_residual_rel_all = mean(std_residual_rel);

end
function [mean_residual,max_residual, std_residual, mean_residual_rel, std_residual_rel] = analyze_deviation_from_line(x_out, y_out)
label= "XX";
% Perform linear regression to fit a straight line
coeffs = polyfit(x_out, y_out, 1); % Linear fit (degree 1 polynomial)
y_fit = polyval(coeffs, x_out);   % Evaluate the fitted line at x_out

% Calculate residuals (deviation from the straight line)
residuals = (y_out - y_fit);
residuals_rel = ((y_out - y_fit))./y_fit;

% Compute metrics for residuals
mean_residual = mean(abs(residuals));
max_residual = max(abs(residuals));
std_residual = std(abs(residuals));
mean_residual_rel = mean(abs(residuals_rel));
std_residual_rel = std(abs(residuals_rel));
end
function [mean_rel, std_rel] = mean_and_std_rel(y, y_std)
mean_rel = mean(y_std./abs(y));
std_rel = std(y_std./y);
end
