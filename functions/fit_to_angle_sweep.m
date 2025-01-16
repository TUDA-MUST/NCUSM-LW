function [coup_angle_sinc, rmse, c_std] = fit_to_angle_sweep(array_angles,lw_amps_all_ch, factor_crop, do_plot)

x_fit = sind(array_angles);
y_fit = lw_amps_all_ch;
y_idx = 1:length(y_fit);

% Fit a sinc function to the data
sincFunc = @(a, b, c, x) a * abs(sinc(b * (x - c)));
%old but error prone method: [max_y, max_idx] = max(y_fit);
max_y = quantile(y_fit,0.85);
max_idx = round(mean(y_idx(y_fit>max_y)));
lowest_y = quantile(y_fit,factor_crop);
%factor_crop = 0.5;
startIdx = find(y_fit>=lowest_y,1);
endIdx = find(y_fit>=lowest_y,1, "last");
y_fit_crop = y_fit(startIdx:endIdx);
x_fit_crop = x_fit(startIdx:endIdx);
[sincModel, gof] = fit(x_fit_crop', y_fit_crop, sincFunc, 'StartPoint', [max_y, 1/4, x_fit(max_idx)],  'Normalize', 'on', "Robust", "Bisquare");
rmse = gof.rmse;
% Fit a parabola to the data
%parabolaModel = fit(x_fit_crop', y_fit_crop', 'poly2');
% Find the position of the maximum in the parabola fit
%parabolaMaxPos = fminbnd(@(x) -parabolaModel(x), min(x_fit_crop), max(x_fit_crop));

% Find the position of the maximum in the sinc fit
sincMaxPos = fminbnd(@(x) -sincModel(x), min(x_fit_crop), max(x_fit_crop));

% Get confidence intervals for the fit parameters
ci = confint(sincModel);
% Extract the uncertainty (half-width of the confidence interval)
c_uncertainty = (ci(2,3) - ci(1,3)) / 2;  % Assuming 'c' is the third parameter
c_std = c_uncertainty/2;

%coup_angle_parab = asind(parabolaMaxPos)
coup_angle_sinc = asind(sincMaxPos);
if do_plot
%figure,
%hold off
%plot(asind(x_fit_crop'), y_fit_crop'./max_y, 'x');
plot(asind(x_fit'), y_fit'./max(y_fit), 'x');
hold on;
%plot(asind(x_fit_crop), sincModel(x_fit_crop)./max_y);
%plot(asind(x_fit_crop), sincModel(x_fit_crop));

plot(asind(x_fit), sincModel(x_fit)./max(y_fit));
plot([asind(sincMaxPos),asind(sincMaxPos)], [0,1], 'k--')
legend('Measured', 'Fit', 'Opt. angle', "Location", "South");
xlabel("TX Angle (°)")
ylabel("Amplitude (norm.)")
axis tight
%ylim([factor_crop,1.025])
end
end