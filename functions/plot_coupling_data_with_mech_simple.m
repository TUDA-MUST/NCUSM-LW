function [x_out, y_out, y_out_std] = plot_coupling_data_with_mech_simple(res, FileName_Bases, eval_array_idx, linAx_x_idx, TX_RX_idx, eval_array_type, y_label, timestamp_mech, res_mech)

figure, hold off
x_out = {};
y_out = {};
y_out_std = {};
for fileIdx = 1:length(FileName_Bases)
    lower_limit = timestamp_mech{fileIdx}(1)-60000;
    upper_limit = timestamp_mech{fileIdx}(end)-18000;
    coup_angle = eval("res(fileIdx)." + eval_array_type + "(:,:,:, eval_array_idx, linAx_x_idx, TX_RX_idx);");
    timestamp = res(fileIdx).Timestamps(1:size(coup_angle,1),1:size(coup_angle,2),1:size(coup_angle,3), eval_array_idx, linAx_x_idx, TX_RX_idx);
    timestamp_std = res(fileIdx).Timestamps_std(1:size(coup_angle,1),1:size(coup_angle,2),1:size(coup_angle,3), eval_array_idx, linAx_x_idx, TX_RX_idx);
    %T_0_mech_meas_string = datestr(datetime('2024-3-1') + milliseconds(timestamp), 'yyyy-mm-dd HH:MM:SS.FFF')
    used_data = timestamp>lower_limit & timestamp<upper_limit;
    timestamp_limits = timestamp(used_data);
    timestamp_std_limits = timestamp_std(used_data);
    coup_angle_limits  = coup_angle;
    coup_angle_limits(~used_data)=NaN;
    stress_interp = interp1(timestamp_mech{fileIdx}, res_mech(fileIdx).stress_mech_meas, timestamp);
    stress_interp(isnan(stress_interp)) = 0;
    stress_interp(~used_data)=NaN;
    if all("tt_corr_times" ==eval_array_type)
        size_coup_angle_limits = size(coup_angle_limits);
        coup_angle_limits_flattened_unwraped = unwrap(coup_angle_limits(:)*40000*2*pi)/2/pi/40000;
        coup_angle_limits  = reshape(coup_angle_limits_flattened_unwraped, size_coup_angle_limits);
    end
    if TX_RX_idx == 1
        mean_nr = 1;
    else
        mean_nr = 50;
    end
    hold on
    coup_angle_limits_rs = coup_angle_limits;
																					  
													   
															  
    coup_angle_limits_rs_median = median(coup_angle_limits_rs, 1, "omitmissing");
    coup_angle_limits_rs_std = std(coup_angle_limits_rs,0, 1, "omitmissing");

    if all("tt_corr_times" ==eval_array_type)
        length_of_measure_range = (0.3575-0.03*tand(51.85));
        c_nom = 867.5;
        coup_angle_limits_rs_gv = length_of_measure_range./(coup_angle_limits_rs_median./2+length_of_measure_range/c_nom); %convert transit time to group velocity
        coup_angle_limits_rs_zero = coup_angle_limits_rs_gv;
        coup_angle_limits_rs_std = c_nom^2./2./length_of_measure_range./( coup_angle_limits_rs_median./2./length_of_measure_range.*c_nom +1).^2.*coup_angle_limits_rs_std;
    else
        coup_angle_limits_rs_zero = coup_angle_limits_rs_median- (coup_angle_limits_rs_median(1) -sign(coup_angle_limits_rs_median(1)).*51.85); %so das es mi +-54 anfängt

    end
    stress_interp_rs = stress_interp;
																		  
										  
    stress_interp_mean = mean(stress_interp_rs, 1, "omitmissing");
    stress_interp_mean_crop = reshape(stress_interp_mean(~isnan(stress_interp_mean)),1,[]);

    coup_angle_limits_rs_zero_crop = reshape(coup_angle_limits_rs_zero(~isnan(coup_angle_limits_rs_zero)),1,[]);
    coup_angle_limits_rs_std_crop = reshape(coup_angle_limits_rs_std(~isnan(coup_angle_limits_rs_std)),1,[]);
    if mean_nr==1
        errorbar(stress_interp_mean_crop./1e6.*1.06, coup_angle_limits_rs_zero_crop, coup_angle_limits_rs_zero_crop.*0+ 0.04, "x") %standard deviation from reference measurement.
        y_out_std{fileIdx} = coup_angle_limits_rs_zero_crop.*0+ 0.04;
    else
        errorbar(stress_interp_mean_crop./1e6.*1.06, coup_angle_limits_rs_zero_crop, coup_angle_limits_rs_std_crop,  "x")
        y_out_std{fileIdx} = coup_angle_limits_rs_std_crop;
    end
    x_out{fileIdx} = stress_interp_mean_crop./1e6.*1.06;
    y_out{fileIdx} = coup_angle_limits_rs_zero_crop;
end
xlabel("Stress (MPa)")
ylabel(y_label)
end