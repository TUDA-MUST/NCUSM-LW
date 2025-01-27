clear all
% Warning: This script might take many hours to run. Evaluation results are already saved in the dataset.
% To view only the plots, execute the script until the desired plot is generated.

do_plot = true; % Set to false to disable plots and improve performance

addpath("functions") % Add the directory containing helper functions

% Four files corresponding to different runs of a hydraulic press
FileName_Bases(1) = "data\us_data\2024-11-19-IFSW_zwei_vent_top1";
FileName_Bases(2) = "data\us_data\2024-11-19-IFSW_zwei_vent_top2";
FileName_Bases(3) = "data\us_data\2024-11-19-IFSW_zwei_vent_top3";
FileName_Bases(4) = "data\us_data\2024-11-19-IFSW_zwei_vent_top4";

for fileIdx_adfs = 1:length(FileName_Bases) % Iterate through all measurement runs
    clear signals_a_ldv_ref % Clear reference signals to avoid interference between runs
    FileName_Base = FileName_Bases(fileIdx_adfs);

    % Load metadata for the current run
    load(FileName_Base + "_meta.mat", "start_timestamp", "script_that_generated_this_data", "outer_Repetition_Names")

    % Determine the number of valid repetition files
    TX_RX_file_extensions = ["TX", "RX1"];
    outer_Repetition_Names = 1:1000; % Check up to 1000 repetitions
    for outRepIdx_first = 1:length(outer_Repetition_Names)
        if ~isfile(FileName_Base + string(outer_Repetition_Names(outRepIdx_first)) + TX_RX_file_extensions(1) + ".mat")
            highest_number = outRepIdx_first;
            break
        end
    end
    outer_Repetition_Names = 1:(outRepIdx_first-1);
    save(FileName_Base + "_meta.mat", "start_timestamp", "script_that_generated_this_data", "outer_Repetition_Names", '-v7.3')

    %% Initialize arrays for data storage
    outer_Repetition_Names_old  = 1; % Compatibility with previous datasets
    inner_Repetition_Names = 1:5; % Number of repetitions per outer cycle
    linAx_x = 1; linAx_x_idx = 1; % Placeholder for future use
    eval_array_n = 0:1; % Sweep parameter for evaluations

    % Arrays for storing measurement data
    TX_coup_angles = NaN(length(inner_Repetition_Names), length(outer_Repetition_Names_old), length(outer_Repetition_Names), length(eval_array_n), length(linAx_x), length(TX_RX_file_extensions));
    TX_c_std_all = TX_coup_angles;
    RX_coup_angles = TX_coup_angles;
    tt_corr_times = TX_coup_angles;
    tt_thres_times = TX_coup_angles;
    RX_max_idxs = TX_coup_angles;
    Timestamps = TX_coup_angles;
    Timestamps_std = TX_coup_angles;
    Temperatures = NaN(length(outer_Repetition_Names), length(eval_array_n), length(TX_RX_file_extensions));
    Temperatures_std = Temperatures;

    % Start looping through all measurement points
    for outRepIdx = 1:length(outer_Repetition_Names)
        for TX_RX_idx = 1:length(TX_RX_file_extensions)

            % Load raw data for each measurement
            load(FileName_Base + string(outer_Repetition_Names(outRepIdx)) + TX_RX_file_extensions(TX_RX_idx) + ".mat", "array_samples", "LDV_datas", "linAx_positions_x", "array_angles", "inner_Repetition_Names", "N_averages", "Fs_LDV", "Fs_array", "acquisition_time_s", "N_Samples_cutoff", "using_rec_array", "using_LDV", "using_arduino_env", "meas_timestamps", "pressures", "humidities", "temperatures", "creation_date")

            for eval_array_idx = 1:length(eval_array_n)
                eval_array = logical(eval_array_n(eval_array_idx)); % Evaluate array-based signals
                eval_ldv = ~eval_array; % Evaluate LDV-based signals

                % Define signal parameters
                f_0 = 40000; % Signal frequency (Hz)
                c = 343; % Speed of sound (m/s)
                lambda = c/f_0; % Wavelength (m)

                % Set sampling rate and measurement coordinates
                if eval_ldv
                    f_s = Fs_LDV; % LDV sampling rate
                    y_tr = linAx_positions_x' ./ 1000 ./ 1.61; % Normalize positions
                    x_tr = y_tr .* 0; % LDV uses a fixed axis
                else
                    f_s = Fs_array; % Array sampling rate
                    x_tr=[-0.00643125; -0.00214375;0.00214375;0.00643125;    -0.00643125;-0.00214375;0.00214375;0.00643125;    -0.00643125;-0.00214375;0.00214375;0.00643125;    -0.00643125;-0.00214375;0.00214375;0.00643125;    -0.00643125;-0.00214375;0.00214375;0.00643125;    -0.00643125;-0.00214375;0.00214375;0.00643125;    -0.00643125;-0.00214375;0.00214375;0.00643125;    -0.00643125;-0.00214375;0.00214375;0.00643125;    -0.00643125;-0.00214375;0.00214375;0.00643125;    -0.00643125;-0.00214375;0.00214375;0.00643125;    -0.00643125;-0.00214375;0.00214375;0.00643125;    -0.00643125;-0.00214375;0.00214375;0.00643125;    -0.00643125;-0.00214375;0.00214375;0.00643125;    -0.00643125;-0.00214375;0.00214375;0.00643125;    -0.00643125;-0.00214375;0.00214375;0.00643125;    -0.00643125;-0.00214375;0.00214375;0.00643125;];
                    y_tr=[-0.03215625;-0.03215625;-0.03215625;-0.03215625;-0.02786875;-0.02786875;-0.02786875;-0.02786875;-0.02358125;-0.02358125;-0.02358125;-0.02358125;    -0.01929375;-0.01929375;-0.01929375;-0.01929375;-0.01500625;-0.01500625;-0.01500625;-0.01500625;-0.01071875;-0.01071875;-0.01071875;-0.01071875;-0.00643125;-0.00643125;-0.00643125;-0.00643125;-0.00214375;-0.00214375;-0.00214375;-0.00214375;0.00214375;0.00214375;0.00214375;0.00214375;0.00643125;0.00643125;0.00643125;0.00643125;0.01071875;0.01071875;0.01071875;0.01071875;0.01500625;0.01500625;0.01500625;0.01500625;0.01929375;0.01929375;0.01929375;0.01929375;0.02358125;0.02358125;0.02358125;0.02358125;0.02786875;0.02786875;0.02786875;0.02786875;0.03215625;0.03215625;0.03215625;0.03215625];
                end

                % Prepare angular parameters for beamforming
                theta = 0:88:88; % Horizontal angles
                if length(array_angles) > 1
                    phi = [57, 0]; % Vertical angles for array
                else
                    phi = -90:0.025:90; % Fine angular resolution
                end
                [theta_mesh, phi_mesh] = meshgrid(theta, phi); % Generate angle pairs
                theta_m = theta_mesh(:).'; % Flatten angle arrays
                phi_m = phi_mesh(:).';

                % Define filtering and processing parameters
                Pulses = 20; % Number of pulses in the signal
                f_lo = 39000; % Lower frequency bound for bandpass filter
                f_up = 41000; % Upper frequency bound for bandpass filter
                N_filter = floor(Pulses / f_0 * f_s); % Filter length

                % Define time window for processing Lamb waves and the Filter length N
                if eval_ldv
                    N = size(LDV_datas,1);
                    time_window_lamb_wave = round(([0.0009 0.0024]) .* f_s);
                else
                    N = size(array_samples,1);
                    time_window_lamb_wave = round(([0.0015 0.0024]) .* f_s);
                end
                time_window_lamb_wave = time_window_lamb_wave(1):time_window_lamb_wave(2);
                [~, bpFilter_f] = FilterBP(f_lo, f_up, N_filter, N, f_s); % Create bandpass filter

                % Initialize result arrays for current run
                dummy_init_one_number_per_run = zeros(size(LDV_datas, 3), size(LDV_datas, 4), size(LDV_datas, 5));
                lw_amps_all_ch = dummy_init_one_number_per_run;
                lw_idxs_all_ch = dummy_init_one_number_per_run;
                lw_angs = dummy_init_one_number_per_run;
                lw_angs_time_mean = dummy_init_one_number_per_run;
                tt_corr = NaN(size(LDV_datas,3),size(LDV_datas,4), size(LDV_datas,5));
                tt_thres = tt_corr;
                Timestamps_inner =  zeros(size(LDV_datas,4), size(LDV_datas,5)); %same as dummy, but without angles
                Timestamps_std_inner = Timestamps_inner;
                lin_axIdx = 1;%for future use

                % Process raw data for time and amplitude analysis
                for outRepIdx_old = 1:length(outer_Repetition_Names_old)
                    for inRepIdx = 1:length(inner_Repetition_Names)
                        Timestamps_inner(inRepIdx, outRepIdx_old) = mean(meas_timestamps(:, :, inRepIdx, outRepIdx_old), "all");
                        Timestamps_std_inner(inRepIdx, outRepIdx_old) = std(meas_timestamps(:, :, inRepIdx, outRepIdx_old), 0, "all");
						
                        for angleIdx = 1:length(array_angles)
							% Do a progress report
                            meas_nr = (angleIdx-1) + (inRepIdx-1)*length(array_angles) + (outRepIdx_old-1)*length(array_angles)*length(inner_Repetition_Names);
                            inner_meas_length = length(outer_Repetition_Names_old)*length(array_angles)*length(inner_Repetition_Names);
                            meas_length = 9* (inner_meas_length*length(TX_RX_file_extensions)*length(outer_Repetition_Names));
                            disp("File Progress: " + string(meas_nr./(length(outer_Repetition_Names_old)*length(array_angles)*length(inner_Repetition_Names)).*100) +"%")
                            disp("Overall Progress: " + string(100*(outRepIdx*inner_meas_length*length(TX_RX_file_extensions) + TX_RX_idx*inner_meas_length + meas_nr)./meas_length) +"%")
                            
                            %% Process the raw data of each single ultrasonic measurement for time and amplitude. Filter, upsample and envelope extraction first
                            if eval_ldv % use either the LDV signal or the array signal in this loop
                                signals_real=LDV_datas(:,:,angleIdx, inRepIdx,outRepIdx_old)';
                            else
                                signals_real=array_samples(:,:,lin_axIdx,angleIdx, inRepIdx,outRepIdx_old)';
                            end
                            % create analytic signal
                            signals_a_f = fft(signals_real, [], 2); %convert to frequency domain
                            signals_a_f(:, floor(N/2+1):end) = 0; % remove redundant half of frequency domain (analytic spectrum)
                            signals_a_f = 2.*signals_a_f; %multipy by two to account for removed energy of redundant half
                            signals_a_f= signals_a_f.*bpFilter_f;
                            signals_a = ifft(signals_a_f, [], 2);

                            if eval_ldv && TX_RX_file_extensions(TX_RX_idx) == "RX1" % calculate transit-time only with ldv in rx-case - for processing time
                                %upsampling:
                                upsample_factor = 20;
                                signals_a_up = ifft(signals_a_f, N*upsample_factor, 2);
                                if ~exist("signals_a_ldv_ref", "var") % use the first signal as a reference
                                    ref_start_idx = (time_window_lamb_wave(1)-100)*upsample_factor;
                                    ref_end_idx = length(signals_a_up);%(time_window_lamb_wave(end))*upsample_factor;
                                    signals_a_ldv_ref = signals_a_up(ref_start_idx:ref_end_idx);
                                end
                                % Matlab correlation time-shift method.
                                [corr, lags] = xcorr(real(signals_a_ldv_ref), real(signals_a_up));
                                [~, idx] = max(corr);
                                time_shift_idx = lags(idx);
                                tt_corr(angleIdx, inRepIdx, outRepIdx_old) = (-time_shift_idx+1-ref_start_idx)./upsample_factor./f_s;
								
                                % Threshold method for transit time method
                                threshold = 0.2;
                                time_shift_idx_thres = find(real(signals_a_up)>threshold*max(abs(signals_a_up)),1,"first");
                                tt_thres(angleIdx, inRepIdx, outRepIdx_old) = time_shift_idx_thres./upsample_factor./f_s;
                            end
                            if eval_array && TX_RX_file_extensions(TX_RX_idx) == "RX1"% Beamforming only with the array in the RX case
                                upsample_factor_array = 4;
                                signals_a_up = ifft(signals_a_f, N*upsample_factor_array, 2);

                                % do beamforming
                                a_hat = exp(-1i*2.*pi.*f_0./c .* ( x_tr * (sind(theta_m) .* cosd(phi_m)) + y_tr * sind(phi_m))).'; %beamforming vector
                                signals_bf_f = a_hat * signals_a_f; %beamforming
                                signals_bf_a_up = ifft(signals_bf_f, N*upsample_factor_array, 2); %upsamplen hiere
                                signals_bf_env = abs(signals_bf_a_up); %create envelope
                                signals_bf_a = signals_bf_a_up;

                                %do a 2d search in the defined bounds for maximum
                                signals_bf_env_at_thetazero = signals_bf_env(theta_m == 0,:);
                                phi_m_at_thetazero = phi_m(theta_m == 0);
                                time_window_lamb_wave_upsampled = time_window_lamb_wave(1)*upsample_factor_array:time_window_lamb_wave(end)*upsample_factor_array;
                                signals_bf_env_in_2d_bounds = signals_bf_env(theta_m == 0,time_window_lamb_wave_upsampled);
                                [all_maxValues, all_maxIndices] = max(signals_bf_env_in_2d_bounds, [], "all");
                                [angle_idx, time_idx_short] = ind2sub(size(signals_bf_env_in_2d_bounds), all_maxIndices);
                                time_idx = time_idx_short + time_window_lamb_wave(1)*upsample_factor_array-1;
                                best_estimate_for_lw_ang = phi_m_at_thetazero(angle_idx);

                                %save
                                lw_amps_all_ch(angleIdx, inRepIdx, outRepIdx_old) = all_maxValues;
                                lw_idxs_all_ch(angleIdx, inRepIdx, outRepIdx_old) = time_idx/upsample_factor_array;
                                lw_angs(angleIdx, inRepIdx, outRepIdx_old) = best_estimate_for_lw_ang;
                                lw_angs_time_mean(angleIdx, inRepIdx, outRepIdx_old) = best_estimate_for_lw_ang; %obsolet
                            end
                            if eval_array && TX_RX_file_extensions(TX_RX_idx) == "TX" % in the TX evaluation case, just calculate the received amplitude of the phased array
                                upsample_factor_array = 4;
                                signals_a_up = ifft(signals_a_f, N*upsample_factor_array, 2);
                                time_window_lamb_wave_upsampled = time_window_lamb_wave(1)*upsample_factor_array:time_window_lamb_wave(end)*upsample_factor_array;
                                [lamb_wave_amplitude, lw_idx] = max(abs(signals_a_up(:,time_window_lamb_wave_upsampled)), [],2); %take the maximum of each channel seperately
                                lw_amps_all_ch(angleIdx, inRepIdx, outRepIdx_old) = median(lamb_wave_amplitude);
                                lw_idxs_all_ch(angleIdx, inRepIdx, outRepIdx_old) = median(lw_idx)/upsample_factor_array;
                            end
                            if eval_ldv && TX_RX_file_extensions(TX_RX_idx) == "TX"% in the TX evaluation case, just calculate the received amplitude of the LDV
                                [lamb_wave_amplitude, lw_idx] = max(abs(signals_a(:,time_window_lamb_wave)), [],2); %no beamforming
                                lw_amps_all_ch(angleIdx, inRepIdx, outRepIdx_old) = sum(lamb_wave_amplitude);
                                lw_idxs_all_ch(angleIdx, inRepIdx, outRepIdx_old) = mean(lw_idx);
                            end


                            %% Now do some plots about this run
                            t = (1:size(signals_a_f, 2))./f_s;
                            t_up = (1:size(signals_a_f, 2)*10)./f_s/10;

                            if do_plot && (angleIdx==42 && eval_ldv)
                                figure(1)

                                hold off,
                                plot(t_up*1000,real(ifft(signals_a_f, N*10, 2))'/2*49),  hold on,
                                xlabel("Time (ms)")
                                ylabel("Surface velocity (mm/s)")
                                text(0.65, 0.55, '(1)')
                                ylim([-0.8 0.8])
                                xlim([0.5,3.9])
                                fi=gcf;

                                set(fi,'Units','centimeters',"OuterPosition", [2 2 8.8/2+1 8.0]);
                                drawnow
                            end
                            if do_plot && (TX_RX_idx==2 && eval_array)
                                t_up_array = (1:size(signals_a_f, 2)*upsample_factor_array)./f_s/upsample_factor_array;
                                figure(3)
                                hold off,
                                plotsignals_real = signals_real;
                                plotsignals_real(:,1:78) = -2+2.*randn(64,78); %there is always an artifact from the not-used T/R-switch in the array electronic. I suppress it, in order to not distract from the important part of the data.
                                [~,bpFilter_f2] = FilterBP(f_lo, f_up, round(N_filter/4), N, f_s);
                                signals_a_f_2= 2.*fft(plotsignals_real, N, 2).*bpFilter_f2;
                                signals_a_f_2(:, floor(N/2+1):end) = 0; % remove redundant half of frequency domain (analytic spectrum)
                                plotsignals3 = ifft(signals_a_f_2, N*10, 2);
                                plotsignals3(:,1:100) = 0;
                                plot(t_up*1000,plotsignals3(1,:)./ max(plotsignals3(1,:)))
                                xlabel("Time (ms)")
                                ylabel("Acoustic pressure (norm.)") %one could calc this into Pa: samples - voltage - mv/Pa
                                ylim([-1.05 1.05])
                                xlim([0.5,3.9])
                                text(1.1, 0.58, '(1)')
                                text(2.8, 0.85, '(2)')
                                %text(4.7, 0.7, '(3)')
                                fi=gcf;
                                set(fi,'Units','centimeters',"OuterPosition", [2 2 8.8/2+1 8.0]);

                                figure(4)
                                hold off,
                                probe_angle = 54;
                                idx = find(phi_m==probe_angle & theta_m ==0);
                                plot(t_up_array*1000,real(signals_bf_a(idx,:))'),  hold on,
                                plot(repmat((time_window_lamb_wave(1)/f_s)*1000, [1,2]), [-10,10], "xk")
                                plot(repmat((time_window_lamb_wave(end)/f_s)*1000, [1,2]), [-10,10], "xk")
                                title("Beamformed signal at " + string(array_angles(angleIdx))),
                                xlabel("Time (ms)")
                                ylabel("Acoustic pressure (arb. unit)") %one could calc this into Pa: samples - voltage - mv/Pa

                                figure(5)
                                hold off,
                                sample_idx = time_idx;%round( 2.215/1000*f_s);
                                y_asdf = abs(signals_bf_a(theta_m==0,sample_idx));
                                plot(phi_m(theta_m==0),y_asdf./max(y_asdf)),  hold on,
                                xlabel("Incident angle (°)")
                                ylabel("Amplitude (norm.)")
                                axis tight
                                ylim([0,1.02])
                                % Parameters for the theoretical array
                                N_t = 16;                 % Number of elements in the array
                                d_t = 0.5;                % Distance between elements in wavelengths
                                phi_t = phi_m(theta_m == 0) ./ 180 * pi;  % Angle range converted to radians
                                k_t = 2 * pi;            % Wave number (assuming wavelength = 1)
                                % Steering Angle
                                steering_angle_t = best_estimate_for_lw_ang;    % Steering angle in degrees
                                steering_angle_rad = steering_angle_t / 180 * pi;  % Steering angle in radians
                                % Calculate the Array Factor with steering
                                AF = abs(sin(N_t * k_t * d_t * (sin(phi_t) - sin(steering_angle_rad)) / 2) ...
                                    ./ (N_t * sin(k_t * d_t * (sin(phi_t) - sin(steering_angle_rad)) / 2)));
                                % Normalize the Array Factor
                                AF = AF / max(AF);
                                % Plot the directivity pattern
                                plot(phi_m(theta_m == 0), AF);  % Plot against phi_m for theta_m == 0
                                plot([best_estimate_for_lw_ang,best_estimate_for_lw_ang],[0,1],"k--")
                                legend("Measured beam pattern", "Theoretical beam pattern", "Incident angle", "Location", "North west")

                                figure(57), hold off,
                                [phi_mesh_N, N_mesh] = meshgrid(phi, double(1:size(signals_bf_env, 2)));
                                surf(phi_mesh_N, N_mesh/f_s*1000/upsample_factor_array, signals_bf_env_at_thetazero'./max(signals_bf_env_at_thetazero,[],'all'));
                                %title('Envelopes of theta=0 direction');
                                xlabel('Angle (°)');
                                ylabel('Time (ms)');
                                axis tight
                                shading interp;
                                view(0,90);
                                colorbar
                                ylim([0.5, 3.9])%
                                hold on
                                drawnow
                            end
                        end
                    end
                end
                
                %% do another loop to evaluate the TX-angle (requires multiple ultrasonic measurement points)
                if length(array_angles)>1
                    factor_crop = 0.5; %Examination shows that 0.8 with every 1° a data point is basically enough (-63 bis -45) und coup_angle_lw_amps_all_ch ist besser
                    outRepIdx_old= 1;
                    inRepIdx = 1;

                    for factor_crop = 0.25
                        for idx_reduction = 1
                            used_idxs = 1:idx_reduction:length(array_angles);
                            coup_angle_lw_amps_all_ch = zeros(length(inner_Repetition_Names), length(outer_Repetition_Names_old));
                            c_std_all = zeros(length(inner_Repetition_Names), length(outer_Repetition_Names_old));
                            for outRepIdx_old = 1:length(outer_Repetition_Names_old)
                                for inRepIdx = 1:length(inner_Repetition_Names)
                                    when_to_do_the_plot_here = eval_ldv && do_plot;
                                    if when_to_do_the_plot_here
                                        figure(2) %plot the result
                                        hold off
                                    end
                                    [coup_angle_lw_amps_all_ch(inRepIdx,outRepIdx_old), rmse, c_std_all(inRepIdx,outRepIdx_old)] = fit_to_angle_sweep(array_angles(used_idxs),lw_amps_all_ch(used_idxs, inRepIdx, outRepIdx_old), factor_crop, when_to_do_the_plot_here);
                                end
                            end
                        end
                    end
                else
                    coup_angle_lw_amps_all_ch = NaN(length(inner_Repetition_Names), length(outer_Repetition_Names_old));
                    c_std_all = NaN(length(inner_Repetition_Names), length(outer_Repetition_Names_old));

                end
                
                %% stuff all the extracted data into nice little arrays
                if contains( TX_RX_file_extensions(TX_RX_idx) , "TX")
                    TX_coup_angles(1:size(coup_angle_lw_amps_all_ch,1),1:size(coup_angle_lw_amps_all_ch,2),outRepIdx,eval_array_idx,linAx_x_idx, TX_RX_idx) = coup_angle_lw_amps_all_ch;
                    TX_c_std_all(1:size(c_std_all,1),1:size(c_std_all,2),outRepIdx,eval_array_idx,linAx_x_idx, TX_RX_idx) = c_std_all;
                else
                    lw_angs_permute = permute(lw_angs,[2,3,1]);
                    RX_coup_angles(1:size(lw_angs_permute,1),1:size(lw_angs_permute,2),outRepIdx,eval_array_idx,linAx_x_idx, TX_RX_idx) = lw_angs_permute;

                    tt_corr_permute = permute(tt_corr,[2,3,1]);
                    tt_corr_times(1:size(tt_corr_permute,1),1:size(tt_corr_permute,2),outRepIdx,eval_array_idx,linAx_x_idx, TX_RX_idx) = tt_corr_permute;

                    tt_thres_permute = permute(tt_thres,[2,3,1]);
                    tt_thres_times(1:size(tt_thres_permute,1),1:size(tt_thres_permute,2),outRepIdx,eval_array_idx,linAx_x_idx, TX_RX_idx) = tt_thres_permute;


                    lw_idxs_all_ch_permute = permute(lw_idxs_all_ch,[2,3,1]);
                    RX_max_idxs(1:size(lw_idxs_all_ch_permute,1),1:size(lw_idxs_all_ch_permute,2),outRepIdx,eval_array_idx,linAx_x_idx, TX_RX_idx) = lw_idxs_all_ch_permute;
                end
                Temperatures(outRepIdx,eval_array_idx, TX_RX_idx) = mean(temperatures);
                Temperatures_std(outRepIdx,eval_array_idx, TX_RX_idx) = std(temperatures);
                Timestamps(1:size(Timestamps_inner,1),1:size(Timestamps_inner,2),outRepIdx,eval_array_idx,linAx_x_idx, TX_RX_idx) = Timestamps_inner;
                Timestamps_std(1:size(Timestamps_inner,1),1:size(Timestamps_inner,2),outRepIdx,eval_array_idx,linAx_x_idx, TX_RX_idx) = Timestamps_std_inner;
            end
        end
    end

    %% save the extracted data from all of the measurements
    %save(FileName_Base + "IFSWeval.mat", "outer_Repetition_Names", "outer_Repetition_Names_old", "inner_Repetition_Names","linAx_x","eval_array_n", "TX_RX_file_extensions", "TX_coup_angles","TX_c_std_all", "RX_coup_angles", "Timestamps", "Timestamps_std", "Temperatures", "Temperatures_std", "RX_max_idxs", "tt_corr_times", "tt_thres_times")


end