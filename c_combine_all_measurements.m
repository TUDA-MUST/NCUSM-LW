addpath("functions\")

FileName_Bases(1) = "data\us_data\2024-11-19-IFSW_zwei_vent_top1";
FileName_Bases(2) = "data\us_data\2024-11-19-IFSW_zwei_vent_top2";
FileName_Bases(3) = "data\us_data\2024-11-19-IFSW_zwei_vent_top3";
FileName_Bases(4) = "data\us_data\2024-11-19-IFSW_zwei_vent_top4";

for fileIdx = 1:length(FileName_Bases)
    res(fileIdx) = load(FileName_Bases(fileIdx) + "IFSWeval.mat", "outer_Repetition_Names", "outer_Repetition_Names_old", "inner_Repetition_Names","linAx_x","eval_array_n", "TX_RX_file_extensions", "TX_coup_angles","TX_c_std_all", "RX_coup_angles", "Timestamps", "Timestamps_std", "Temperatures", "Temperatures_std", "RX_max_idxs", "tt_corr_times", "tt_thres_times");
end

mech_folder  = "data\press_data\";

FileName_Bases_m(1) = "ACM_mfp_1_1.ASC";
FileName_Bases_m(2) = "ACM_mfp_1_1  5 mm weg.ASC";
FileName_Bases_m(3) = "ACM_mfp_1_1  5 mm weg 2 .ASC";
FileName_Bases_m(4) = "ACM_mfp_1_1  5 mm weg 3 .ASC";

for fileIdx = 1:length(FileName_Bases_m)
    res_mech(fileIdx) = load(mech_folder+"data_of_"+ FileName_Bases_m(fileIdx)+".mat", "time", "zylinderweg","zylinder_kraft","w1_hof","w2_mensa","dms_1","spielzahl","stress_mech_meas", "T_0_mech_meas_string", "T_0_mech_meas_timestamp");
    timestamp_mech{fileIdx} = res_mech(fileIdx).T_0_mech_meas_timestamp + res_mech(fileIdx).time*1000;
end

sim_folder = "data\sim_data\";
sim_filename = "Sim_map_stress_vel.mat";
load(sim_folder+ sim_filename , "s_all_sorted", "phase_vel_sim", "group_vel_sim", "coup_angle_sim")


%the timestamps for the mechanical measurement are less. so for every timestamp with a mechanical measurement there also exists an ultrasonic measurement close by. but not  vice versa.
%But there are far less ultrasonic measurements than mechanical measurements!
%common time limits is mechanical timestamp- limits, and the common time points  are from the us-data

stresses_sim = [flipud(-s_all_sorted')' s_all_sorted];
angles_sim_TX = -[flipud(-coup_angle_sim')'+2.*coup_angle_sim(1) coup_angle_sim];
angles_sim_Rx = [flipud(-coup_angle_sim')'+2.*coup_angle_sim(1) coup_angle_sim];
group_vel_sim_tt = [flipud(group_vel_sim')' -group_vel_sim+2.*group_vel_sim(1)];

MPa_per_degree = 100/(angles_sim_Rx(stresses_sim==-100)-angles_sim_Rx(stresses_sim==0.1));
MPa_per_mps_group_velocity = 100/(group_vel_sim_tt(stresses_sim==-100)-group_vel_sim_tt(stresses_sim==0.1));

% TX coupling angles
plot_coupling_data_with_mech_singlemeas(res, FileName_Bases,1, 1, 1, "TX_coup_angles", "TX-coupling angle \alpha_{c,tx} (°)", timestamp_mech, res_mech);
plot(stresses_sim, angles_sim_TX)
xlim([-122 0])
legend("Measurement", "Simulation", "Location", "north west")

% RX coupling angles
plot_coupling_data_with_mech_simple(res, FileName_Bases,2, 1, 2, "RX_coup_angles", "RX-coupling angle \alpha_{c,rx} (°)", timestamp_mech, res_mech);
plot(stresses_sim, angles_sim_Rx)
xlim([-122 0])
legend("Measurement", "Simulation", "Location", "north east")

% Transit time correlation
plot_coupling_data_with_mech_simple(res, FileName_Bases, 1,1, 2, "tt_corr_times", "Measured group velocity (m/s)", timestamp_mech, res_mech);
plot(stresses_sim,group_vel_sim_tt)
xlim([-122 0])
legend("Measurement", "Simulation", "Location", "north west")

%% Metrics: Repeatability 

[x_out_tx, y_out_tx, y_out_std_tx] = plot_coupling_data_with_mech_singlemeas(res, FileName_Bases,1, 1, 1, "TX_coup_angles", "TX-coupling angle (°)", timestamp_mech, res_mech);
[x_out_rx, y_out_rx, y_out_std_rx] = plot_coupling_data_with_mech_singlemeas(res, FileName_Bases,2, 1, 2, "RX_coup_angles", "RX-coupling angle (°)", timestamp_mech, res_mech);
[x_out_tt, y_out_tt, y_out_std_tt] = plot_coupling_data_with_mech_singlemeas(res, FileName_Bases, 1,1, 2, "tt_corr_times", "Measured Group velocity (m/s)", timestamp_mech, res_mech);

[~,~, mean_residual_all_tx, std_residual_all_tx, ~,~, ~] = calculate_metrics(x_out_tx, y_out_tx,y_out_std_tx, 'TX');
[~,~, mean_residual_all_rx, std_residual_all_rx, ~,~, ~] = calculate_metrics(x_out_rx, y_out_rx, y_out_std_rx, 'RX');
[~,~, mean_residual_all_tt, std_residual_all_tt, ~,~, ~] = calculate_metrics(x_out_tt, y_out_tt, y_out_std_tt, 'TT');

data_repeatability = [mean_residual_all_tx.*MPa_per_degree mean_residual_all_rx.*MPa_per_degree mean_residual_all_tt.*abs(MPa_per_mps_group_velocity)];
errors_repeatability = [std_residual_all_tx*MPa_per_degree std_residual_all_rx.*MPa_per_degree std_residual_all_tt.*abs(MPa_per_mps_group_velocity)];

% Create the grouped bar plot
figure;
b = bar(data_repeatability);
ylabel('Repeatability (MPa)');
xticks(1:3); % Set x-tick marks for groups
xticklabels({'TX coup. angle', 'RX coup. angle', 'Transit-time'}); % Label the groups
grid on; % Add grid for better readability

%% Metrics: Linearity 

[x_out_tx, y_out_tx, y_out_std_tx] = plot_coupling_data_with_mech_simple(res, FileName_Bases,1, 1, 1, "TX_coup_angles", "TX-coupling angle (°)", timestamp_mech, res_mech);
[x_out_rx, y_out_rx, y_out_std_rx] = plot_coupling_data_with_mech_simple(res, FileName_Bases,2, 1, 2, "RX_coup_angles", "RX-coupling angle (°)", timestamp_mech, res_mech);
[x_out_tt, y_out_tt, y_out_std_tt] = plot_coupling_data_with_mech_simple(res, FileName_Bases, 1,1, 2, "tt_corr_times", "Measured Group velocity (m/s)", timestamp_mech, res_mech);

[~, ~, ~, ~, ~, ~, max_residual_all_tx] = calculate_metrics(x_out_tx, y_out_tx,y_out_std_tx, 'TX');
[~, ~, ~, ~, ~, ~, max_residual_all_rx] = calculate_metrics(x_out_rx, y_out_rx, y_out_std_rx, 'RX');
[~, ~, ~, ~, ~, ~, max_residual_all_tt] = calculate_metrics(x_out_tt, y_out_tt, y_out_std_tt, 'TT');

data_linearity = [max_residual_all_tx.*MPa_per_degree max_residual_all_rx.*MPa_per_degree max_residual_all_tt.*abs(MPa_per_mps_group_velocity)];
% Create the grouped bar plot
figure;
b = bar(data_linearity);
ylabel('Linearity (MPa)');
xticks(1:3); % Set x-tick marks for groups
xticklabels({'TX coup. angle', 'RX coup. angle', 'Transit-time'}); % Label the groups
grid on; % Add grid for better readability




