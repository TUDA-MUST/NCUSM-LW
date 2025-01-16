
% The data files of the four different runs are evaluated one by one, by
% each uncommenting one of the four following lines and commenting the
% rest.
%filename = 'ACM_mfp_1_1.ASC';T_0_mech_meas_string = "2024-11-19 10:52:33"; T_0_mech_meas_timestamp = milliseconds(datetime(T_0_mech_meas_string, 'Format', 'yyyy-MM-dd HH:mm:ss.SSS')-datetime('2024-3-1'));
filename = 'ACM_mfp_1_1  5 mm weg.ASC';T_0_mech_meas_string = "2024-11-19 10:58:27"; T_0_mech_meas_timestamp = milliseconds(datetime(T_0_mech_meas_string, 'Format', 'yyyy-MM-dd HH:mm:ss.SSS')-datetime('2024-3-1'));
%filename = 'ACM_mfp_1_1  5 mm weg 2 .ASC';T_0_mech_meas_string = "2024-11-19 11:07:32"; T_0_mech_meas_timestamp = milliseconds(datetime(T_0_mech_meas_string, 'Format', 'yyyy-MM-dd HH:mm:ss.SSS')-datetime('2024-3-1'));
%filename = 'ACM_mfp_1_1  5 mm weg 3 .ASC';T_0_mech_meas_string = "2024-11-19 11:22:08"; T_0_mech_meas_timestamp = milliseconds(datetime(T_0_mech_meas_string, 'Format', 'yyyy-MM-dd HH:mm:ss.SSS')-datetime('2024-3-1'));


folder  = "data\press_data\";
filename_here = filename;

% Read the file, skipping the header lines (first 29 lines)
data = readtable(folder+filename_here, 'FileType', 'text', 'Delimiter', '\t', 'HeaderLines', 37);


% The columns of the data array correspond to:
% Column 1: Time (s)
% Column 2: Zylinder-Kraft_150KN_GTM_Dose (kN)
% Column 3: Zylinderweg_600mm (mm)
% Column 4: W1_Hof_IL-300 (mm)
% Column 5: W2_Mensa_IL-300 (mm)
% Column 6: DMS_1 (µm/m)
% Column 7: Spielzahl

% Extract the columns for analysis
time = str2double(strrep(data{:,1}, ',', '.'));                % Time in seconds
zylinder_kraft = str2double(strrep(data{:,2}, ',', '.'));  % Convert Force
zylinderweg = str2double(strrep(data{:,3}, ',', '.'));      % Convert extension
w1_hof = str2double(strrep(data{:,4}, ',', '.'));  % Convert W1_Hof_IL-300 (mm)         
w2_mensa = str2double(strrep(data{:,5}, ',', '.'));    % Convert W2_Mensa_IL-300 (mm)     
dms_1 = str2double(strrep(data{:,6}, ',', '.'));     % Convert DMS_1 (µm/m)       
spielzahl = str2double(strrep(data{:,7}, ',', '.'));   % Convert Spielzahl  

% Plot the data
figure;
subplot(3,2,1);
plot(time, zylinder_kraft);
title('Zylinder-Kraft (kN)');
xlabel('Time (s)');
ylabel('Force (kN)');

subplot(3,2,2);
plot(time, zylinderweg);
title('Zylinderweg (mm)');
xlabel('Time (s)');
ylabel('Displacement (mm)');

subplot(3,2,3);
plot(time, w1_hof);
title('W1 Hof IL-300 (mm)');
xlabel('Time (s)');
ylabel('Displacement (mm)');

subplot(3,2,4);
plot(time, w2_mensa);
title('W2 Mensa IL-300 (mm)');
xlabel('Time (s)');
ylabel('Displacement (mm)');

subplot(3,2,5);
plot(time, dms_1);
title('DMS 1 (\mum/m)');
xlabel('Time (s)');
ylabel('Strain (\mum/m)');

subplot(3,2,6);
plot(time, spielzahl);
title('Spielzahl');
xlabel('Time (s)');
ylabel('Count');



figure, 
plot(zylinder_kraft, dms_1)

ylim([0,890])
ylabel('Strain (\mum/m)');
yyaxis right
xlim([0,14])
plot(zylinder_kraft, zylinderweg)
ylim([0,45])
ylabel('Displacement (mm)');
xlabel('Force (kN)');


% save in a more compact format (was already done)
%save(folder+"data_of_"+ filename_here+".mat", "time", "zylinderweg","zylinder_kraft","w1_hof","w2_mensa","dms_1","spielzahl","stress_mech_meas", "T_0_mech_meas_string", "T_0_mech_meas_timestamp")
