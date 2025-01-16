
Folder = "..\data\sim_data\over_stress\";
Folder_out = "..\data\sim_data\";
DirList = dir(fullfile(Folder, '*.mat'));
datas = cell(1,length(DirList));
for idx = 1:length(DirList)
    datas{1,idx}= load(strcat(Folder,string(DirList(idx).name)));
end

f_interp = 0:0.0001:max(datas{1,1}.ResultTable.Mode0_fd_MHzmm);
ph_V_all = zeros(length(f_interp), length(DirList));
gr_V_all = zeros(length(f_interp), length(DirList));
s_all = zeros(1, length(DirList));
for idx = 1:length(DirList)
    ph_V = datas{1,idx}.ResultTable.Mode0_PhaseVelocity_m_s;
    gr_V = datas{1,idx}.ResultTable.Mode0_GroupVelocity_m_s;
    f =datas{1,idx}.ResultTable.Mode0_fd_MHzmm;
    ph_V_all(:,idx) = interp1(f(~isnan(f)),ph_V(~isnan(ph_V)),f_interp,'spline');
    gr_V_all(:,idx) = interp1(f(~isnan(f)),gr_V(~isnan(gr_V)),f_interp,'spline');
    s_all(:,idx) = datas{1,idx}.S.Prestress_Sigma11;
end


c=346;
f_interp_Hz = f_interp./datas{1,idx}.I.Thickness*1000;
f_40kHz_idx = 201;%due to rounding errors this does not work so well: find(f_interp_Hz==round(40.*0.5./0.4));%sheet thickness was wrong;

[s_all_sorted, sort_idx] = sort(s_all./1000000);  % Sort s_all and get sorting indices
ph_V_all_sorted = ph_V_all(f_40kHz_idx, sort_idx);% Reorder ph_V_all(f_40kHz_idx,:) according to the sorted s_all

phase_vel_sim = ph_V_all(f_40kHz_idx,sort_idx);
group_vel_sim = gr_V_all(f_40kHz_idx,sort_idx);
coup_angle_sim = asind(c./ph_V_all(f_40kHz_idx,sort_idx));
% 
save(Folder_out+"Sim_map_stress_vel.mat", "s_all_sorted", "phase_vel_sim", "group_vel_sim", "coup_angle_sim")

figure,
plot(s_all_sorted, gr_V_all(f_40kHz_idx,sort_idx))
xlabel("Stress (MPa)")
ylabel("Group velocity (m/s)")
xlim([0 500])
ylim([807 873])

figure,
plot(s_all_sorted, ph_V_all(f_40kHz_idx,sort_idx))
xlabel("Stress (MPa)")
ylabel("Phase velocity (m/s)")
xlim([0 500])
ylim([432 479])
