addpath("..\functions\")
defineTUDcolors

Folder = "..\data\sim_data\changes\";
DirList = dir(fullfile(Folder, '*.mat'));
datas = cell(1,length(DirList));
for idx = 1:length(DirList)
    datas{1,idx}= load(strcat(Folder,string(DirList(idx).name)));
end



f_interp = 0:0.0001:max(datas{1,1}.ResultTable.Mode0_fd_MHzmm);
ph_V_all = zeros(length(f_interp), length(DirList));
gr_V_all = zeros(length(f_interp), length(DirList));
s_all = zeros(1, length(DirList));
idx = 1;
ph_V0 = datas{1,idx}.ResultTable.Mode0_PhaseVelocity_m_s;
ph_V1 = datas{1,idx}.ResultTable.Mode1_PhaseVelocity_m_s;
ph_V2 = datas{1,idx}.ResultTable.Mode2_PhaseVelocity_m_s;
ph_V3 = datas{1,idx}.ResultTable.Mode3_PhaseVelocity_m_s;
ph_V4 = datas{1,idx}.ResultTable.Mode4_PhaseVelocity_m_s;
ph_V5 = datas{1,idx}.ResultTable.Mode5_PhaseVelocity_m_s;
ph_V6 = datas{1,idx}.ResultTable.Mode6_PhaseVelocity_m_s;
f =datas{1,idx}.ResultTable.Mode0_fd_MHzmm;
f_Hz0 = datas{1,idx}.ResultTable.Mode0_fd_MHzmm./datas{1,idx}.I.Thickness*1000.*0.4./0.5;%steel thickness was wrong during simulation data creation
f_Hz1 = datas{1,idx}.ResultTable.Mode1_fd_MHzmm./datas{1,idx}.I.Thickness*1000.*0.4./0.5;
f_Hz2 = datas{1,idx}.ResultTable.Mode2_fd_MHzmm./datas{1,idx}.I.Thickness*1000.*0.4./0.5;
f_Hz3 = datas{1,idx}.ResultTable.Mode3_fd_MHzmm./datas{1,idx}.I.Thickness*1000.*0.4./0.5;
f_Hz4 = datas{1,idx}.ResultTable.Mode4_fd_MHzmm./datas{1,idx}.I.Thickness*1000.*0.4./0.5;
f_Hz5 = datas{1,idx}.ResultTable.Mode5_fd_MHzmm./datas{1,idx}.I.Thickness*1000.*0.4./0.5;
f_Hz6 = datas{1,idx}.ResultTable.Mode6_fd_MHzmm./datas{1,idx}.I.Thickness*1000.*0.4./0.5;
s_all(:,idx) = datas{1,idx}.S.Prestress_Sigma11;

f_interp_Hz = f_interp./datas{1,idx}.I.Thickness*1000;

figure,  hold on
plot(f_Hz0, ph_V0, "Color", TUDcolorsRGB(1,:))
hold on
plot(f_Hz2, ph_V2, "Color", TUDcolorsRGB(2,:))
legend("Asymmetric","Symmetric", "location","north west")
plot(f_Hz4, ph_V4,"--", "Color", TUDcolorsRGB(1,:), 'HandleVisibility', 'off')
plot(f_Hz5, ph_V5,"--", "Color", TUDcolorsRGB(2,:), 'HandleVisibility', 'off')
xlabel("Frequency (kHz)")

plot([40,40], [0,10000],"k--",'HandleVisibility', 'off')
ylabel("Phase velocity (m/s)")

set(gca, 'XScale', 'log')
xlim([5.*0.4./0.5,10000.*0.4./0.5])
xticks([10 40 100 1000])
xticklabels(["10^1", "4\cdot10^1", "10^2", "10^3"])

ylim([0 8000])
yticks([0 440 2000 4000 5407 6000 8000])






