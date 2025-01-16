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
%figure,  hold on
for idx = 1:length(DirList)
    ph_V = datas{1,idx}.ResultTable.Mode0_PhaseVelocity_m_s;
    gr_V = datas{1,idx}.ResultTable.Mode0_GroupVelocity_m_s;
    f =datas{1,idx}.ResultTable.Mode0_fd_MHzmm;
    ph_V_all(:,idx) = interp1(f(~isnan(f)),ph_V(~isnan(ph_V)),f_interp,'spline');
    gr_V_all(:,idx) = interp1(f(~isnan(f)),gr_V(~isnan(gr_V)),f_interp,'spline');
    s_all(:,idx) = datas{1,idx}.S.Prestress_Sigma11;
    %plot(f,ph_V)
end


figure,
c=346; %speed of sound
y_lim_left = [-38,39];
y_lim_right = y_lim_left./y_lim_left(2).*16.9;

f_interp_Hz = f_interp./datas{1,2}.I.Thickness*1000.*0.4./0.5;%steel thickness was wrong during simulation data creation
plot(f_interp_Hz, ph_V_all(:,2)-ph_V_all(:,1))
hold on
plot(f_interp_Hz, gr_V_all(:,2)-gr_V_all(:,1))
yyaxis right
coup_angle_diff = -asind(c./ph_V_all(:,2)) + asind(c./ph_V_all(:,1));
plot(f_interp_Hz(f_interp_Hz>24.7), coup_angle_diff(f_interp_Hz>24.7), "Color", TUDcolorsRGB(4,:))
ax = gca;
ax.YAxis(2).Color = TUDcolorsRGB(4,:);
ylim(y_lim_right)
ylabel("Coupling angle difference (°)")
yyaxis left
legend("Phase vel. diff.", "Group vel. diff.", "Coup. ang. diff.")

plot(f_interp_Hz,f_interp_Hz.*0, "k", 'HandleVisibility', 'off')
plot([40,40],y_lim_left, "k--", 'HandleVisibility', 'off')
xlabel("Frequency (kHz)")
ylabel("Phase/group velocity difference (m/s)")
set(gca, 'XScale', 'log')
axis tight
xlim([5.*0.4./0.5,10000.*0.4./0.5])

