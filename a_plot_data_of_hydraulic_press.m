folder  = "data\press_data\";

FileNames(1) = "ACM_mfp_1_1.ASC";
FileNames(2) = "ACM_mfp_1_1  5 mm weg.ASC";
FileNames(3) = "ACM_mfp_1_1  5 mm weg 2 .ASC";
FileNames(4) = "ACM_mfp_1_1  5 mm weg 3 .ASC";

for fileidx1 = 1:length(FileNames)
    filename_here = FileNames(fileidx1);
    datas_mechifsw(fileidx1) = load(folder+"data_of_"+ filename_here+".mat", "time", "zylinderweg","zylinder_kraft","w1_hof","w2_mensa","dms_1","spielzahl","stress_mech_meas", "T_0_mech_meas_string", "T_0_mech_meas_timestamp");
end


figure, hold off, xlabel("t"),hold on
for fileidx1 = 2%1:length(FileNames)
    plot(datas_mechifsw(fileidx1).time(1:1953), datas_mechifsw(fileidx1).dms_1(1:1953).*210E9./1e6./1e6.*1.06)
    ylabel('Stress (MPa)');
    ylim([-13e1+1 0])
    yyaxis right
    plot(datas_mechifsw(fileidx1).time(1:1953), -datas_mechifsw(fileidx1).zylinderweg(1:1953))
    ylabel('Displacement (mm)');
    ylim([-35 0])
end
xlabel('Time (s)');
xlim([1,399])


figure, hold off
for fileidx1 = 2%1:length(FileNames)
    plot(datas_mechifsw(fileidx1).zylinder_kraft(1:1953), datas_mechifsw(fileidx1).dms_1(1:1953).*210E9./1e6./1e6.*1.06)
end
axis tight
ylabel('Stress (MPa)');
xlim([0,10])
ylim([-13e1+1 0])
xlabel('Force (kN)');


%