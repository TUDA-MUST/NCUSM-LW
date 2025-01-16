%creates some plots that need data from different runs of the experiment

%four files for four different runs of the hydraulic press
FileName_Bases(1) = "data\us_data\2024-11-19-IFSW_zwei_vent_top1";
FileName_Bases(2) = "data\us_data\2024-11-19-IFSW_zwei_vent_top2";
FileName_Bases(3) = "data\us_data\2024-11-19-IFSW_zwei_vent_top3";
FileName_Bases(4) = "data\us_data\2024-11-19-IFSW_zwei_vent_top4";

for fileIdx = 1:length(FileName_Bases)
    res(fileIdx) = load(FileName_Bases(fileIdx) + "IFSWeval.mat", "outer_Repetition_Names", "outer_Repetition_Names_old", "inner_Repetition_Names","linAx_x","eval_array_n", "TX_RX_file_extensions", "TX_coup_angles", "RX_coup_angles", "Timestamps", "Timestamps_std", "Temperatures", "Temperatures_std", "RX_max_idxs", "tt_corr_times", "tt_thres_times");
end

linAx_x_idx = 1;
TX_RX_idx = 1;
for eval_array_idx = 1%1:length(res(fileIdx).eval_array_n)
    figure, hold off,
    for fileIdx = 1:length(FileName_Bases)
        coup_angle = res(fileIdx).TX_coup_angles(:,:,:, eval_array_idx,linAx_x_idx, TX_RX_idx);
        timestamp = res(fileIdx).Timestamps(1:size(coup_angle,1),1:size(coup_angle,2),1:size(coup_angle,3), eval_array_idx,linAx_x_idx, TX_RX_idx);
        timestamp_std = res(fileIdx).Timestamps_std(1:size(coup_angle,1),1:size(coup_angle,2),1:size(coup_angle,3),eval_array_idx,linAx_x_idx, TX_RX_idx);
        % temperatures = res(fileIdx).Temperatures(1:size(coup_angle,1),1:size(coup_angle,2),1:size(coup_angle,3), eval_array_idx,linAx_x_idx, TX_RX_idx);
        % temperatures_std = res(fileIdx).Temperatures_std(1:size(coup_angle,1),1:size(coup_angle,2),1:size(coup_angle,3),eval_array_idx,linAx_x_idx, TX_RX_idx);
        errorbar(timestamp(:), coup_angle(:), 0.*coup_angle(:),0.*coup_angle(:), +timestamp_std(:), -timestamp_std(:), "x")
        %errorbar(fileIdx, mean(coup_angle, "all","omitmissing"), std(coup_angle, 0, "all","omitmissing"))
        hold on

    end
    title("TX, eval array is " + string(res(fileIdx).eval_array_n(eval_array_idx)))
    xlabel("Time (ms)")
    ylabel("Coupling angle (°)")
end

linAx_x_idx = 1;
TX_RX_idx = 2;
for eval_array_idx = 2%1:length(res(fileIdx).eval_array_n)
    figure, hold off,
    for fileIdx = 1:length(FileName_Bases)
        coup_angle = res(fileIdx).RX_coup_angles(:,:,:, eval_array_idx,linAx_x_idx, TX_RX_idx);
        timestamp = res(fileIdx).Timestamps(1:size(coup_angle,1),1:size(coup_angle,2),1:size(coup_angle,3), eval_array_idx,linAx_x_idx, TX_RX_idx);
        timestamp_std = res(fileIdx).Timestamps_std(1:size(coup_angle,1),1:size(coup_angle,2),1:size(coup_angle,3),eval_array_idx,linAx_x_idx, TX_RX_idx);
        errorbar(timestamp(:), coup_angle(:), 0.*coup_angle(:),0.*coup_angle(:), +timestamp_std(:), -timestamp_std(:), "x")
        %errorbar(fileIdx, mean(coup_angle, "all","omitmissing"), std(coup_angle, 0, "all","omitmissing"))
        hold on
    end
    title("RX, eval array is " + string(res(fileIdx).eval_array_n(eval_array_idx)))
    xlabel("Time (ms)")
    ylabel("Coupling angle (°)")
end


