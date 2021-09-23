clear;
allMouse = {'cd017','cd036','cd037','cd042','cd044'};
savePath = 'D:\labData\excitatory\tuning\masterData\';

figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.1, 0.95, 0.9]);
nCol = length(allMouse);
nRow = 3;
margins = [0.04, 0.04];

for currCol = 1:length(allMouse)
    load([savePath '\' allMouse{currCol} '\' 'prePostTuning.mat']);

    subplot_tight(nRow,nCol, nCol*0+currCol,margins)
    plot(peakCountPreAll,'LineWidth',3); hold on; plot(peakCountPostAll,'LineWidth',3); 
    ylimm = ylim; xlim([1 17])
    plot([targIdx targIdx],ylimm,'--g');plot([foilIdx foilIdx],ylimm,'--r');
    legend('pre','post');ylim(ylimm); ylabel('counts')
    xticks([1 5 9 13 17]); xticklabels({'4', '8', '16', '32', '64'})
    title(['peak ' allMouse{currCol}])

    subplot_tight(nRow,nCol, nCol*1+currCol,margins)
    plot(sum(responsiveTonePre,2),'LineWidth',3); hold on; plot(sum(responsiveTonePost,2),'LineWidth',3); 
    ylimm = ylim; xlim([1 17]); 
    plot([targIdx targIdx],ylimm,'--g');plot([foilIdx foilIdx],ylimm,'--r');
    legend('pre','post');ylim(ylimm); ylabel('counts')
    xticks([1 5 9 13 17]); xticklabels({'4', '8', '16', '32', '64'})
    title('resp neurons')

    meanPeakActPre = mean(mean(peakActPre(:,:,responsiveFlagPre),2),3);
    meanPeakActPost = mean(mean(peakActPost(:,:,responsiveFlagPost),2),3);
    subplot_tight(nRow,nCol, nCol*2+currCol,margins)
    plot(meanPeakActPre,'LineWidth',3); hold on; plot(meanPeakActPost,'LineWidth',3); 
    ylimm = ylim; xlim([1 17])
    plot([targIdx targIdx],ylimm,'--g');plot([foilIdx foilIdx],ylimm,'--r');
    legend('pre','post');ylim(ylimm); ylabel('avg zscore')
    xticks([1 5 9 13 17]); xticklabels({'4', '8', '16', '32', '64'})
    title('mean zscore')
end
saveas(gcf,[savePath '\allMouse\f1_avgChange.png']);

