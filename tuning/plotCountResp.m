clear;
allMouse = {'cd017','cd036','cd037','cd042','cd044'};
%allMouse = {'cd017'};
savePath = 'D:\labData\excitatory\tuning\masterData\';

figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.1, 0.95, 0.9]);
nCol = length(allMouse);
nRow = 4;
margins = [0.04, 0.04];

for currCol = 1:length(allMouse)
    load([savePath '\' allMouse{currCol} '\' 'prePostTuning.mat']);

    tPre = (targIdx == preTuning.tuningPeak) & responsiveFlagPre;
    fPre = (foilIdx == preTuning.tuningPeak) & responsiveFlagPre;
    tPost = (targIdx == postTuning.tuningPeak) & responsiveFlagPost;
    fPost = (foilIdx == postTuning.tuningPeak) & responsiveFlagPost;
    targResp{1} = (responsiveTonePre(targIdx,:) & responsiveFlagPre);
    targResp{2} = (responsiveTonePost(targIdx,:) & responsiveFlagPost );
    foilResp{1} = (responsiveTonePre(foilIdx,:) & responsiveFlagPre);
    foilResp{2} = (responsiveTonePost(foilIdx,:) & responsiveFlagPost );
    tPreZ = mean(mean(peakActPre(targIdx,:,responsiveFlagPre),2),3);
    tPostZ = mean(mean(peakActPost(targIdx,:,responsiveFlagPost),2),3);
    fPreZ = mean(mean(peakActPre(foilIdx,:,responsiveFlagPre),2),3);
    fPostZ = mean(mean(peakActPost(foilIdx,:,responsiveFlagPost),2),3);
    
    tempData = {[sum(responsiveFlagPre) sum(responsiveFlagPost)],...
        [sum(tPre) sum(tPost); sum(fPre) sum(fPost)], ...
    [sum(targResp{1}) sum(targResp{2}); sum(foilResp{1}) sum(foilResp{2})],...
    [tPreZ, tPostZ; fPreZ, fPostZ] };
    h = {subplot_tight(nRow,nCol, nCol*0+currCol,margins),...
        subplot_tight(nRow,nCol, nCol*1+currCol,margins),...
        subplot_tight(nRow,nCol, nCol*2+currCol,margins),...
        subplot_tight(nRow,nCol, nCol*3+currCol,margins)};
    func_barPlot(tempData, 'axes', h, 'legend',...
       {'tone responsive',{'pre','post'},{'pre','post'},{'pre','post'}}, 'title',...
       {['all responsive neuron ' allMouse{currCol}],'peak','t/f responsive','zscore'},...
       'xticklabels',{{'pre','post'},{'target','foil'},{'target','foil'},{'target','foil'}});
    
end
saveas(gcf,[savePath '\allMouse\f2_avgChange.png']);

