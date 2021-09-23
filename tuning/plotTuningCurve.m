clear;

allMouse = 'cd044';
savePath = 'D:\labData\excitatory\tuning\masterData\';
load([savePath '\' allMouse '\' 'prePostTuning.mat']);

%%
margins = [0.001;0.02];

meanPeakActPre = squeeze(mean(peakActPre(:,:,responsiveFlagPre),2));
meanPeakActPost = squeeze(mean(peakActPost(:,:,responsiveFlagPost),2));

[~,preidx] = sort(preTuning.tuningPeak(responsiveFlagPre));
[~,postidx] = sort(postTuning.tuningPeak(responsiveFlagPost));

tempMax = prctile([meanPeakActPre(:);meanPeakActPost(:)],98);
tempMin = prctile([meanPeakActPre(:);meanPeakActPost(:)],2);

plotTuningCurveChange({meanPeakActPre,meanPeakActPost},{preidx,postidx},[tempMin tempMax])


%%
meanPeakActPre = squeeze(mean(peakActPre(:,:,respCellFlag),2));
meanPeakActPost = squeeze(mean(peakActPost(:,:,respCellFlag),2));

temp = preTuning.tuningPeak;
temp((~responsiveFlagPre) & responsiveFlagPost) = 18;
temp = temp(respCellFlag);
[~,tempidx] = sort(temp);

plotTuningCurveChangeSplit({meanPeakActPre,meanPeakActPost},{tempidx,tempidx},[tempMin tempMax], [targIdx, foilIdx])


%%
meanPeakActPre = squeeze(mean(peakActPre(:,:,tPre | fPre | middlePre),2));
meanPeakActPost = squeeze(mean(peakActPost(:,:,tPre | fPre | middlePre),2));

[~,preidx] = sort(preTuning.tuningPeak(tPre | fPre | middlePre));
plotTuningCurveChangeSplit({meanPeakActPre,meanPeakActPost},{preidx,preidx},[tempMin tempMax], [targIdx, foilIdx])

figure; histogram(-tuningChange(tPre | fPre | middlePre& responsiveFlagPost))
%%
meanPeakActPre = squeeze(mean(peakActPre(:,:,tSidePre),2));
meanPeakActPost = squeeze(mean(peakActPost(:,:,tSidePre),2));

[~,preidx] = sort(preTuning.tuningPeak(tSidePre));

plotTuningCurveChangeSplit({meanPeakActPre,meanPeakActPost},{preidx,preidx},[tempMin tempMax], [targIdx, foilIdx])

figure; histogram(-tuningChange(tSidePre& responsiveFlagPost))

%%
meanPeakActPre = squeeze(mean(peakActPre(:,:,fSidePre),2));
meanPeakActPost = squeeze(mean(peakActPost(:,:,fSidePre),2));

[~,preidx] = sort(preTuning.tuningPeak(fSidePre));

plotTuningCurveChangeSplit({meanPeakActPre,meanPeakActPost},{preidx,preidx},[tempMin tempMax], [targIdx, foilIdx])

figure; histogram(-tuningChange(fSidePre& responsiveFlagPost))
%%

allMouse = {'cd017','cd036','cd037','cd042','cd044'};
savePath = 'D:\labData\excitatory\tuning\masterData\';

meanMid = [];
meanLeft = []; 
meanRight = [];
for i = 1:length(allMouse)
    load([savePath '\' allMouse{i} '\' 'prePostTuning.mat']);
    
    meanMid = [meanMid -tuningChange(tPre | fPre | middlePre& responsiveFlagPost)];
    if targIdx < foilIdx
        meanLeft = [meanLeft -tuningChange(tSidePre& responsiveFlagPost)];
        meanRight = [meanRight -tuningChange(fSidePre& responsiveFlagPost)];
    else 
        meanLeft = [meanLeft -tuningChange(fSidePre& responsiveFlagPost)];
        meanRight = [meanRight -tuningChange(tSidePre& responsiveFlagPost)];
    end

    
end
margins = [0.2, 0.05];
figure; subplot_tight(1,3,1,margins);histogram(meanLeft); subplot_tight(1,3,2,margins);histogram(meanMid)
subplot_tight(1,3,3,margins);histogram(meanRight)

%%
function plotTuningCurveChange(prepost,prepostidx,tempMinMax)
    margins = [0.001;0.02];
    figure; subplot_tight(1,2,1,margins); imagesc(prepost{1}(:,prepostidx{1})'); caxis(tempMinMax); axis off
    subplot_tight(1,2,2,margins); imagesc(prepost{2}(:,prepostidx{2})'); caxis(tempMinMax); axis off

end


function plotTuningCurveChangeSplit(prepost,prepostidx,tempMinMax,targfoilidx)
    if targfoilidx(2) < targfoilidx(1); targfoilidx = targfoilidx([2 1]); end
    margins = [0.01;0.02];
    marginsBet = 0.02;
    tempLen = [targfoilidx(1)-1, 4, 17 - targfoilidx(2)];
    plotLenSum = 0.5 - margins(1)*4-marginsBet;
    plotLen = plotLenSum ./ 17 * tempLen;
    
    f = figure;
    allStart = [0 0.5+marginsBet];
    for i = 1:2
        start = allStart(i);
        h1 = subplot(1,6,1+3*(i-1)); 
        h2 = subplot(1,6,2+3*(i-1)); 
        h3 = subplot(1,6,3+3*(i-1));
        set(h1, 'Position', [start + margins(1), margins(2), plotLen(1), 1-margins(2)*2])
        set(h2, 'Position', [start + margins(1)*2 + plotLen(1), margins(2), plotLen(2), 1-margins(2)*2])
        set(h3, 'Position', [start + margins(1)*3 + plotLen(1) + plotLen(2), margins(2), plotLen(3), 1-margins(2)*2])
        set(f,'CurrentAxes',h1)
        imagesc(h1,prepost{i}(1:targfoilidx(1)-1,prepostidx{i})'); caxis(tempMinMax); axis off
        set(f,'CurrentAxes',h2)
        imagesc(h2,prepost{i}(targfoilidx(1):targfoilidx(2),prepostidx{i})'); caxis(tempMinMax); axis off
        set(f,'CurrentAxes',h3)
        imagesc(h3,prepost{i}(targfoilidx(2)+1:end,prepostidx{i})'); caxis(tempMinMax); axis off
    end
   

end
