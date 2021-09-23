clear;
allMouse = {'cd017','cd036','cd037','cd042','cd044'};
savePath = 'D:\labData\excitatory\tuning\masterData\';
for currCol = 1:length(allMouse)
    load([savePath '\' allMouse{currCol} '\' 'prePostTuning.mat']);
    disp(['mouse ' allMouse{currCol}]);

    regPreAct = [squeeze(mean(peakActPre(targIdx,:,:),2))';...
        squeeze(mean(peakActPre(foilIdx,:,:),2))'];
    regPreActAll = [squeeze(mean(peakActPre(:,:,:),2))];
    regPrePeak = [preTuning.tuningPeak];
    tuningPeakVec = zeros(17,nNeuron);
    for i = 1:17; tuningPeakVec(i,(preTuning.tuningPeak==i)) = 1; end 
    regPrePeakVec = [tuningPeakVec];
    regPreTone = [responsiveTonePre];

    tuningChange = - preTuning.tuningPeak + postTuning.tuningPeak;
    tuningChangeLog = tuningChange;
    tuningChangeLog(tuningChange>0) = 1;
    tuningChangeLog(tuningChange<0) = -1;
    
    targPre = squeeze(peakActPre(targIdx,:,:));
    targPost = squeeze(peakActPost(targIdx,:,:));
    targDiff = - targPre + targPost;

    preAcc = mean(decoderPre.cellAcc(:,:,2),2) - 0.5;
    postAcc = mean(decoderPost.cellAcc(:,:,2),2) - 0.5;
    decodeDiff = - preAcc + postAcc;


    getRegress({tuningChange, tuningChangeLog,decodeDiff},...
        {regPreAct, regPreActAll, regPrePeak, regPrePeakVec, regPreTone,...
        [regPreAct;regPreActAll;regPrePeak;regPrePeakVec;regPreTone]},...
        {staySel, staySel,decodeFlag},...
        {'tuning change', 'tuning change logical', 'decode diff'},...
        {'pre T/F', 'pre act' , 'pre peak','pre peak vec','pre tone','all'});

end


predictFlag = (signifFlagPre | signifFlagPost);
preAcc = zeros(1,nNeuron);preAcc(signifFlagPre) = cellWPre;
postAcc = zeros(1,nNeuron);postAcc(signifFlagPost) = cellWPost;

function getRegress(data, regressor, dataFlag,dataName,regressName)
    for j = 1:length(data)
        disp(['Fit ' dataName{j}])
        for i = 1:length(regressor)
            results = fitlm(regressor{i}(:,dataFlag{j})',data{j}(dataFlag{j})');
            disp([regressName{i} ': r-square = ' num2str(results.Rsquared.Adjusted) ])
        end
    end
end