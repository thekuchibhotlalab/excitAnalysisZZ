%% load the data
clear;
experimental = {'cd017','cd036','cd037','cd042','cd044'};
for i = 1:length(experimental)
    expDatapath{i} = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' experimental{i} '\corrected_plane0'];
    tempData = load([expDatapath{i} '\' experimental{i} '_baseline_data_pairwiseCorr.mat']);
    expData{i} = tempData;
    tempData = load([expDatapath{i} '\' experimental{i} '_behaviorSI.mat']);
    expDataSI{i} = tempData;

    % find what days these baseline sessions correspond to 
    configPath = 'C:\Users\zzhu34\Documents\tempdata\excitData\config\mouse';
    filename = [configPath '\' experimental{i} '_config.csv'];
    configTable = readtable(filename);
    baselineSessionFlag = strcmp(configTable.BehavType,'Baseline');
    expData{i}.baselineDay = configTable.Day(baselineSessionFlag);
    % find the day that is common between behavior and baseline, and the
    % corresponding index in 
    [expData{i}.commonDay,expData{i}.commonDayIdx,~] = intersect(expData{i}.baselineDay,expDataSI{i}.daysList ); 
end

ctrl = {'cd019','cd041'};
for i = 1:length(ctrl)
    ctrlDatapath{i} = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' ctrl{i} '\corrected_plane0'];
    tempData = load([ctrlDatapath{i} '\' ctrl{i} '_baseline_data_pairwiseCorr.mat']);
    ctrlData{i} = tempData;
    tempData = load([ctrlDatapath{i} '\' ctrl{i} '_behaviorSI.mat']);
    ctrlDataSI{i} = tempData;

    % find what days these baseline sessions correspond to 
    configPath = 'C:\Users\zzhu34\Documents\tempdata\excitData\config\mouse';
    filename = [configPath '\' ctrl{i} '_config.csv'];
    configTable = readtable(filename);
    baselineSessionFlag = strcmp(configTable.BehavType,'Baseline');
    ctrlData{i}.baselineDay = configTable.Day(baselineSessionFlag);
    if strcmp(ctrl{i},'cd019') % for cd019, behavior start on day 2, so use -1 to align with behavior session
        ctrlData{i}.baselineDay = ctrlData{i}.baselineDay - 1;
    end
    
    % find the day that is common between behavior and baseline, and the
    % corresponding index in 
    [ctrlData{i}.commonDay,ctrlData{i}.commonDayIdx,~] = intersect(ctrlData{i}.baselineDay,expDataSI{i}.daysList );
    
end
allMouse = [experimental ctrl];
%% find the common days between all the animals (control and experimental)
tempDays = 0:25;
for i = 1:length(allMouse)
    if i <= length(experimental)
        tempDays = intersect(tempDays, expData{i}.commonDay);
    else
        tempDays = intersect(tempDays, ctrlData{i -  length(experimental)}.commonDay);
    end
end
commonDays = tempDays; % note: the starting day of each animal is 1, so commondays can be used directly as index. if any animal does not start at day==1 then should not use as index.
%% figure 1 -- get significant matrix in the correlation
zscoreThreshold = 4;
[covMatSignif, cellPairIdx, cellPairCount] = getSignificantMatrix (experimental, expData,zscoreThreshold);
% now do the same with control animals
[covMatSignif_ctrl, cellPairIdx_ctrl, cellPairCount_ctrl] = ...
    getSignificantMatrix (ctrl, ctrlData,zscoreThreshold);

%% figure 1.1 -- plot the significance matrix
% experimental animals
for mouse = 1:length(experimental)
    [value, sortIdx] =sort(nanmean(cellPairCount{mouse},2),'descend');
    figure; 
    for day = 1:length(expData{mouse}.covMat)
        subplot_tight(4,5,day);
        fullMat = covMatSignif{mouse}(:,:,day) + covMatSignif{mouse}(:,:,day)';
        fullMat = fullMat + diag(ones(size(fullMat,1)));
        imagesc(fullMat(sortIdx,sortIdx)); axis off;
    end
    suptitle(experimental{mouse})   
end
% control animals 
for mouse = 1:length(ctrl)
    [value, sortIdx] =sort(nanmean(cellPairCount_ctrl{mouse},2),'descend');
    figure; 
    for day = 1:length(ctrlData{mouse}.covMat)
        subplot_tight(4,5,day);
        fullMat = covMatSignif_ctrl{mouse}(:,:,day) + covMatSignif_ctrl{mouse}(:,:,day)';
        fullMat = fullMat + diag(ones(size(fullMat,1)));
        imagesc(fullMat(sortIdx,sortIdx)); axis off;
    end
    suptitle(ctrl{mouse})   
end
%% figure 1.2 -- plot the correlation matric 
expDataCorrected = removeBadSession(experimental, expData);
for mouse = 1:length(experimental)
    [value, sortIdx] =sort(nanmean(cellPairCount{mouse},2),'descend');
    figure; 
    for day = 1:length(expDataCorrected{mouse}.covMat)
        subplot_tight(4,5,day);
        try
            corrMat = corrcov(expDataCorrected{mouse}.covMat{day});
        catch
            corrMat = expDataCorrected{mouse}.covMat{day};
        end
        imagesc(corrMat(sortIdx,sortIdx)); axis off; caxis([0 0.2])
    end
    suptitle(experimental{mouse})   
end
for mouse = 1:length(ctrl)
    [value, sortIdx] =sort(nanmean(cellPairCount_ctrl{mouse},2),'descend');
    figure; 
    for day = 1:length(ctrlData{mouse}.covMat)
        subplot_tight(4,5,day);
        corrMat = corrcov(ctrlData{mouse}.covMat{day});
        imagesc(corrMat(sortIdx,sortIdx)); axis off; caxis([0 0.2])
    end
    suptitle(ctrl{mouse})   
end
%% figure 1.3 -- plot the significance count disctribution across days
% subplot 1 - mean significant correlation per neuron
figure; subplot(2,2,1);hold on;
tempMean = [];  tempMean_ctrl = [];
for mouse = 1:length(experimental)
    hexp{mouse} = cdfplot(nanmean(cellPairCount{mouse}(:,commonDays),2) / size(cellPairCount{mouse},1));
    set( hexp{mouse}, 'Color', matlabColors(1,0.6), 'LineWidth', 1.5);
    tempMean(:,mouse) = nanmean(cellPairCount{mouse}(:,commonDays),1) / size(cellPairCount{mouse},1);
end
for mouse = 1:length(ctrl)
    hctrl{mouse} = cdfplot(nanmean(cellPairCount_ctrl{mouse}(:,commonDays),2) / size(cellPairCount_ctrl{mouse},1));
    set( hctrl{mouse}, 'Color', matlabColors(2,0.6), 'LineWidth', 1.5);
    tempMean_ctrl(:,mouse) = nanmean(cellPairCount_ctrl{mouse}(:,commonDays),1) / size(cellPairCount_ctrl{mouse},1);
end
%legend([experimental ctrl])
legend([hexp{1} hctrl{1}],'experimental','control','Location','Best'); title('perc of significant pair, average of all days')
xlabel('percent of significant correlation'); ylabel('neuron'); xlim([0 0.6])
% subplot 2 - percent of significant correlation over learning
subplot(4,2,2); hold on;
fn_plotMeanSampleLine(commonDays',tempMean',{'Color',matlabColors(1),'LineWidth',1.5},...
    {'Color',matlabColors(1,0.3),'LineWidth',0.5});
fn_plotMeanSampleLine(commonDays',tempMean_ctrl',{'Color',matlabColors(2),'LineWidth',1.5},...
    {'Color',matlabColors(2,0.3),'LineWidth',0.5});
xlim([commonDays(1) commonDays(end)]); xticklabels([]); ylim([0 0.3]);ylabel('perc pairs')
title('mean perc of significant pair, across days')
subplot(4,2,4); hold on;
[f_line_exp, ~] =fn_plotMeanErrorbar(commonDays',tempMean',matlabColors(1),...
    {'Color',matlabColors(1),'LineWidth',1.5},{'faceAlpha',0.3,'LineStyle','none'});
[f_line_ctrl, ~] = fn_plotMeanErrorbar(commonDays',tempMean_ctrl',matlabColors(2),...
    {'Color',matlabColors(2),'LineWidth',1.5},{'faceAlpha',0.3,'LineStyle','none'});
xlim([commonDays(1) commonDays(end)]); ylim([0 0.3])
xlabel('days'); ylabel('perc pairs'); 
legend([f_line_exp f_line_ctrl],'experimental','control','Location','Best');
% subplot 3 - first day of learning
subplot(2,2,3);hold on; selectDay = commonDays(1);
for mouse = 1:length(experimental)
    hexp{mouse} = cdfplot(cellPairCount{mouse}(:,selectDay) / size(cellPairCount{mouse},1));
    set( hexp{mouse}, 'Color', matlabColors(1,0.6), 'LineWidth', 1.5);
end
for mouse = 1:length(ctrl)
    hctrl{mouse} = cdfplot(cellPairCount_ctrl{mouse}(:,selectDay) / size(cellPairCount_ctrl{mouse},1));
    set( hctrl{mouse}, 'Color', matlabColors(2,0.6), 'LineWidth', 1.5);
end
legend([hexp{1} hctrl{1}],'experimental','control','Location','Best'); 
xlabel('perc of significant correlation pair'); ylabel('# of neuron'); xlim([0 0.6])
title(['perc of significant pair, day' int2str(selectDay)])
% subplot 4 - last day of learning
subplot(2,2,4);hold on; selectDay = commonDays(end-7);
for mouse = 1:length(experimental)
    hexp{mouse} = cdfplot(cellPairCount{mouse}(:,selectDay) / size(cellPairCount{mouse},1));
    set( hexp{mouse}, 'Color', matlabColors(1,0.6), 'LineWidth', 1.5);
end
for mouse = 1:length(ctrl)
    hctrl{mouse} = cdfplot(cellPairCount_ctrl{mouse}(:,selectDay) / size(cellPairCount_ctrl{mouse},1));
    set( hctrl{mouse}, 'Color', matlabColors(2,0.6), 'LineWidth', 1.5);
end
legend([hexp{1} hctrl{1}],'experimental','control','Location','Best'); 
xlabel('perc of significant correlation pair'); ylabel('# of neuron'); xlim([0 0.6])
title(['perc of significant pair, day' int2str(selectDay)])

%% figure 1.4 -- plot the significance pair count, all neurons all days
figure; 
for mouse = 1:length(experimental)
    subplot_tight(1,length(experimental),mouse)
    imagesc(cellPairCount{mouse}); axis off;     
    caxis([prctile(cellPairCount{mouse}(:),5) prctile(cellPairCount{mouse}(:),95)])
end
suptitle('# of significant pairs, all neurons all days, experimental') 
figure; 
for mouse = 1:length(ctrl)
    subplot_tight(1,length(ctrl),mouse)
    imagesc(cellPairCount_ctrl{mouse}); axis off;
    caxis([prctile(cellPairCount_ctrl{mouse}(:),5) prctile(cellPairCount_ctrl{mouse}(:),95)])
end
suptitle('# of significant pairs, all neurons all days, ctrl') 

figure; 
for mouse = 1:length(experimental)
    subplot_tight(1,length(experimental),mouse)
    [value, sortIdx] =sort(nanmean(expDataCorrected{mouse}.allPair,2),'descend');
    imagesc(expDataCorrected{mouse}.allPair(sortIdx,:)); axis off;     
    caxis([prctile(expDataCorrected{mouse}.allPair(:),5) prctile(expDataCorrected{mouse}.allPair(:),95)])
end
suptitle('correlation value, experimental') 

figure; 
for mouse = 1:length(experimental)
    subplot_tight(1,length(experimental),mouse)
    [value, sortIdx] =sort(nanmean(expDataCorrected{mouse}.allPairH,2),'descend');
    imagesc(expDataCorrected{mouse}.allPairH(sortIdx,:)); axis off;     
    caxis([0 1])
end
suptitle('correlation significance, experimental') 

figure; 
for mouse = 1:length(ctrl)
    subplot_tight(1,length(ctrl),mouse)
    [value, sortIdx] =sort(nanmean(ctrlData{mouse}.allPair,2),'descend');
    imagesc(ctrlData{mouse}.allPair(sortIdx,:)); axis off;     
    caxis([prctile(ctrlData{mouse}.allPair(:),5) prctile(ctrlData{mouse}.allPair(:),95)])
end
suptitle('correlation value, ctrl') 

figure; 
for mouse = 1:length(ctrl)
    subplot_tight(1,length(ctrl),mouse)
    [value, sortIdx] =sort(nanmean(ctrlData{mouse}.allPairH,2),'descend');
    imagesc(ctrlData{mouse}.allPairH(sortIdx,:)); axis off;     
    caxis([0 1])
end
suptitle('correlation significance, ctrl') 
%% figure 1.5 -- plot the correlation pairs value distribution
figure;  
subplot(2,2,1);hold on;
tempMean = [];  tempMean_ctrl = []; tempMeanSig = [];  tempMean_ctrlSig = [];
for mouse = 1:length(experimental)
    hexp{mouse} = cdfplot(nanmean(expDataCorrected{mouse}.allPair(:,commonDays),2));
    set( hexp{mouse}, 'Color', matlabColors(1,0.6), 'LineWidth', 1.5);
    tempMean(:,mouse) = nanmean(expDataCorrected{mouse}.allPair(:,commonDays),1);
    for j = 1:length(commonDays)
        sigFlag = expDataCorrected{mouse}.allPairH(:,j)==1;
        tempMeanSig(j,mouse) = nanmean(expDataCorrected{mouse}.allPair(sigFlag,commonDays(j)),1);
    end
end
for mouse = 1:length(ctrl)
    hctrl{mouse} = cdfplot(nanmean(ctrlData{mouse}.allPair(:,commonDays),2));
    set( hctrl{mouse}, 'Color', matlabColors(2,0.6), 'LineWidth', 1.5);
    tempMean_ctrl(:,mouse) = nanmean(ctrlData{mouse}.allPair(:,commonDays),1);
    for j = 1:length(commonDays)
        sigFlag = ctrlData{mouse}.allPairH(:,j)==1;
        tempMean_ctrlSig(j,mouse) = nanmean(ctrlData{mouse}.allPair(sigFlag,commonDays(j)),1);
    end
end
legend([hexp{1} hctrl{1}],'experimental','control','Location','Best'); 
xlabel('pairwise correlation'); ylabel('# of neuron'); xlim([-0.03 0.2]);
title('pairwise correlation, all day')
% subplot 2
subplot(4,2,2); hold on;
fn_plotMeanSampleLine(commonDays',tempMean',{'Color',matlabColors(1),'LineWidth',1.5},...
    {'Color',matlabColors(1,0.3),'LineWidth',0.5});
fn_plotMeanSampleLine(commonDays',tempMean_ctrl',{'Color',matlabColors(2),'LineWidth',1.5},...
    {'Color',matlabColors(2,0.3),'LineWidth',0.5});
xlim([commonDays(1) commonDays(end)]); xticklabels([]); ylim([0 0.05]);ylabel('perc pairs')
title('mean perc of significant pair, across days')
subplot(4,2,4); hold on;
[f_line_exp, ~] = fn_plotMeanSampleLine(commonDays',tempMeanSig',{'Color',matlabColors(1),'LineWidth',1.5},...
    {'Color',matlabColors(1,0.3),'LineWidth',0.5});
[f_line_ctrl, ~] = fn_plotMeanSampleLine(commonDays',tempMean_ctrlSig',{'Color',matlabColors(2),'LineWidth',1.5},...
    {'Color',matlabColors(2,0.3),'LineWidth',0.5});
xlim([commonDays(1) commonDays(end)]); ylim([0.05 0.1])
xlabel('days'); ylabel('perc pairs'); 
legend([f_line_exp f_line_ctrl],'experimental','control','Location','Best');
% subplot 3 - day 1
subplot(2,2,3);hold on; selectDay = commonDays(1);
for mouse = 1:length(experimental)
    hexp{mouse} = cdfplot(expDataCorrected{mouse}.allPair(:,selectDay));
    set( hexp{mouse}, 'Color', matlabColors(1,0.6), 'LineWidth', 1.5);
end
for mouse = 1:length(ctrl)
    hctrl{mouse} = cdfplot(ctrlData{mouse}.allPair(:,selectDay));
    set( hctrl{mouse}, 'Color', matlabColors(2,0.6), 'LineWidth', 1.5);
end
legend([hexp{1} hctrl{1}],'experimental','control','Location','Best'); 
xlabel('pairwise correlation'); ylabel('# of neuron'); xlim([-0.03 0.2])
title(['pairwise correlation, day' int2str(selectDay)])
% subplot 4 - day 1
subplot(2,2,4);hold on; selectDay = commonDays(end);
for mouse = 1:length(experimental)
    hexp{mouse} = cdfplot(expDataCorrected{mouse}.allPair(:,selectDay));
    set( hexp{mouse}, 'Color', matlabColors(1,0.6), 'LineWidth', 1.5);
end
for mouse = 1:length(ctrl)
    hctrl{mouse} = cdfplot(ctrlData{mouse}.allPair(:,selectDay));
    set( hctrl{mouse}, 'Color', matlabColors(2,0.6), 'LineWidth', 1.5);
end
legend([hexp{1} hctrl{1}],'experimental','control','Location','Best'); 
xlabel('pairwise correlation'); ylabel('# of neuron'); xlim([-0.03 0.2])
title(['pairwise correlation, day' int2str(selectDay)])
%% figure 1.6 -- Plot the disctribution of stable (many day significant) vs. transient (few day significant) pairs
expMultidayPair = {}; ctrlMultidayPair = {}; expN = []; ctrlN = [];
expMultidayPairValue = {}; ctrlMultidayPairValue = {}; expHML = {}; ctrlHML = {};
expHMLflag = {}; ctrlHMLflag = {};
figure; 
for mouse = 1:length(experimental)
    expMultidayPair{mouse} = nansum(expDataCorrected{mouse}.allPairH(:,commonDays),2);
    [N, edges, bins] = histcounts(expMultidayPair{mouse},[0.5 4.5 9.5 15.5]); 
    expHMLflag{mouse} = bins;
    tempPair = expDataCorrected{mouse}.allPair(:,commonDays); 
    tempPair(expDataCorrected{mouse}.allPairH(:,commonDays)==0) = nan;
    for i = 1:3; expHML{mouse,i} = tempPair(bins==i,:); end
    [N, edges, bins] = histcounts(expMultidayPair{mouse},commonDays(1)-1.5:commonDays(end)+0.5); 
    expN(:,mouse) = cumsum(N)./sum(N);
end
for mouse = 1:length(ctrl)
    ctrlMultidayPair{mouse} = nansum(ctrlData{mouse}.allPairH(:,commonDays),2);
    [N, edges, bins] = histcounts(ctrlMultidayPair{mouse},[0.5 4.5 9.5 15.5]); 
    ctrlHMLflag{mouse} = bins;
    tempPair = ctrlData{mouse}.allPair(:,commonDays); 
    tempPair(ctrlData{mouse}.allPairH(:,commonDays)==0) = nan;
    for i = 1:3; ctrlHML{mouse,i} = tempPair(bins==i,:); end
    [N, edges, bins] = histcounts(ctrlMultidayPair{mouse},commonDays(1)-1.5:commonDays(end)+0.5);
    ctrlN(:,mouse) = cumsum(N)/sum(N);
end
subplot(1,2,1); hold on; 
[f_line_exp, ~] = fn_plotMeanSampleLine([0;commonDays]',expN',{'Color',matlabColors(1),'LineWidth',1.5},...
    {'Color',matlabColors(1,0.3),'LineWidth',0.5});
[f_line_ctrl, ~] = fn_plotMeanSampleLine([0;commonDays]',ctrlN',{'Color',matlabColors(2),'LineWidth',1.5},...
    {'Color',matlabColors(2,0.3),'LineWidth',0.5});
xlabel('# of days that a pair is significant'); ylabel('# pairs'); legend([f_line_exp,f_line_ctrl],'experimental','ctrl','Location','Best')
title('Distribution of stability of pairwise correlation')
subplot(1,2,2); hold on;
[f_line_exp, ~] =fn_plotMeanErrorbar([0;commonDays]',expN',matlabColors(1),...
    {'Color',matlabColors(1),'LineWidth',1.5},{'faceAlpha',0.3,'LineStyle','none'});
[f_line_ctrl, ~] = fn_plotMeanErrorbar([0;commonDays]',ctrlN',matlabColors(2),...
    {'Color',matlabColors(2),'LineWidth',1.5},{'faceAlpha',0.3,'LineStyle','none'});
xlabel('# of days that a pair is significant'); ylabel('# pairs'); legend([f_line_exp,f_line_ctrl],'experimental','ctrl','Location','Best')
title('Distribution of stability of pairwise correlation')
%% figure 1.7 -- plot the 
figure;
expHMLAvg_all = []; expHML_dayStart = []; expHML_dayEnd = [];
subplot(2,3,1)
selectEdge = (-0.1-0.00025):0.0005:(1+0.00025);
selectBin = -0.1:0.0005:1;
for mouse = 1:length(experimental)
    temp = [];
    for i = 1:3
        [N, edges, bins] = histcounts(nanmean(expHML{mouse,i},2),selectEdge);
        expHMLAvg_all(:,mouse,i) = cumsum(N)./ sum(N);
        [N, edges, bins] = histcounts(expHML{mouse,i}(:,commonDays(1)),selectEdge);
        expHML_dayStart (:,mouse,i) = cumsum(N)./ sum(N);
        [N, edges, bins] = histcounts(expHML{mouse,i}(:,commonDays(end)),selectEdge);
        expHML_dayEnd (:,mouse,i) = cumsum(N)./ sum(N);
    end
end
[f_line_exp_low, ~] =fn_plotMeanSampleLine(selectBin,expHMLAvg_all(:,:,1)',...
    {'Color',matlabColors(1),'LineWidth',1.5},{'Color',matlabColors(1,0.3),'LineWidth',0.5});
[f_line_exp_mid, ~] =fn_plotMeanSampleLine(selectBin,expHMLAvg_all(:,:,2)',...
    {'Color',matlabColors(2),'LineWidth',1.5},{'Color',matlabColors(2,0.3),'LineWidth',0.5});
[f_line_exp_high, ~] =fn_plotMeanSampleLine(selectBin,expHMLAvg_all(:,:,3)',...
    {'Color',matlabColors(3),'LineWidth',1.5},{'Color',matlabColors(3,0.3),'LineWidth',0.5});
xlim([0.02 0.16]); xlabel('mean correlation value'); ylabel('# of neurons')
legend([f_line_exp_low,f_line_exp_mid,f_line_exp_high],'transient','intermediate','stable','Location','Best')
title('Experimental, mean correlation of all day')

subplot(2,3,2)
[f_line_exp, ~] =fn_plotMeanSampleLine(selectBin,expHML_dayStart(:,:,1)',...
    {'Color',matlabColors(1),'LineWidth',1.5},{'Color',matlabColors(1,0.3),'LineWidth',0.5});
[f_line_exp, ~] =fn_plotMeanSampleLine(selectBin,expHML_dayStart(:,:,2)',...
    {'Color',matlabColors(2),'LineWidth',1.5},{'Color',matlabColors(2,0.3),'LineWidth',0.5});
[f_line_exp, ~] =fn_plotMeanSampleLine(selectBin,expHML_dayStart(:,:,3)',...
    {'Color',matlabColors(3),'LineWidth',1.5},{'Color',matlabColors(3,0.3),'LineWidth',0.5});
xlim([0.02 0.16]); xlabel('correlation value'); ylabel('# of neurons')
title('Experimental, correlation value of day 1')

subplot(2,3,3)
[f_line_exp, ~] =fn_plotMeanSampleLine(selectBin,expHML_dayEnd(:,:,1)',...
    {'Color',matlabColors(1),'LineWidth',1.5},{'Color',matlabColors(1,0.3),'LineWidth',0.5});
[f_line_exp, ~] =fn_plotMeanSampleLine(selectBin,expHML_dayEnd(:,:,2)',...
    {'Color',matlabColors(2),'LineWidth',1.5},{'Color',matlabColors(2,0.3),'LineWidth',0.5});
[f_line_exp, ~] =fn_plotMeanSampleLine(selectBin,expHML_dayEnd(:,:,3)',...
    {'Color',matlabColors(3),'LineWidth',1.5},{'Color',matlabColors(3,0.3),'LineWidth',0.5});
xlim([0.02 0.16]); xlabel('correlation value'); ylabel('# of neurons')
title('Experimental, correlation value of day 15')

ctrlHMLAvg_all = []; ctrlHML_dayStart = []; ctrlHML_dayEnd = [];
subplot(2,3,4)
selectEdge = (-0.1-0.00025):0.0005:(1+0.00025);
selectBin = -0.1:0.0005:1;
for mouse = 1:length(ctrl)
    temp = [];
    for i = 1:3
        [N, edges, bins] = histcounts(nanmean(ctrlHML{mouse,i},2),selectEdge);
        ctrlHMLAvg_all(:,mouse,i) = cumsum(N)./ sum(N);
        [N, edges, bins] = histcounts(ctrlHML{mouse,i}(:,commonDays(1)),selectEdge);
        ctrlHML_dayStart (:,mouse,i) = cumsum(N)./ sum(N);
        [N, edges, bins] = histcounts(ctrlHML{mouse,i}(:,commonDays(end)),selectEdge);
        ctrlHML_dayEnd (:,mouse,i) = cumsum(N)./ sum(N);
    end
end
[f_line_exp, ~] =fn_plotMeanSampleLine(selectBin,ctrlHMLAvg_all(:,:,1)',...
    {'Color',matlabColors(1),'LineWidth',1.5},{'Color',matlabColors(1,0.3),'LineWidth',0.5});
[f_line_exp, ~] =fn_plotMeanSampleLine(selectBin,ctrlHMLAvg_all(:,:,2)',...
    {'Color',matlabColors(2),'LineWidth',1.5},{'Color',matlabColors(2,0.3),'LineWidth',0.5});
[f_line_exp, ~] =fn_plotMeanSampleLine(selectBin,ctrlHMLAvg_all(:,:,3)',...
    {'Color',matlabColors(3),'LineWidth',1.5},{'Color',matlabColors(3,0.3),'LineWidth',0.5});
xlim([0.02 0.16]); xlabel('mean correlation value'); ylabel('# of neurons')
title('Control, mean correlation of all day')
subplot(2,3,5)
[f_line_exp, ~] =fn_plotMeanSampleLine(selectBin,ctrlHML_dayStart(:,:,1)',...
    {'Color',matlabColors(1),'LineWidth',1.5},{'Color',matlabColors(1,0.3),'LineWidth',0.5});
[f_line_exp, ~] =fn_plotMeanSampleLine(selectBin,ctrlHML_dayStart(:,:,2)',...
    {'Color',matlabColors(2),'LineWidth',1.5},{'Color',matlabColors(2,0.3),'LineWidth',0.5});
[f_line_exp, ~] =fn_plotMeanSampleLine(selectBin,ctrlHML_dayStart(:,:,3)',...
    {'Color',matlabColors(3),'LineWidth',1.5},{'Color',matlabColors(3,0.3),'LineWidth',0.5});
xlim([0.02 0.16]); xlabel('correlation value'); ylabel('# of neurons')
title('Control, correlation value of day 1')
subplot(2,3,6)
[f_line_exp, ~] =fn_plotMeanSampleLine(selectBin,ctrlHML_dayEnd(:,:,1)',...
    {'Color',matlabColors(1),'LineWidth',1.5},{'Color',matlabColors(1,0.3),'LineWidth',0.5});
[f_line_exp, ~] =fn_plotMeanSampleLine(selectBin,ctrlHML_dayEnd(:,:,2)',...
    {'Color',matlabColors(2),'LineWidth',1.5},{'Color',matlabColors(2,0.3),'LineWidth',0.5});
[f_line_exp, ~] =fn_plotMeanSampleLine(selectBin,ctrlHML_dayEnd(:,:,3)',...
    {'Color',matlabColors(3),'LineWidth',1.5},{'Color',matlabColors(3,0.3),'LineWidth',0.5});
xlim([0.02 0.16]); xlabel('correlation value'); ylabel('# of neurons')
title('Control, correlation value of day 15')
%% figure 1.7 -- Compute the rate of change between consecutive days in experimental vs. ctrl
figure; expChangeRate = struct();ctrlChangeRate = struct();
for mouse = 1:length(experimental)
    for j = 1:length(commonDays)-1
        nNeuron = size(expDataCorrected{mouse}.allPairH,1);
        expChangeRate.forwardChange{mouse}(:,j) = expDataCorrected{mouse}.allPairH(:,commonDays(j)+1)...
            -expDataCorrected{mouse}.allPairH(:,commonDays(j));
        if ~all(isnan(expChangeRate.forwardChange{mouse}(:,j)))
            expChangeRate.forwardChangeSum(mouse,j) = nanmean(abs(expChangeRate.forwardChange{mouse}(:,j)));
            expChangeRate.forwardChangeGainSum(mouse,j) = nansum(expChangeRate.forwardChange{mouse}(:,j)==1)/nNeuron;
            expChangeRate.forwardChangeLossSum(mouse,j) = nansum(expChangeRate.forwardChange{mouse}(:,j)==-1)/nNeuron;
            expChangeRate.forwardChangeStaySum(mouse,j)  = nansum(expChangeRate.forwardChange{mouse}(:,j)==0)/nNeuron;
        else
            expChangeRate.forwardChangeSum(mouse,j) = nan; 
            expChangeRate.forwardChangeGainSum(mouse,j)  = nan;
            expChangeRate.forwardChangeLossSum(mouse,j)  = nan;
            expChangeRate.forwardChangeStaySum(mouse,j)  = nan;
        end
    end
end
for mouse = 1:length(ctrl)
    for j = 1:length(commonDays)-1
        nNeuron = size(ctrlData{mouse}.allPairH,1);
        ctrlChangeRate.forwardChange{mouse}(:,j) = ctrlData{mouse}.allPairH(:,commonDays(j)+1)...
            -ctrlData{mouse}.allPairH(:,commonDays(j));
        if ~all(isnan(ctrlChangeRate.forwardChange{mouse}(:,j)))
            ctrlChangeRate.forwardChangeSum(mouse,j) = nanmean(abs(ctrlChangeRate.forwardChange{mouse}(:,j)));
            ctrlChangeRate.forwardChangeGainSum(mouse,j) = nansum(ctrlChangeRate.forwardChange{mouse}(:,j)==1)/nNeuron;
            ctrlChangeRate.forwardChangeLossSum(mouse,j) = nansum(ctrlChangeRate.forwardChange{mouse}(:,j)==-1)/nNeuron;
            ctrlChangeRate.forwardChangeStaySum(mouse,j)  = nansum(ctrlChangeRate.forwardChange{mouse}(:,j)==0)/nNeuron;
        else
            ctrlChangeRate.forwardChangeSum(mouse,j) = nan; 
            ctrlChangeRate.forwardChangeGainSum(mouse,j)  = nan;
            ctrlChangeRate.forwardChangeLossSum(mouse,j)  = nan;
            ctrlChangeRate.forwardChangeStaySum(mouse,j)  = nan;
        end
    end
end
subplot(2,2,1)
fn_plotMeanSampleLine(commonDays(1:end-1)',expChangeRate.forwardChangeSum,{'Color',matlabColors(1),'LineWidth',1.5},...
    {'Color',matlabColors(1,0.3),'LineWidth',0.5});
fn_plotMeanSampleLine(commonDays(1:end-1)',ctrlChangeRate.forwardChangeSum,{'Color',matlabColors(2),'LineWidth',1.5},...
    {'Color',matlabColors(2,0.3),'LineWidth',0.5});
xlim([commonDays(1) commonDays(end-1)]); title('perc of pair change');xlabel('days'); ylabel('perc of pair')
subplot(2,2,2)
fn_plotMeanSampleLine(commonDays(1:end-1)',expChangeRate.forwardChangeGainSum,{'Color',matlabColors(1),'LineWidth',1.5},...
    {'Color',matlabColors(1,0.3),'LineWidth',0.5});
fn_plotMeanSampleLine(commonDays(1:end-1)',ctrlChangeRate.forwardChangeGainSum,{'Color',matlabColors(2),'LineWidth',1.5},...
    {'Color',matlabColors(2,0.3),'LineWidth',0.5});
xlim([commonDays(1) commonDays(end-1)]); title('perc of pair gain signif');xlabel('days'); ylabel('perc of pair')
subplot(2,2,3)
fn_plotMeanSampleLine(commonDays(1:end-1)',expChangeRate.forwardChangeLossSum,{'Color',matlabColors(1),'LineWidth',1.5},...
    {'Color',matlabColors(1,0.3),'LineWidth',0.5});
fn_plotMeanSampleLine(commonDays(1:end-1)',ctrlChangeRate.forwardChangeLossSum,{'Color',matlabColors(2),'LineWidth',1.5},...
    {'Color',matlabColors(2,0.3),'LineWidth',0.5});
xlim([commonDays(1) commonDays(end-1)]); title('perc of pair lose signif');xlabel('days'); ylabel('perc of pair')
subplot(2,2,4)
fn_plotMeanSampleLine(commonDays(1:end-1)',expChangeRate.forwardChangeStaySum,{'Color',matlabColors(1),'LineWidth',1.5},...
    {'Color',matlabColors(1,0.3),'LineWidth',0.5});
fn_plotMeanSampleLine(commonDays(1:end-1)',ctrlChangeRate.forwardChangeStaySum,{'Color',matlabColors(2),'LineWidth',1.5},...
    {'Color',matlabColors(2,0.3),'LineWidth',0.5});
xlim([commonDays(1) commonDays(end-1)]); title('perc of pair stay'); xlabel('days'); ylabel('perc of pair')

%% figure 1.8 -- Correlation of significant correlation pairs / correlation value across days
figure;
for mouse = 1:length(experimental)
    subplot(2,3,mouse);
    imagesc(corr(expDataCorrected{mouse}.allPairH)); caxis([0 0.5]); title(experimental{mouse})
end
figure;
for mouse = 1:length(ctrl)
    subplot(1,2,mouse);
    imagesc(corr(ctrlData{mouse}.allPairH));caxis([0 0.5]);  title(ctrl{mouse})
end

figure;
for mouse = 1:length(experimental)
    subplot(2,3,mouse);
    imagesc(corr(expDataCorrected{mouse}.allPair)); caxis([0 0.8]); title(experimental{mouse})
end
figure;
for mouse = 1:length(ctrl)
    subplot(1,2,mouse);
    imagesc(corr(ctrlData{mouse}.allPair));caxis([0 0.8]);  title(ctrl{mouse})
end
%% figure 1.9 -- Correlation of significant correlation pairs, aligned to each day
selectDay = -7:7;
[shiftMatMask] = getAlignToDayMatrix(commonDays,selectDay);
expCorr = zeros(length(commonDays),length(commonDays),length(experimental));
expCorrAligned = zeros(size(expCorr));
for mouse = 1:length(experimental)
    expCorr(:,:,mouse)= corr(expDataCorrected{mouse}.allPairH(:,commonDays));
    meanCommonDay = round((commonDays(1)+commonDays(end))/2);
    for j = 1:length(commonDays)
        expCorrAligned(j,:,mouse) = circshift(expCorr(j,:,mouse),meanCommonDay-j);
    end
    temp = expCorrAligned(:,:,mouse); temp(shiftMatMask==0) = nan;
    expCorrAligned(:,:,mouse) = temp;
end
expCorrAlignedMean = squeeze(nanmean(expCorrAligned,1));
ctrlCorr = zeros(length(commonDays),length(commonDays),length(ctrl));
ctrlCorrAligned = zeros(size(ctrlCorr));
for mouse = 1:length(ctrl)
    ctrlCorr(:,:,mouse)= corr(ctrlData{mouse}.allPairH(:,commonDays));
    meanCommonDay = round((commonDays(1)+commonDays(end))/2);
    for j = 1:length(commonDays)
        ctrlCorrAligned(j,:,mouse) = circshift(ctrlCorr(j,:,mouse),meanCommonDay-j);
    end
    temp = ctrlCorrAligned(:,:,mouse); temp(shiftMatMask==0) = nan;
    ctrlCorrAligned(:,:,mouse) = temp;
end
ctrlCorrAlignedMean = squeeze(nanmean(ctrlCorrAligned,1));
figure; subplot(1,3,1)
[f_line_exp, ~] =fn_plotMeanSampleLine(selectDay,expCorrAlignedMean',...
    {'Color',matlabColors(1),'LineWidth',1.5},{'Color',matlabColors(1,0.3),'LineWidth',1.5});
[f_line_ctrl, ~] = fn_plotMeanSampleLine(selectDay,ctrlCorrAlignedMean',...
    {'Color',matlabColors(2),'LineWidth',1.5},{'Color',matlabColors(2,0.3),'LineWidth',1.5});
xlim([selectDay(1) selectDay(end)]); ylim([0 1])
xlabel('days'); ylabel('corr between significant pairs'); title('all days')
legend([f_line_exp f_line_ctrl],'experimental','control','Location','Best');
subplot(1,3,2);
expMean = squeeze(nanmean(expCorrAligned(1:7,5:11,:),1));
ctrlMean = squeeze(nanmean(ctrlCorrAligned(1:7,5:11,:),1));
[f_line_exp, ~] =fn_plotMeanSampleLine(-3:3,expMean',...
    {'Color',matlabColors(1),'LineWidth',1.5},{'Color',matlabColors(1,0.3),'LineWidth',1.5});
[f_line_ctrl, ~] = fn_plotMeanSampleLine(-3:3,ctrlMean',...
    {'Color',matlabColors(2),'LineWidth',1.5},{'Color',matlabColors(2,0.3),'LineWidth',1.5});
xlim([-3 3]); ylim([0 1])
xlabel('days'); ylabel('corr between significant pairs'); title('day1-7')
legend([f_line_exp f_line_ctrl],'experimental','control','Location','Best');
subplot(1,3,3);
expMean = squeeze(nanmean(expCorrAligned(8:14,5:11,:),1));
ctrlMean = squeeze(nanmean(ctrlCorrAligned(8:14,5:11,:),1));
[f_line_exp, ~] =fn_plotMeanSampleLine(-3:3,expMean',...
    {'Color',matlabColors(1),'LineWidth',1.5},{'Color',matlabColors(1,0.3),'LineWidth',1.5});
[f_line_ctrl, ~] = fn_plotMeanSampleLine(-3:3,ctrlMean',...
    {'Color',matlabColors(2),'LineWidth',1.5},{'Color',matlabColors(2,0.3),'LineWidth',1.5});
xlim([-3 3]); ylim([0 1])
xlabel('days'); ylabel('corr between significant pairs'); title('day8-14')
legend([f_line_exp f_line_ctrl],'experimental','control','Location','Best');

%% figure 1.10 -- Correlation of significant correlation values, aligned to each day


%% figure 2.1 -- 
expHML_SI = {}; expHML_SImean = {};
for mouse = 1:length(experimental)
    
    selectSI = expDataSI{mouse}.SSI(expData{mouse}.selectNeuron,commonDays);
    nNeuron = size(selectSI,1);
    % Use the line below to test if the code works fine. The output pair1 and pair 2 should correspond to cell 1 and 2 for that pair
    % selectSI = repmat((1:nNeuron)',[1 length(commonDays)]); 
    selectSIPair1_temp = repmat(reshape(selectSI, [1 nNeuron length(commonDays)]),[nNeuron 1 1]);
    selectSIPair2_temp = repmat(reshape(selectSI, [nNeuron 1 length(commonDays)]),[1 nNeuron 1]);
    selectSIPair1 = []; selectSIPair2 = [];
    for day = 1:length(commonDays)
        temp = selectSIPair1_temp(:,:,day); temp(expData{mouse}.tridIdx==0) = [];
        selectSIPair1(:,day) = temp;
        temp = selectSIPair2_temp(:,:,day); temp(expData{mouse}.tridIdx==0) = [];
        selectSIPair2(:,day) = temp;
    end
    % now selectSIPair1 and selectSIPair2 should be the SI of cell 1 and
    % cell 2 for each neural pair
    
    for i = 1:3
        HMLflag = expHMLflag{mouse}==i;
        signifMask = isnan(expHML{mouse,i});
        binSI1 = selectSIPair1(HMLflag,:); binSI1(signifMask) = nan;
        binSI2 = selectSIPair2(HMLflag,:); binSI2(signifMask) = nan;
        expHML_SI{mouse,i,1} = binSI1; expHML_SImean{mouse,i,1} = nanmean(expHML_SI{mouse,i,1},2);
        expHML_SI{mouse,i,2} = binSI2; expHML_SImean{mouse,i,2} = nanmean(expHML_SI{mouse,i,2},2);
    end
end

figure; 
for mouse = 1:length(experimental)
    subplot(2,3,mouse); hold on; 
    for i = 1:3
        scatter(expHML_SImean{mouse,i,1},expHML_SImean{mouse,i,2},5,matlabColors(i,0.3),'filled');
    end
end
%%
figure; day = 1;
for mouse = 1:length(experimental)
    subplot(2,3,mouse); hold on; 
    h = {};
    for i = 1:3
        diffSSI = expHML_SI{mouse,i,1} - expHML_SI{mouse,i,2};
        h{i} = cdfplot(abs(diffSSI(:,day)));
        set( h{i}, 'Color', matlabColors(i,0.6), 'LineWidth', 2);
    end
    title(experimental{mouse}); xlim([0 0.8]); xlabel('difference in stimulus selectivity');ylabel('# pairs')
    legend([h{1} h{2} h{3}], 'transient','intermediate','stable')
end

%%

for mouse = 1:length(experimental)
    figure;
    selectSSI = expDataSI{mouse}.SSI(expData{mouse}.selectNeuron,commonDays);
    selectT = expDataSI{mouse}.ttestT(expData{mouse}.selectNeuron,commonDays);
    selectF = expDataSI{mouse}.ttestF(expData{mouse}.selectNeuron,commonDays);
    

    nNeuron = size(selectSSI,1);
    %selectSSI = repmat((1:nNeuron)',[1 length(commonDays)]); 
    selectSIPair1_temp = repmat(reshape(selectSSI, [1 nNeuron length(commonDays)]),[nNeuron 1 1]);
    selectSIPair2_temp = repmat(reshape(selectSSI, [nNeuron 1 length(commonDays)]),[1 nNeuron 1]);
    pairDiff = selectSIPair1_temp - selectSIPair2_temp; % here is the difference of SSI of cell1 - cell2
    selectPairDiff = [];
    for i = 1:length(commonDays)
        temp = pairDiff(:,:,i);
        temp(expData{mouse}.tridIdx==0) = [];
        selectPairDiff(:,day) = temp;
        subplot_tight(3,5,i);
        scatter(abs(temp), expData{mouse}.allPair(:,i),1,'filled');
        xlim([-0.2 2]); ylim([-0.05 0.2])
        axis off
    end
    
    
end

%% old figures
for i = 1:length(expData)
    tic;
    nDays = min(size(expDataSI{i}.SI,2), size(expData{i}.allPair,2));
    n = [];nSSI = [];
    for j = 1:nDays
        signifIdx = find(expData{i}.allPairH(:,j));
        nTemp = zeros(length(signifIdx),2);
        nSSITemp = zeros(length(signifIdx),2);
        for k = 1:length(signifIdx)
            [temp1,temp2] = find(expData{i}.tridIdx==signifIdx(k));
            nTemp(k,:) = [temp1 temp2];
            tempSSI = expDataSI{i}.SSI(expData{i}.selectNeuron,j);
            nSSITemp(k,:) = [tempSSI(temp1) tempSSI(temp2)];
        end
        n{j} = nTemp; nSSI{j} = nSSITemp; 
    end
    % get SSI of neuron pairs that are active all days 
    nAllDay = [];nSSIAllDay = [];
    signifIdxAllday = find(sum(expData{i}.allPairH(:,1:nDays),2) == nDays);
    nTemp = zeros(length(signifIdxAllday),2);
    nSSITemp = zeros(length(signifIdxAllday),2);
    for k = 1:length(signifIdxAllday)
        [temp1,temp2] = find(expData{i}.tridIdx==signifIdxAllday(k));
        nAllDay(k,:) = [temp1 temp2];
        tempSSI = expDataSI{i}.SSI(expData{i}.selectNeuron,:);
        nSSIAllDay(k,:,:) = tempSSI([temp1 temp2],:);
    end
    save([expDatapath{i} '\' experimental{i} '_baseline_corr_SI.mat'],...
        'n','nSSI','nDays','nAllDay','nSSIAllDay');
    toc;
end

for i = 1:length(ctrlData)
    tic;
    nDays = min(size(ctrlDataSI{i}.SI,2), size(ctrlData{i}.allPair,2));
    n = [];nSSI = [];
    for j = 1:nDays
        signifIdx = find(ctrlData{i}.allPairH(:,j));
        nTemp = zeros(length(signifIdx),2);
        nSSITemp = zeros(length(signifIdx),2);
        for k = 1:length(signifIdx)
            [temp1,temp2] = find(ctrlData{i}.tridIdx==signifIdx(k));
            nTemp(k,:) = [temp1 temp2];
            tempSSI = ctrlDataSI{i}.SSI(ctrlData{i}.selectNeuron,j);
            nSSITemp(k,:) = [tempSSI(temp1) tempSSI(temp2)];
        end
        n{j} = nTemp; nSSI{j} = nSSITemp; 
    end
    % get SSI of neuron pairs that are active all days 
    nAllDay = [];nSSIAllDay = [];
    signifIdxAllday = find(sum(ctrlData{i}.allPairH(:,1:nDays),2) == nDays);
    nTemp = zeros(length(signifIdxAllday),2);
    nSSITemp = zeros(length(signifIdxAllday),2);
    for k = 1:length(signifIdxAllday)
        [temp1,temp2] = find(ctrlData{i}.tridIdx==signifIdxAllday(k));
        nAllDay(k,:) = [temp1 temp2];
        tempSSI = ctrlDataSI{i}.SSI(ctrlData{i}.selectNeuron,:);
        nSSIAllDay(k,:,:) = tempSSI([temp1 temp2],:);
    end
        
    save([ctrlDatapath{i} '\' ctrl{i} '_baseline_corr_SI.mat'],...
        'n','nSSI','nDays','nAllDay','nSSIAllDay');
    toc;
end
%% load new data
clear;
experimental = {'cd017','cd036','cd037','cd042','cd044'};
for i = 1:length(experimental)
    expDatapath{i} = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' experimental{i}];
    tempData = load([expDatapath{i} '\' experimental{i} '_baseline_data_pairwiseCorr.mat']);
    expData{i} = tempData;
    tempData = load([expDatapath{i} '\' experimental{i} '_baseline_corr_SI.mat']);
    expData{i} = catstruct(expData{i},tempData);
end

ctrl = {'cd019','cd041'};
for i = 1:length(ctrl)
    ctrlDatapath{i} = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' ctrl{i}];
    tempData = load([ctrlDatapath{i} '\' ctrl{i} '_baseline_data_pairwiseCorr.mat']);
    ctrlData{i} = tempData;
    tempData = load([ctrlDatapath{i} '\' ctrl{i} '_baseline_corr_SI.mat']);
    ctrlData{i} = catstruct(ctrlData{i},tempData);
end
%%
for i = 1:length(experimental)
    figure; hold on;
    tempColors = [linspace(0,1,expData{i}.nDays);zeros(2,expData{i}.nDays)]; 
    for j = 1:expData{i}.nDays
        h = cdfplot(sum(abs(expData{i}.nSSI{j}),2)); 
        set( h, 'Color', tempColors(:,j), 'LineWidth', 1.5);
    end
    title(experimental{i})
end

for i = 1:length(ctrl)
    figure; hold on;
    tempColors = [linspace(0,1,ctrlData{i}.nDays);zeros(2,ctrlData{i}.nDays)]; 
    for j = 1:ctrlData{i}.nDays
        h = cdfplot(abs(ctrlData{i}.nSSI{j}(:,1)-ctrlData{i}.nSSI{j}(:,2))); 
        set( h, 'Color', tempColors(:,j), 'LineWidth', 1.5);
    end
    title(ctrl{i})
end
%% Pairs with significant correlation on all days
expPair = []; expSSIdiff = []; expSSIsum = []; maxDay = 15;
for i = 1:length(experimental)
    signifIdxAllday = find(sum(expData{i}.allPairH(:,1:expData{i}.nDays),2) == expData{i}.nDays);
    temp = expData{i}.allPair(signifIdxAllday,1:maxDay);
    corrMatExp(:,:,i) = corr(temp);
    expPair = cat(1,expPair,temp);
    temp = squeeze(abs(expData{i}.nSSIAllDay(:,1,1:maxDay)-expData{i}.nSSIAllDay(:,2,1:maxDay)));
    expSSIdiff = cat(1,expSSIdiff,temp);
    temp = squeeze(abs(expData{i}.nSSIAllDay(:,1,1:maxDay))+abs(expData{i}.nSSIAllDay(:,2,1:maxDay)));
    expSSIsum = cat(1,expSSIsum,temp);
end
figure; subplot(1,3,1); imagesc(expPair); caxis([0 0.3]);yticklabels([])
subplot(1,3,2); imagesc(expSSIdiff); caxis([0 0.8]);yticklabels([])
subplot(1,3,3); imagesc(expSSIsum); caxis([0 1.4]);yticklabels([])
[~,peakIdx] = max(expSSIdiff);
highFlag = sum(expSSIsum,2) > 15; lowFlag = sum(expSSIsum,2) <= 15;
figure; plot(mean(expPair(highFlag,:),1)); hold on; plot(mean(expPair(lowFlag,:),1));
%[basis, varExp, proj, covMat] = fn_pca(expPairAll);

ctrlPair = []; ctrlSSIdiff = [];ctrlSSIsum = [];maxDay = 14;
for i = 1:length(ctrl)
    signifIdxAllday = find(sum(ctrlData{i}.allPairH(:,1:ctrlData{i}.nDays),2) == ctrlData{i}.nDays);
    temp = ctrlData{i}.allPair(signifIdxAllday,1:maxDay);
    corrMatCtrl(:,:,i) = corr(temp);
    ctrlPair = cat(1,ctrlPair,temp);
    temp = squeeze(abs(ctrlData{i}.nSSIAllDay(:,1,1:maxDay)-ctrlData{i}.nSSIAllDay(:,2,1:maxDay)));
    ctrlSSIdiff = cat(1,ctrlSSIdiff,temp);
    temp = squeeze(abs(ctrlData{i}.nSSIAllDay(:,1,1:maxDay))+abs(ctrlData{i}.nSSIAllDay(:,2,1:maxDay)));
    ctrlSSIsum = cat(1,ctrlSSIsum,temp);
end
figure; subplot(1,3,1); imagesc(ctrlPair); caxis([0 0.3]);yticklabels([])
subplot(1,3,2); imagesc(ctrlSSIdiff); caxis([0 0.8]);yticklabels([])
subplot(1,3,3); imagesc(ctrlSSIsum); caxis([0 1.4]);yticklabels([])


highFlag = sum(ctrlSSIsum,2) > 15; lowFlag = sum(ctrlSSIsum,2) <= 15;
figure; plot(mean(ctrlPair(highFlag,:),1)); hold on; plot(mean(ctrlPair(lowFlag,:),1));
%[basis, varExp, proj, covMat] = fn_pca(ctrlPairAll);


figure; for i = 1:5; subplot(2,3,i); imagesc(corrMatExp(:,:,i));caxis([0.3 0.8]); end
subplot(2,3,6); imagesc(mean(corrMatExp,3))

figure; for i = 1:2; subplot(1,3,i); imagesc(corrMatCtrl(:,:,i));caxis([0.3 0.8]); end
subplot(1,3,3); imagesc(mean(corrMatCtrl,3))

%% Figures
maxDay = 15; nExp = [];
for i = 1:length(experimental)
    for j = 1:maxDay
        nExp(i,j) = size(expData{i}.n{j},1);
    end
    nExp = nExp ./ repmat(nExp(:,1),1,maxDay);
end

maxDay = 14; nCtrl =[];
for i = 1:length(ctrl)
    for j = 1:maxDay
        nCtrl(i,j) = size(ctrlData{i}.n{j},1);
    end
    nCtrl = nCtrl ./ repmat(nCtrl(:,1),1,maxDay);
end
figure; hold on;

plot(mean(nExp),'Color',matlabColors(1),'LineWidth',2);
plot(mean(nCtrl),'Color',matlabColors(2),'LineWidth',2);
plot(nExp','Color',(matlabColors(1)*0.2 + [1 1 1] * 0.8));
plot(nCtrl','Color',(matlabColors(2)*0.2 + [1 1 1] * 0.8));

xlabel('day'); ylabel('# of pairs');title('# of significant pairs, normalized')
legend('behav','ctrl')

figure; hold on;plot(mean(expPair,1))
plot(mean(ctrlPair,1))

maxDay = 15;
colorScale = repmat(linspace(1,0,maxDay),3,1) .* repmat([0;0;0],1,maxDay) ...
    + repmat(linspace(0,1,maxDay),3,1) .* repmat([1;0.5;0.5],1,maxDay) ;
figure; hold on;
for i = 1:maxDay
h = cdfplot(expPair(:,i)); set( h, 'Color', colorScale(:,i), 'LineWidth', 1.5);
end
xlim([0 0.5]);xlabel('correlation'); title('behavior')
maxDay = 14;
colorScale = repmat(linspace(1,0,maxDay),3,1) .* repmat([0;0;0],1,maxDay) ...
    + repmat(linspace(0,1,maxDay),3,1) .* repmat([1;0.5;0.5],1,maxDay) ;
figure; hold on;
for i = 1:maxDay
h = cdfplot(ctrlPair(:,i)); set( h, 'Color', colorScale(:,i), 'LineWidth', 1.5);
end
xlim([0 0.5]); xlabel('correlation'); title('ctrl')

figure; hold on;
maxDay = 15;plot(1:maxDay,mean(expPair,1),'Color',matlabColors(1));
maxDay = 14;plot(1:maxDay,mean(ctrlPair,1),'Color',matlabColors(2));
maxDay = 15;
f = fn_plotFillErrorbar(1:maxDay,mean(expPair,1), ...
    std(expPair,1),matlabColors(1),'LineStyle','none');
f.FaceAlpha = 0.3;
maxDay = 14;
f = fn_plotFillErrorbar(1:maxDay,mean(ctrlPair,1), ...
    std(ctrlPair,1),matlabColors(2),'LineStyle','none');
f.FaceAlpha = 0.3;
xlim([1 maxDay]); legend('behavior','ctrl'); xlabel('day'); ylabel('mean corr')


function [covMatSignif, cellPairIdx, cellPairCount] = getSignificantMatrix (animalList, animalData, zscoreThreshold)
    covMatSignif = {}; cellPairIdx = {}; cellPairCount = {};
    animalData = removeBadSession(animalList, animalData);

    for mouse = 1:length(animalList)
        tic; 
        nNeuron = size(animalData{mouse}.tridIdx,1);
        nDays = length(animalData{mouse}.baselineDay);
        cellPairCount{mouse} = zeros(nNeuron,nDays);
        covMatSignif{mouse} = zeros(nNeuron,nNeuron,nDays);
        for day = 1:nDays
            signifIdx = find(animalData{mouse}.allPairH(:,day) ==1);
            if exist('zscoreThreshold')
                signifIdx = find(animalData{mouse}.allPairZscore(:,day) > zscoreThreshold);
            end

            cellPairIdx{mouse}{day} = zeros(2,length(signifIdx));
            for i = 1:length(signifIdx)
                [idx1, idx2] = find(animalData{mouse}.tridIdx == signifIdx(i));
                cellPairIdx{mouse}{day}(:,i) = [idx1;idx2];
            end
            if ~isempty(cellPairIdx{mouse}{day}) % process the exception of bad days
                for i = 1:nNeuron
                    cellPairCount{mouse}(i,day) = ...
                    sum(cellPairIdx{mouse}{day}(1,:)==i) + sum(cellPairIdx{mouse}{day}(2,:)==i);
                end
                signifMatrix = sum(animalData{mouse}.tridIdx== reshape(signifIdx,1,1,[]),3);
            else
                cellPairCount{mouse}(:,day) = nan;
                signifMatrix = nan(nNeuron,nNeuron);
            end
            covMatSignif{mouse}(:,:,day) = signifMatrix; 
        end 
        t = toc; disp([ animalList{mouse} ' find significant pairs. Time = ' num2str(t,'%.2f') 'secs'])
    end
end

function [animalData] = removeBadSession(animalList, animalData)
    for mouse = 1:length(animalList)
        nDays = length(animalData{mouse}.baselineDay);
        for day = 1:nDays
            if ~usableSession(animalList{mouse},day)
                animalData{mouse}.covMat{day}(:,:) = nan;
                animalData{mouse}.allPair(:,day) = nan;
                animalData{mouse}.allPairZscore(:,day) = nan;
                animalData{mouse}.allPairH(:,day) = nan;     
            end
        end
    end
end

function useFlag = usableSession(mouse,day)
    useFlag = true;
    if strcmp(mouse,'cd017') && ((day==6)||(day==16))
        useFlag = false;
    end
    if strcmp(mouse,'cd037') && (day==12)
        useFlag = false;
    end
    if strcmp(mouse,'cd044') && (day==17)
        useFlag = false;
    end
end

function [shiftMat] = getAlignToDayMatrix(commonDays,selectDay)
    shiftMat = ones(length(commonDays), length(selectDay));
    upleftMat = ones(abs(selectDay(1)),abs(selectDay(1)));
    upleftMat = fliplr(1-triu(upleftMat));

    bottomrightMat = ones(abs(selectDay(end)),abs(selectDay(end)));
    bottomrightMat = fliplr(triu(bottomrightMat)');

    shiftMat(1:abs(selectDay(1)),1:abs(selectDay(1))) = upleftMat;
    shiftMat(end-abs(selectDay(end))+1:end, end-abs(selectDay(end))+1:end) = bottomrightMat;

end

function convertMat = convertTridMat(tridIdx, flatValues)
    convertMat = zeros(size(tridIdx));
    for i = 1:length(flatValues)
        convertMat(tridIdx==flatValues(i))=flatValues(i);
    end
    diagMat = diag(ones(1,size(convertMat,1)));
    convertMat(diagMat==1) = nan;
    convertMat = convertMat + convertMat';
end

