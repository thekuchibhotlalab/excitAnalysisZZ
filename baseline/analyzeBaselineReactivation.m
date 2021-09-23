%% load the data
clear;
experimental = {'cd017','cd036','cd037','cd042','cd044'};
for i = 1:length(experimental)
    expDatapath{i} = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' experimental{i} '\corrected_plane0'];
    tempData = load([expDatapath{i} '\' experimental{i} '_baseline_data_pairwiseCorr.mat']);
    expData{i} = tempData;
    tempData = load([expDatapath{i} '\' experimental{i} '_behaviorSI.mat']);
    expDataSI{i} = tempData;

    tempData = load([expDatapath{i} '\' experimental{i} '_baselineS.mat']);
    expDataS{i} = tempData.baselineS;
    
    if strcmp(experimental(i),'cd037')
        load('C:\Users\zzhu34\Documents\tempdata\excitData\cd037\TC\cd037_nanflag_plane0.mat','nanflag');
        expDataS{i}(nanflag,:,:) = [];
        for j = 1:length(expDataSI{i}.tPeakTrial)
            expDataSI{i}.tPeakTrial{j}(nanflag,:) = [];
            expDataSI{i}.fPeakTrial{j}(nanflag,:) = [];
        end
    end
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
    tempData = load([ctrlDatapath{i} '\' ctrl{i} '_baselineS.mat']);
    ctrlDataS{i} = tempData.baselineS;

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
%%
expDecoder = {};expReact = {}; 
expReactMat = zeros(length(experimental),length(commonDays),size(expDataS{1},2));
expReactMatZ = zeros(length(experimental),length(commonDays),size(expDataS{1},2));
expReactMatZShuf = zeros(length(experimental),length(commonDays),size(expDataS{1},2));
shuffleRep = 20; 
for mouse = 1:length(experimental)
    tic; disp([experimental{mouse} 'started.'])
    for day = 1:length(commonDays)
        
        expDecoder{mouse, day} = getDecoderWeight(expDataSI{mouse}.tPeakTrial{day},expDataSI{mouse}.fPeakTrial{day});
        expReact{mouse, day} = runDecorder(expDecoder{mouse, day},expDataS{mouse}(:,:,day));
         
        % shuffle to determine the mean and std of the baseline
        shuffleS = shuffle(expDataS{mouse}(:,:,day), 2,shuffleRep);
        shuffleS = reshape(shuffleS,size(shuffleS,1),[]);    
        tempResults = runDecorder(expDecoder{mouse, day},shuffleS);
        tempRep = tempResults.proj;
        shuffleMean = nanmean(tempRep); shuffleStd = std(tempRep);
        % zscore the decoding projection value
        tempProj = expReact{mouse, day}.proj;
        expReactMat(mouse,day,:) = tempProj;
        expReactMatZ(mouse,day,:) = (tempProj-nanmedian(tempProj)) ./ std(tempProj);
        expReactMatZShuf(mouse,day,:) = (tempProj - shuffleMean) / shuffleStd;
    end
    t = toc;
    disp([experimental{mouse} ' is finished. Time = ' num2str(t,'%.2f') ' sec'])
end
%
ctrlDecoder = {};ctrlReact = {}; 
ctrlReactMat = zeros(length(ctrl),length(commonDays),size(ctrlDataS{1},2));
ctrlReactMatZ = zeros(length(ctrl),length(commonDays),size(ctrlDataS{1},2));
ctrlReactMatZShuf = zeros(length(ctrl),length(commonDays),size(ctrlDataS{1},2));
shuffleRep = 20; 
for mouse = 1:length(ctrl)
    tic; disp([ctrl{mouse} 'started.'])
    for day = 1:length(commonDays)
        
        ctrlDecoder{mouse, day} = getDecoderWeight(ctrlDataSI{mouse}.tPeakTrial{day},ctrlDataSI{mouse}.fPeakTrial{day});
        ctrlReact{mouse, day} = runDecorder(ctrlDecoder{mouse, day},ctrlDataS{mouse}(:,:,day));
         
        % shuffle to determine the mean and std of the baseline
        shuffleS = shuffle(ctrlDataS{mouse}(:,:,day), 2,shuffleRep);
        shuffleS = reshape(shuffleS,size(shuffleS,1),[]);    
        tempResults = runDecorder(ctrlDecoder{mouse, day},shuffleS);
        tempRep = tempResults.proj;
        shuffleMean = nanmean(tempRep); shuffleStd = std(tempRep);
        % zscore the decoding projection value
        tempProj = ctrlReact{mouse, day}.proj;
        ctrlReactMat(mouse,day,:) = tempProj;
        ctrlReactMatZ(mouse,day,:) = (tempProj-nanmedian(tempProj)) ./ std(tempProj);
        ctrlReactMatZShuf(mouse,day,:) = (tempProj - shuffleMean) / shuffleStd;
    end
    t = toc;
    disp([ctrl{mouse} ' is finished. Time = ' num2str(t,'%.2f') ' sec'])
end
%% select significant reactivation based on zscore value
threshold = 5;

expReactT = expReactMatZShuf > threshold; expReactT_on = diff(expReactT,1,3)==1;
expReactTSum = nansum(expReactT_on,3); 
expReactF = expReactMatZShuf < -threshold; expReactF_on = diff(expReactF,1,3)==1;
expReactFSum = nansum(expReactF_on,3);
ctrlReactT = ctrlReactMatZShuf > threshold; ctrlReactT_on =  diff(ctrlReactT,1,3)==1;
ctrlReactTSum = nansum(ctrlReactT_on,3); ctrlReactTSum(1,1) = nan;
ctrlReactF = ctrlReactMatZShuf < -threshold; ctrlReactF_on =  diff(ctrlReactF,1,3)==1;
ctrlReactFSum = nansum(ctrlReactF_on,3); ctrlReactFSum(1,1) = nan;

ylimm = [0 180];
figure; subplot(1,2,1); hold on
[hT, ~] = fn_plotMeanSampleLine(commonDays',expReactTSum,{'Color',matlabColors(1),'LineWidth',1.5},...
    {'Color',matlabColors(1,0.3),'LineWidth',0.5});
[hF, ~] = fn_plotMeanSampleLine(commonDays',expReactFSum,{'Color',matlabColors(2),'LineWidth',1.5},...
    {'Color',matlabColors(2,0.3),'LineWidth',0.5});
legend([hT hF],'target','foil','Location','Best'); xlim([1 15]); ylim(ylimm)
xlabel('days'); ylabel('# reactivation event');title('Experimental');
subplot(1,2,2); hold on 
[hT, ~] = fn_plotMeanSampleLine(commonDays',ctrlReactTSum,{'Color',matlabColors(1),'LineWidth',1.5},...
    {'Color',matlabColors(1,0.3),'LineWidth',0.5});
[hF, ~] = fn_plotMeanSampleLine(commonDays',ctrlReactFSum,{'Color',matlabColors(2),'LineWidth',1.5},...
    {'Color',matlabColors(2,0.3),'LineWidth',0.5});
legend([hT hF],'target','foil','Location','Best'); xlim([1 15]); ylim(ylimm)
xlabel('days'); ylabel('# reactivation event');title('Control');
%fn_plotMeanSampleLine(commonDays',ctrlReactF(2,:),{'Color',matlabColors(3),'LineWidth',1.5},...
%    {'Color',matlabColors(3,0.3),'LineWidth',0.5});
%% separate out animal groups and plot
ylimm = [0 200];
figure;
subplot(3,2,1)
[hT, ~] = fn_plotMeanSampleLine(commonDays',expReactTSum([1 2 5],:),{'Color',matlabColors(1),'LineWidth',1.5},...
    {'Color',matlabColors(1,0.3),'LineWidth',0.5});
[hT, ~] = fn_plotMeanSampleLine(commonDays',expReactFSum([1 2 5],:),{'Color',matlabColors(2),'LineWidth',1.5},...
    {'Color',matlabColors(2,0.3),'LineWidth',0.5});
title('experimental');xlim([1 15]); ylim(ylimm);xlabel('days'); ylabel('# reactivation event');
subplot(3,2,3)
[hT, ~] = fn_plotMeanSampleLine(commonDays',expReactTSum(3,:),{'Color',matlabColors(1),'LineWidth',1.5},...
    {'Color',matlabColors(1,0.3),'LineWidth',0.5});
[hT, ~] = fn_plotMeanSampleLine(commonDays',expReactFSum(3,:),{'Color',matlabColors(2),'LineWidth',1.5},...
    {'Color',matlabColors(2,0.3),'LineWidth',0.5});
xlim([1 15]); ylim(ylimm);xlabel('days'); ylabel('# reactivation event');
subplot(3,2,5)
[hT, ~] = fn_plotMeanSampleLine(commonDays',expReactTSum(4,:),{'Color',matlabColors(1),'LineWidth',1.5},...
    {'Color',matlabColors(1,0.3),'LineWidth',0.5});
[hT, ~] = fn_plotMeanSampleLine(commonDays',expReactFSum(4,:),{'Color',matlabColors(2),'LineWidth',1.5},...
    {'Color',matlabColors(2,0.3),'LineWidth',0.5});
xlim([1 15]); ylim(ylimm);xlabel('days'); ylabel('# reactivation event');
subplot(3,2,2)
fn_plotMeanSampleLine(commonDays',ctrlReactTSum(1,:),{'Color',matlabColors(1),'LineWidth',1.5},...
    {'Color',matlabColors(1,0.3),'LineWidth',0.5});
fn_plotMeanSampleLine(commonDays',ctrlReactFSum(1,:),{'Color',matlabColors(2),'LineWidth',1.5},...
    {'Color',matlabColors(2,0.3),'LineWidth',0.5});
xlim([1 15]); ylim(ylimm);
xlabel('days'); ylabel('# reactivation event');title('Control');
subplot(3,2,4)
fn_plotMeanSampleLine(commonDays',ctrlReactTSum(2,:),{'Color',matlabColors(1),'LineWidth',1.5},...
    {'Color',matlabColors(1,0.3),'LineWidth',0.5});
fn_plotMeanSampleLine(commonDays',ctrlReactFSum(2,:),{'Color',matlabColors(2),'LineWidth',1.5},...
    {'Color',matlabColors(2,0.3),'LineWidth',0.5});
xlim([1 15]); ylim(ylimm);xlabel('days'); ylabel('# reactivation event');
%% plot example session
animal = 1; day = 11; xlimm = [1 500];
figure; 
subplot(3,1,1);
imagesc(expDataS{animal}(:,:,day)); caxis([0 0.1])
tempT = expDataS{animal}(:,:,day); tempF = expDataS{animal}(:,:,day); 
tFlag = squeeze(expReactT(animal,day,:)); 
fFlag = squeeze(expReactF(animal,day,:)); 
tempT(:,~tFlag) = 0; tempF(:,~fFlag) = 0;
xlim(xlimm);
subplot(3,1,2);
imagesc(tempT); caxis([0 0.1]);xlim(xlimm);
subplot(3,1,3);
imagesc(tempF); caxis([0 0.1]);xlim(xlimm);
%plot(squeeze(expReactT(1,1,:)))
%%
function decoder = getDecoderWeight(tActTrial,fActTrial,xvalidT,xvalidF)
    
    decoder.nNeuron = size(tActTrial,1);
    decoder.meanT = nanmean(tActTrial,2); decoder.meanF = nanmean(fActTrial,2);
    
    covT = cov(tActTrial'); covF = cov(fActTrial');
    varT = diag(diag(covT)); varF = diag(diag(covF)); 
    decoder.varMean = (varT + varF) / 2;

    decoder.w = (decoder.varMean) \ (decoder.meanT - decoder.meanF);
    decoder.c = decoder.w' * (decoder.meanT + decoder.meanF) / 2; 

    results = runDecorder(decoder, tActTrial,fActTrial);
    decoder = catstruct(decoder,results);

    if exist('xvalidT')
        xvalid = runDecorder(decoder, xvalidT,xvalidF);decoder.xvalid = xvalid;
    end
    
end

function results = runDecorder(decoder, T,F)
    if ~isempty(decoder.w)
        if exist('F')
            results.projT = decoder.w' * T - decoder.c; results.predT = results.projT > 0;
            results.projF = decoder.w' * F - decoder.c; results.predF = results.projF <= 0;
            results.acc = nanmean([results.predT results.predF]);
        else
            results.proj = decoder.w' * T - decoder.c; 
            results.predT = results.proj > 0; results.predF = results.proj <= 0;
        end
    else
        if exist('F')
            results.projT = nan(1,size(T,2)); results.projF = nan(1,size(F,2));
            results.acc = nan;
        else
            results.proj = nan(1,size(T,2));
            results.predT = nan(1,size(T,2));
        end
        
    end
end

function [shuffleData] = shuffle(data, dim,nRep)
    datadim = ndims(data);
    nSamp = size(data,dim);
    shuffleData = repmat(data,[ones(1,datadim), nRep]);
    if ~exist('dim'); dim = 1; end
    for i = 1:nRep
        
        if dim==1
            for j = 1:size(data,2)
                randIdx = randperm(nSamp);
                shuffleData(:,j,i) = data(randIdx,j);
            end
        elseif dim==2
            for j = 1:size(data,1)
                randIdx = randperm(nSamp);
                shuffleData(j,:,i) = data(j,randIdx);
            end
        end
        
    end
    shuffleData = squeeze(shuffleData);
end