function getBehaviorSI(mouse)
tic;
% load data
datapath = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' mouse];
savingPath = [datapath '\corrected_plane0\']; mkdir(savingPath);
load([datapath '\ar1_allday_fixBadSession\' mouse '_calman_ar1_foo90_pars_allday_s_plane0.mat'],'s');
configPath = 'C:\Users\zzhu34\Documents\tempdata\excitData\config\mouse';
filename = [configPath '\' mouse '_config.csv'];
configTable = readtable(filename);

% get the number of days
daysList = sort(unique(configTable.Day(strcmp(configTable.BehavType,'Behavior'))),'ascend');
daysList = reshape(daysList,1,[]);

% get SI of each day and test significance
peakDetection = false; peakFrames = 1:6;
for i = daysList
    sessionNum = reshape(find((configTable.Day == i) & strcmp(configTable.BehavType,'Behavior')),1,[]);
    
    
    
    tPeakMeanTrial = []; fPeakMeanTrial = []; tPeak = []; fPeak = [];
    for j = sessionNum

        tempS = s{j};tempS = [nan(size(tempS,1),1) tempS];
        behavMat =  importdata([configTable.behavpath{j} '\' configTable.BehavFile{j}]);
        [tempTPM, tempFPM, tempT, tempF,tempTPT,tempFPT] = getPeakAct(behavMat,tempS,peakDetection,peakFrames);       
        tPeakMeanTrial = cat(2,tPeakMeanTrial,tempTPT);
        fPeakMeanTrial = cat(2,fPeakMeanTrial,tempFPT);       
        tPeak = cat(3,tPeak,tempT);fPeak = cat(3,fPeak,tempF);
    end
    tAct{i} = tPeak; fAct{i} = fPeak;
    tPeakTrial{i} = tPeakMeanTrial;fPeakTrial{i} = fPeakMeanTrial;
    
    peakDetection = false;peakFrames = 1:6;smoothWindow = 0;
    for k = 1:size(tPeak,1)
        tempAct = squeeze(tPeak(k,:,:)); [h, tempAuc, tempIdx] = testResponsive(tempAct,peakDetection,peakFrames,smoothWindow);
        ttestT(k,i) = h;rocAucT(k,i) = tempAuc; 
        tempAct = squeeze(fPeak(k,:,:)); [h, tempAuc] = testResponsive(tempAct,peakDetection,peakFrames,smoothWindow);
        ttestF(k,i) = h;rocAucF(k,i) = tempAuc;
    end
    tPeakMeanTrial = mean(tPeakMeanTrial,2); fPeakMeanTrial = mean(fPeakMeanTrial,2);
    tempSI = (tPeakMeanTrial- fPeakMeanTrial) ./ (abs(tPeakMeanTrial) +  abs(fPeakMeanTrial));
    tempSI = reshape(tempSI,size(tempSI,1)*size(tempSI,2),[]);
    tempSI(isnan(tempSI)) = 0;
    SSI(:,i) = tempSI;
    SI(:,i) = abs(tempSI);
end
save([savingPath '\' mouse '_behaviorSI.mat'],'SI','SSI','ttestT','ttestF',...
    'rocAucT','rocAucF','tAct','fAct','tPeakTrial','fPeakTrial','daysList');
t = toc;
disp([mouse ' is finished. Time = ' num2str(t,'%.2f') ' sec'])
end


function [tPeakMean, fPeakMean, tAct, fAct,tPeakTrial,fPeakTrial] = getPeakAct(behavMatrix,act,peakDetection,peakFrames)
    tFrame = floor(behavMatrix(behavMatrix(:,4)==1 | behavMatrix(:,4)==2,12)/2);
    fFrame = floor(behavMatrix(behavMatrix(:,4)==3 | behavMatrix(:,4)==4,12)/2);
    selectFrame = -5 : 30; toneFrame = abs(selectFrame(1)); selectToneFrame = toneFrame + peakFrames;
    tAct = []; fAct = [];
    for k = 1:length(tFrame);tAct(:,:,k) = act(:,tFrame(k) + selectFrame);end
    for k = 1:length(fFrame);fAct(:,:,k) = act(:,fFrame(k) + selectFrame);end
    if peakDetection
        tActMean = mean(tAct,3); tPeakMean = max(smoothdata(tActMean(:,selectToneFrame),2,'gaussian',3),[],2);
        fActMean = mean(fAct,3); fPeakMean = max(smoothdata(fActMean(:,selectToneFrame),2,'gaussian',3),[],2);
        tPeakTrial = squeeze(max(smoothdata(tAct(:,selectToneFrame,:),2,'gaussian',3),[],2));
        fPeakTrial = squeeze(max(smoothdata(fAct(:,selectToneFrame,:),2,'gaussian',3),[],2));
        
    else
        tPeakTrial = squeeze(mean(tAct(:,selectToneFrame,:),2));tPeakMean = mean(tPeakTrial,2);
        fPeakTrial = squeeze(mean(fAct(:,selectToneFrame,:),2));fPeakMean = mean(fPeakTrial,2);
    end
end

function [h, tempAuc, peakIdx] = testResponsive(tempAct,peakDetection,peakFrames,smoothWindow)
    tempBaseline = tempAct(1:3,:); tempBaselineMean = mean(tempBaseline,1);
    if smoothWindow ~= 0; tempMeanAct = smoothdata(mean(tempAct,2),'gaussian',smoothWindow);
    else; tempMeanAct = mean(tempAct,2); end
    toneFrame = 5; % number of pretone frames
    selectToneFrame = toneFrame + peakFrames;
    [~,tempPeakIdx] = max(tempMeanAct(selectToneFrame));
    peakIdx = tempPeakIdx + selectToneFrame(1)-1; 
    if peakDetection; tempPeak = tempAct(peakIdx,:); 
    else; tempPeak = mean(tempAct(selectToneFrame,:),1); end
    [h,p] = ttest(tempPeak,tempBaselineMean,'tail','right');
    % do an roc for each tone 
    rocAct = [tempBaseline(:)' tempPeak]; 
    rocAct(rocAct<1e-4) = randn(1,sum(rocAct<1e-4))*1e-4;% introduce small noise to spks to avoid domination of 0s
    rocLabel = [zeros(1,numel(tempBaseline)) ones(1,length(tempPeak))];
    [tpr, fpr, threshold] = roc(rocLabel, rocAct);
    tempAuc = trapz([0 fpr 1],[0 tpr 1]);
end