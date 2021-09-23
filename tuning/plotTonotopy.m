%% load data
clear; 
mouse = 'cd044';
switch mouse
    case 'cd017'
        preName = 'cd017_000_001';
        postName = 'cd017_016_007';
        target = 9514;
        foil = 16000;
    case 'cd036'
        preName = 'cd036_000_003';
        postName = 'cd036_017_001';
        target = 13454;
        foil = 22627;
    case 'cd037'
        preName = 'cd037_000_001';
        postName = 'cd037_017_003';
        target = 11314;
        foil = 6727;
    case 'cd042'
        preName = 'cd042_000_003';
        postName = 'cd042_016_001';
        target = 11314;
        foil = 6727;
    case 'cd044'
        preName = 'cd044_000_001';
        postName = 'cd044_016_001';
        target = 9514;
        foil = 16000;
    case 'cd019'
        preName = 'cd019_000_001';
        postName = 'cd019_016_007';
        target = 6727;
        foil = 4000;
        
end
preTuning = load(['D:\labData\excitatory\tuning\' preName '_Tuning\population\tuning.mat']);
postTuning = load(['D:\labData\excitatory\tuning\' postName '_Tuning\population\tuning.mat']);
savePath = 'D:\labData\excitatory\tuning\masterData\';
zscoreFlag = true;
mkdir([savePath mouse '\tonotopy']);

%%
pretoneFrames = 10;
rbmap = imresize(redbluecmap, [256, 3]);  % original color map contain just 11 colors, this increase it to 64
rbmap = min(max(rbmap, 0), 1);
nPlanes = 2;
nNeuron = size(preTuning.TCpretone_reorderCorr,4);

tones = [4000,4757,5657,6727,8000,...
    9514,11314,13454,16000,19027,22627,...
    26909,32000,38055,45255,53817,64000];
targIdx = find(tones==target);
foilIdx = find(tones==foil);


responsiveFlagPre= preTuning.anovaPeakCorr(2,:) > 0 & sum(preTuning.anovaSignifToneCorr,1) > 0;
responsiveFlagPost= postTuning.anovaPeakCorr(2,:) > 0 & sum(postTuning.anovaSignifToneCorr,1) > 0;
respCellFlag = responsiveFlagPre | responsiveFlagPost;

totalSel = sum(respCellFlag);
gainSel = (respCellFlag & ~responsiveFlagPre);
loseSel = (respCellFlag & ~responsiveFlagPost);
staySel = (responsiveFlagPre & responsiveFlagPost);

responsiveTonePre = preTuning.anovaSignifToneCorr;
responsiveTonePost = postTuning.anovaSignifToneCorr;

peakCountPreAll = sum(preTuning.tuningPeak(logical(responsiveFlagPre)) == (1:17)',2);
peakCountPostAll = sum(postTuning.tuningPeak(logical(responsiveFlagPost)) == (1:17)',2);
% figure;
% plot(peakCountPreAll,'LineWidth',3); hold on; plot(peakCountPostAll,'LineWidth',3); 
% ylimm = ylim; xlim([1 17])
% plot([targIdx targIdx],ylimm,'--g');plot([foilIdx foilIdx],ylimm,'--r');
% legend('pre','post');
% xlabel('tone'); ylabel('counts')


allToneIdx = 1:17;
for i = 1:nPlanes
    [tempCount, tempSignif] = getPeakCountPlane(preTuning.tuningPeak,preTuning.neuronPlane(i:i+1),responsiveFlagPre,allToneIdx);
    peakCountPre{i} = tempCount; peakSignifPre{i} = tempSignif;

    [tempCount, tempSignif] = getPeakCountPlane(postTuning.tuningPeak,postTuning.neuronPlane(i:i+1),responsiveFlagPost,allToneIdx);
    peakCountPost{i} = tempCount; peakSignifPost{i} = tempSignif;
end

neuronPlane = preTuning.neuronPlane;
img = preTuning.refImg;
xlen = size(preTuning.refImg{1},1);
ylen = size(preTuning.refImg{1},2);
rois = preTuning.roisBound;
roiCentroid = cell(size(rois));
for i = 1:length(rois)
    roiCentroid{i} = zeros(length(rois{i}),2);
    for j = 1:length(rois{i})
        temp = polyshape (rois{i}{j});
        [x,y] = centroid (temp);
        roiCentroid{i}(j,:) = [x,y]; 
    end
end
neuronEachPlane = preTuning.neuronEachPlane;

preTuning.TCpretone_zscoreCorr = zscoreTC(preTuning.TCpretone_reorderCorr);
postTuning.TCpretone_zscoreCorr = zscoreTC(postTuning.TCpretone_reorderCorr);
if zscoreFlag
    [peakActPre,~] = getPeakAct(preTuning.TCpretone_zscoreCorr,pretoneFrames,preTuning.peakFrames);
    [peakActPost,~] = getPeakAct(postTuning.TCpretone_zscoreCorr,pretoneFrames,postTuning.peakFrames);
else
    [peakActPre,~] = getPeakAct(preTuning.TCpretone_reorderCorr,pretoneFrames,preTuning.peakFrames);
    [peakActPost,~] = getPeakAct(postTuning.TCpretone_reorderCorr,pretoneFrames,postTuning.peakFrames);
    
end


%% get decoder weight, 5:5 cross-validation
tempT = squeeze(peakActPre(targIdx,:,:));
tempF = squeeze(peakActPre(foilIdx,:,:));
decoderPre = runDecoder(tempT, tempF,0.5);
[h,~] = ttest(decoderPre.cellAcc(:,:,2)',decoderPre.cellAccShuf(:,:,1)','tail','right');
signifFlagPre = (responsiveFlagPre & logical(h));

meanAccPre = mean(decoderPre.cellAcc(signifFlagPre,:,2),2);
cellWPre = mean(decoderPre.cellW(signifFlagPre,:),2);
[~,tempidx] = sort(abs(mean(cellWPre,2)),'descend');

%figure; plot(abs(mean(cellW(tempidx,:),2))); hold on; plot(std(cellW(tempidx,:),0,2));

tempT = squeeze(peakActPost(targIdx,:,:));
tempF = squeeze(peakActPost(foilIdx,:,:));
decoderPost = runDecoder(tempT, tempF,0.5);
[h,~] = ttest(decoderPost.cellAcc(:,:,2)',decoderPost.cellAccShuf(:,:,1)','tail','right');
signifFlagPost = (responsiveFlagPost & logical(h));

meanAccPost = mean(decoderPost.cellAcc(signifFlagPost,:,2),2);
cellWPost = mean(decoderPost.cellW(signifFlagPost,:),2);
decodeFlag = (signifFlagPre | signifFlagPost);

save([savePath '\' mouse '\' 'prePostTuning.mat']);
%% Make the plot
f = figure; nColumns = 5;  nRows = 3;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.1, 0.95, 0.9]);
margins = [0.03, 0.005];
planeToShow = 2;
for i = 1:nPlanes
    neuronPlaneIndex{i} = neuronPlane(i)+1:neuronPlane(i+1);
end
%% 1.1 plot tonotopic map 

h = {subplot_tight(nRows,nColumns,1,margins), subplot_tight(nRows,nColumns,nColumns+1,margins)};
colormapIndex = round(linspace(1,64,17));
func_tonotopyPlot(img{planeToShow}, rois{planeToShow},jet,colormapIndex,...
    {preTuning.tuningPeak, postTuning.tuningPeak},'axes', h,...
    'title', {'peak pre', 'peak post'}, ...
    'planeFlag', neuronPlaneIndex{planeToShow},...
    'plotFlag', {responsiveFlagPre, responsiveFlagPost});
%% 1.2 plot peak diff
h = {subplot_tight(nRows,nColumns,nColumns*2+1,margins)};
tuningChange = - preTuning.tuningPeak + postTuning.tuningPeak;
tempMax = prctile(tuningChange,95);
tempMin = prctile(tuningChange,5);
changeIdx = normalizeIndex({tuningChange}, tempMax, tempMin);

func_tonotopyPlot(img{planeToShow}, rois{planeToShow},rbmap,...
    1:256,changeIdx,'axes', h, 'planeFlag',...
    neuronPlaneIndex{planeToShow}, 'plotFlag', ...
    {staySel},'title', {'peak diff'});

%% 2.1 targ and foil activity in pre, post learning
h = {subplot_tight(nRows,nColumns,2,margins), subplot_tight(nRows,nColumns,nColumns+2,margins),...
    subplot_tight(nRows,nColumns,3,margins), subplot_tight(nRows,nColumns,nColumns+3,margins), ...
    subplot_tight(nRows,nColumns,nColumns*2+2,margins), subplot_tight(nRows,nColumns,nColumns*2+3,margins)};
targPre = squeeze(peakActPre(targIdx,:,:));
targPost = squeeze(peakActPost(targIdx,:,:));
targDiff = - targPre + targPost;
foilPre = squeeze(peakActPre(foilIdx,:,:));
foilPost = squeeze(peakActPost(foilIdx,:,:));
foilDiff = - foilPre + foilPost;

tempIdx = normalizeIndex({mean(targPre), mean(targPost),...
    mean(foilPre), mean(foilPost), mean(targDiff), mean(foilDiff)});

func_tonotopyPlot(img{planeToShow}, rois{planeToShow},rbmap,1:256,...
   tempIdx, 'axes', h, 'planeFlag',neuronPlaneIndex{planeToShow},...
   'title', {'targAct pre', 'targAct post','foilAct pre', 'foilAct post',...
   'targAct diff','foilAct diff'});
%% 2.2 foil activity in pre, post learning
% h = {subplot_tight(4,4,14), subplot_tight(4,4,15)};
% 
% func_tonotopyPlot(img{planeToShow}, rois{planeToShow},jet,colormapIndex,...
%    {preTuning.tuningPeak, preTuning.tuningPeak}, 'axes', h,...
%     'planeFlag',neuronPlaneIndex{planeToShow},'colorMod', {mean(targDiff), mean(foilDiff)});



%% plot4.1: barplot of significant neurons 
% barPlot({[sum(signifFlagPre) sum(signifFlagPost)], ...
%     [mean(mean(decoderPre.cellAcc(signifFlagPre,:,2),2),1),...
%     mean(mean(decoderPost.cellAcc(signifFlagPost,:,2),2),1)],...
%     [sum(cellWPre>0) sum(cellWPost>0);sum(cellWPre<0) sum(cellWPost<0)]},...
%     {'significant cell','significant cell',{'pre','post'}},...
%     {'# of decoding neuron','mean accuracy','by t/f preference'},...
%     {{'pre','post'},{'pre','post'},{'target-prefer','foil-prefer'}});
%% 3.1 decoder accuracy 
preIdx = mean(decoderPre.cellAcc(:,:,2),2) - 0.5;
postIdx = mean(decoderPost.cellAcc(:,:,2),2) - 0.5;

tempIdx = normalizeIndex({preIdx,postIdx, -preIdx + postIdx});

h = {subplot_tight(nRows,nColumns,4,margins), subplot_tight(nRows,nColumns,nColumns+4,margins),...
    subplot_tight(nRows,nColumns,nColumns*2+4,margins)};
func_tonotopyPlot(img{planeToShow}, rois{planeToShow},rbmap,...
    1:256,tempIdx,'axes', h, 'planeFlag',...
    neuronPlaneIndex{planeToShow}, 'plotFlag', ...
    {signifFlagPre,signifFlagPost,decodeFlag},...
    'title', {'decodeAcc pre','decodeAcc post','decodeAcc diff'});



%% 3.2 decoder weights
preIdx = zeros(1,nNeuron);preIdx(signifFlagPre) = cellWPre;
postIdx = zeros(1,nNeuron);postIdx(signifFlagPost) = cellWPost;

tempIdx = normalizeIndex({preIdx,postIdx, -preIdx + postIdx});

h = {subplot_tight(nRows,nColumns,5,margins), subplot_tight(nRows,nColumns,nColumns+5,margins),...
    subplot_tight(nRows,nColumns,nColumns*2+5,margins)};
func_tonotopyPlot(img{planeToShow}, rois{planeToShow},rbmap,...
    1:256,tempIdx,'axes', h, 'planeFlag',...
    neuronPlaneIndex{planeToShow}, 'plotFlag', ...
    {signifFlagPre,signifFlagPost,decodeFlag},...
    'title', {'decodeW pre','decodeW post','decodeW diff'});
saveas(gcf,[savePath '\' mouse '\tonotopy\f0_tonotopy_plane' int2str(planeToShow) '.png']);
%% all functions
function f1 = barPlot(data,legends, titles,xticklabelss)
    nData = length(data);
    f1 = figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.025, 0.3, 0.95, 0.4]);
    for i = 1:nData
       subplot(1,nData,i)
       bar(data{i});
       legend(legends{i},'Location','Best');
       title(titles{i})
       xticklabels(xticklabelss{i})
    end
end

function [peakAct,baseAct] = getPeakAct(TC,pretoneFrames,peakFrames)
nNeuron = size(TC,4);
peakAct = zeros(size(TC,2),size(TC,3),size(TC,4));
baseAct = zeros(size(TC,2),size(TC,3),size(TC,4));
for i = 1:nNeuron
    peakAct(:,:,i) = squeeze(TC(pretoneFrames+peakFrames(i),:,:,i));
    baseAct(:,:,1) = squeeze(TC(pretoneFrames,:,:,i));
end
end

function pretoneTC = zscoreTC(pretoneTC)
    nNeuron = size(pretoneTC,4);
    temp = reshape(pretoneTC,8500,nNeuron);
    temp = zscore(temp, 0,1);
    pretoneTC = reshape(temp,50,17,10,nNeuron);
end

function index1 = normalizeIndex(index1, tempMax, tempMin)
if nargin == 1
    tempMax = prctile(cell2mat(index1),95);tempMin = prctile(cell2mat(index1),5);
end
tempLim = max([abs(tempMax) abs(tempMin)]);
for i = 1:length(index1)
    tempIdx = index1{i}/tempLim;
    tempIdx = round(tempIdx*127.5 + 128.5);
    tempIdx(tempIdx > 256) = 256;
    tempIdx(tempIdx < 1) = 1;
    index1{i} = tempIdx;
end
end

function [peakCount, peakSignif] = getPeakCountPlane(peak,neuronPlane,responsiveFlag,allToneIdx)
    tempIdx = neuronPlane(1)+1:neuronPlane(2);
    temp = peak(tempIdx);  
    peakCount= sum(temp(responsiveFlag(tempIdx))==allToneIdx',2);
    peakSignif = find(peakCount > sum(responsiveFlag(tempIdx))*0.005);
end

function decoder = runDecoder(data1,data2, percTrain)
    nNeuron = size(data1,2);
    rep = 100;
    nTrial = size(data1,1);
    nTrain = round (percTrain * nTrial);
    nTest = nTrial - nTrain;
    decoder.cellW = zeros(nNeuron,rep);
    decoder.cellAcc = zeros(nNeuron,rep,2);
    decoder.popW = zeros(nNeuron,rep);
    decoder.popAcc = zeros(rep,2);

    decoder.cellWShuf = zeros(nNeuron,rep);
    decoder.cellAccShuf = zeros(nNeuron,rep,2);
    decoder.popWShuf = zeros(nNeuron,rep);
    decoder.popAccShuf = zeros(rep,2);
    for i = 1:rep
        idxT = randperm(nTrial);
        idxF = randperm(nTrial);
        trainT = data1(idxT(1:nTrain),:);
        trainF = data2(idxF(1:nTrain),:);
        testT = data1(idxT(nTrain+1:end),:);
        testF = data2(idxF(nTrain+1:end),:);

        
        [tempCellW, tempCellAcc, tempPopW, tempPopAcc] = getDecoderW(trainT,trainF,testT,testF);
        decoder.cellW(:,i) = tempCellW; decoder.cellAcc(:,i,:) = tempCellAcc;
        decoder.popW(:,i) = tempPopW; decoder.popAcc(i,:) = tempPopAcc;

        allData = [data1;data2];
        idxR = randperm(nTrial*2);
        trainT = allData(idxR(1:nTrain),:);
        trainF = allData(idxR(nTrain+1:nTrain*2),:);
        testT = allData(idxR(nTrain*2+1:nTrain*2+nTest),:);
        testF = allData(idxR(nTrain*2+nTest+1:end),:);
        [tempCellW, tempCellAcc, tempPopW, tempPopAcc] = getDecoderW(trainT,trainF,testT,testF);
        decoder.cellWShuf(:,i) = tempCellW; decoder.cellAccShuf(:,i,:) = tempCellAcc;
        decoder.popWShuf(:,i) = tempPopW; decoder.popAccShuf(i,:) = tempPopAcc;
    end
end


function [cellW, cellAcc, popW, popAcc] = getDecoderW(trainT,trainF,testT,testF)
    nTrain = size(trainT,1);
    nTest = size(testT,1);

    nNeuron = size(trainT,2);
    cellW = zeros(nNeuron,1);
    cellAcc = zeros(nNeuron,2);
    popAcc = zeros(1,2);
    for j = 1:nNeuron    
        tempW = (mean(trainT(:,j)) - mean(trainF(:,j))) / (var(trainT(:,j))/2 + var(trainF(:,j))/2);
        cellW(j) = tempW;
        predictT = tempW * trainT(:,j); predictF = tempW * trainF(:,j); 
        c = mean(predictT)/2 + mean(predictF)/2;
        correctT = predictT>c;correctF = predictF<c;
        cellAcc(j,1) = (sum(correctT) + sum(correctF)) / (2*nTrain);
        predictT = tempW * testT(:,j); correctT = predictT>c;
        predictF = tempW * testF(:,j); correctF = predictF<c;
        cellAcc(j,2) = (sum(correctT) + sum(correctF)) / (2*(nTest));
    end
    warning ('off')
    tempW = (mean(trainT) - mean(trainF)) * inv(cov(trainT)/2 + cov(trainF)/2);
    popW = tempW;
    predictT = trainT * tempW'; predictF = trainF * tempW'; 
    c = mean(predictT)/2 + mean(predictF)/2;
    correctT = predictT>c;correctF = predictF<c;
    popAcc(1) = (sum(correctT) + sum(correctF)) / (2*nTrain);
    predictT = testT * tempW'; predictF = testT * tempW'; 
    correctT = predictT>c;correctF = predictF<c;
    popAcc(2) = (sum(correctT) + sum(correctF)) / (2*(nTest));

end