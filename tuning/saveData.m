clear; 
mouse = 'cd017';
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
colormapIndex = round(linspace(1,64,17));
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

totalSel = sum(respCellFlag);
gainSel = (respCellFlag & ~responsiveFlagPre);
loseSel = (respCellFlag & ~responsiveFlagPost);
staySel = (responsiveFlagPre & responsiveFlagPost);

%% classify neurons based on initial peak location
tuningChange = preTuning.tuningPeak - postTuning.tuningPeak;
tPre = (targIdx == preTuning.tuningPeak) & responsiveFlagPre;
fPre = (foilIdx == preTuning.tuningPeak) & responsiveFlagPre;

tPost = (targIdx == postTuning.tuningPeak) & responsiveFlagPost;
fPost = (foilIdx == postTuning.tuningPeak) & responsiveFlagPost;

if targIdx < foilIdx
    middlePre = (targIdx < preTuning.tuningPeak...
        & foilIdx > preTuning.tuningPeak & responsiveFlagPre);
    tSidePre = (targIdx > preTuning.tuningPeak...
        & foilIdx > preTuning.tuningPeak & responsiveFlagPre);
    fSidePre = (targIdx < preTuning.tuningPeak...
        & foilIdx < preTuning.tuningPeak & responsiveFlagPre);
else
    middlePre = (targIdx > preTuning.tuningPeak...
        & foilIdx < preTuning.tuningPeak & responsiveFlagPre);
    tSidePre = (targIdx < preTuning.tuningPeak...
        & foilIdx < preTuning.tuningPeak & responsiveFlagPre);
    fSidePre = (targIdx > preTuning.tuningPeak...
        & foilIdx > preTuning.tuningPeak & responsiveFlagPre);
end

peak.pre = preTuning.tuningPeak;
peak.post = postTuning.tuningPeak;
peak.tuningChange = tuningChange;
peak.tuningUp = sum(tuningChange(staySel)<0);
peak.tuningDown = sum(tuningChange(staySel)>0);
peak.tuningStay = sum(tuningChange(staySel)==0);

peak.dist.tPreDist = abs(targIdx - preTuning.tuningPeak);
peak.dist.fPreDist = abs(foilIdx - preTuning.tuningPeak);
peak.dist.tPostDist = abs(targIdx - postTuning.tuningPeak);
peak.dist.fPostDist = abs(foilIdx - postTuning.tuningPeak);

if targIdx > foilIdx
    peak.tPre.stay = sum(tuningChange(tPre & staySel)==0);
    peak.tPre.switch = sum(tuningChange(tPre & staySel)>0);
    peak.tPre.away = sum(tuningChange(tPre & staySel)<0);

    peak.fPre.stay = sum(tuningChange(fPre & staySel)==0);
    peak.fPre.switch = sum(tuningChange(fPre & staySel)<0);
    peak.fPre.away = sum(tuningChange(fPre & staySel)>0);

    peak.middlePre.stay = sum(tuningChange(middlePre & staySel)==0);
    peak.middlePre.T = sum(tuningChange(middlePre & staySel)<0);
    peak.middlePre.F = sum(tuningChange(middlePre & staySel)>0);

    peak.tSidePre.stay = sum(tuningChange(tSidePre & staySel)<=0);
    peak.tSidePre.closeT = sum(tuningChange(tSidePre & staySel)>0 & (postTuning.tuningPeak(tSidePre & staySel) > mean([targIdx foilIdx]) )  );
    peak.tSidePre.closeF = sum(tuningChange(tSidePre & staySel)>0 & (postTuning.tuningPeak(tSidePre & staySel) < mean([targIdx foilIdx]) )  );

    peak.fSidePre.stay = sum(tuningChange(fSidePre & staySel)>=0);
    peak.fSidePre.closeT = sum(tuningChange(fSidePre & staySel)<0 & (postTuning.tuningPeak(fSidePre & staySel) > mean([targIdx foilIdx]) )  );
    peak.fSidePre.closeF = sum(tuningChange(fSidePre & staySel)<0 & (postTuning.tuningPeak(fSidePre & staySel) < mean([targIdx foilIdx]) )  );

    peak.noPre.closeT = sum(gainSel & (postTuning.tuningPeak > mean([targIdx foilIdx]) )  );
    peak.noPre.closeF = sum(gainSel & (postTuning.tuningPeak < mean([targIdx foilIdx]) )  );

else
    peak.tPre.stay = sum(tuningChange(tPre & staySel)==0);
    peak.tPre.switch = sum(tuningChange(tPre & staySel)<0);
    peak.tPre.away = sum(tuningChange(tPre & staySel)>0);

    peak.fPre.stay = sum(tuningChange(fPre & staySel)==0);
    peak.fPre.switch = sum(tuningChange(fPre & staySel)>0);
    peak.fPre.away = sum(tuningChange(fPre & staySel)<0);

    peak.middlePre.stay = sum(tuningChange(middlePre & staySel)==0);
    peak.middlePre.T = sum(tuningChange(middlePre & staySel)>0);
    peak.middlePre.F = sum(tuningChange(middlePre & staySel)<0);

    peak.tSidePre.stay = sum(tuningChange(tSidePre & staySel)>=0);
    peak.tSidePre.closeT = sum(tuningChange(tSidePre & staySel)<0 & (postTuning.tuningPeak(tSidePre & staySel) < mean([targIdx foilIdx]) )  );
    peak.tSidePre.closeF = sum(tuningChange(tSidePre & staySel)<0 & (postTuning.tuningPeak(tSidePre & staySel) > mean([targIdx foilIdx]) )  );

    peak.fSidePre.stay = sum(tuningChange(fSidePre & staySel)<=0);
    peak.fSidePre.closeT = sum(tuningChange(fSidePre & staySel)>0 & (postTuning.tuningPeak(fSidePre & staySel) < mean([targIdx foilIdx]) )  );
    peak.fSidePre.closeF = sum(tuningChange(fSidePre & staySel)>0 & (postTuning.tuningPeak(fSidePre & staySel) > mean([targIdx foilIdx]) )  );

    peak.noPre.closeT = sum(gainSel & (postTuning.tuningPeak < mean([targIdx foilIdx]) )  );
    peak.noPre.closeF = sum(gainSel & (postTuning.tuningPeak > mean([targIdx foilIdx]) )  );

end
peak.tPre.lose = sum(tPre & loseSel);
peak.fPre.lose = sum(fPre & loseSel);
peak.middlePre.lose = sum(middlePre & loseSel);
peak.tSidePre.lose = sum(tSidePre & loseSel);
peak.fSidePre.lose = sum(fSidePre & loseSel);

peak.nResp = sum(respCellFlag);


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

%% get decoding weight
prepostTuning.decode.cellWPre = cellWPre;
prepostTuning.decode.cellWPost = cellWPost;

totalDecode = (signifFlagPre|signifFlagPost);
gainDecode = totalDecode & ~signifFlagPre;
loseDecode = totalDecode & ~signifFlagPost;
stayDecode = (signifFlagPre&signifFlagPost);

cellWPre = mean(decoderPre.cellW(decodeFlag,:),2);
cellWPost = mean(decoderPost.cellW(decodeFlag,:),2);

prepostTuning.decode.allCellWPre = cellWPre;
prepostTuning.decode.allCellWPost = cellWPost;

cellPrePostAbs = abs(cellWPre) - abs(cellWPost);
cellPrePost = cellWPre - cellWPost;

tFlagPre = mean(decoderPre.cellW(:,:),2) > 0 & signifFlagPre'; tFlagPre = tFlagPre(decodeFlag);
fFlagPre = mean(decoderPre.cellW(:,:),2) < 0 & signifFlagPre'; fFlagPre = fFlagPre(decodeFlag);
tFlagPost = mean(decoderPost.cellW(:,:),2) > 0 & signifFlagPost'; tFlagPost = tFlagPost(decodeFlag);
fFlagPost = mean(decoderPost.cellW(:,:),2) < 0 & signifFlagPost'; fFlagPost = fFlagPost(decodeFlag);

decoder.prePost.tFlagPre = tFlagPre;
decoder.prePost.fFlagPre = fFlagPre;
decoder.prePost.tFlagPost = tFlagPost;
decoder.prePost.fFlagPost = fFlagPost;
%% quantify change in decoder weight
rep = 100;
shufData = [];
shufDataPre = [];
shufDataAbs = [];
for i = 1:1000
    randidx1 = randperm(sum(respCellFlag)-sum(decodeFlag));
    temp = mean(decoderPre.cellW(~decodeFlag,:),2); tempPre = temp(randidx1);
    randidx2 = randperm(sum(respCellFlag)-sum(decodeFlag));
    temp = mean(decoderPost.cellW(~decodeFlag,:),2); tempPost = temp(randidx2);
    shufData = [shufData;tempPre-tempPost];
    shufDataAbs = [shufDataAbs;abs(tempPre)-abs(tempPost)];
    shufDataPre = [shufDataPre;tempPre];
end
downFlag = cellPrePost > mean(shufData)+std(shufData)*2;
upFlag = cellPrePost < mean(shufData)-std(shufData)*2;

tChange.up = sum(tFlagPre&upFlag);
tChange.down = sum(tFlagPre&downFlag);
tChange.stay = sum(tFlagPre&(~upFlag & ~downFlag));
tChange.downSwitch = sum(tFlagPre&downFlag&fFlagPost);
decoder.tChange = tChange;

fChange.up = sum(fFlagPre&upFlag);
fChange.down = sum(fFlagPre&downFlag);
fChange.stay = sum(fFlagPre&(~upFlag & ~downFlag));
fChange.upSwitch = sum(fFlagPre&upFlag&tFlagPost);
decoder.fChange = fChange;

nChange.up = sum((~tFlagPre & ~fFlagPre)&upFlag);
nChange.down =  sum((~tFlagPre & ~fFlagPre)&downFlag);
nChange.stay = sum((~tFlagPre & ~fFlagPre)&(~upFlag & ~downFlag));

nChange.upT =  sum((~tFlagPre & ~fFlagPre)&upFlag&tFlagPost);
nChange.upF = sum((~tFlagPre & ~fFlagPre)&upFlag&fFlagPost);
nChange.downT = sum((~tFlagPre & ~fFlagPre)&downFlag&tFlagPost);
nChange.downF = sum((~tFlagPre & ~fFlagPre)&downFlag&fFlagPost);
decoder.nChange = nChange;
decoder.nDecode = sum(decodeFlag);

%% plot4.9: change in absolute decoder weight based on their initial decoding weight

downFlagAbs = cellPrePostAbs > mean(shufDataAbs)+std(shufDataAbs)*2;
upFlagAbs = cellPrePostAbs < mean(shufDataAbs)-std(shufDataAbs)*2;

tChangeAbs.up = sum(tFlagPre&upFlagAbs);
tChangeAbs.down = sum(tFlagPre&downFlagAbs);
tChangeAbs.stay = sum(tFlagPre&(~upFlagAbs & ~downFlagAbs));
decoder.tChangeAbs = tChangeAbs;

fChangeAbs.up = sum(fFlagPre&upFlagAbs);
fChangeAbs.down = sum(fFlagPre&downFlagAbs);
fChangeAbs.stay = sum(fFlagPre&(~upFlagAbs & ~downFlagAbs));
decoder.fChangeAbs = fChangeAbs;

nChangeAbs.up = sum((~tFlagPre & ~fFlagPre)&upFlagAbs);
nChangeAbs.down =  sum((~tFlagPre & ~fFlagPre)&downFlagAbs);
nChangeAbs.stay = sum((~tFlagPre & ~fFlagPre)&(~upFlagAbs & ~downFlagAbs));

nChangeAbs.upT =  sum((~tFlagPre & ~fFlagPre)&upFlagAbs&tFlagPost);
nChangeAbs.upF = sum((~tFlagPre & ~fFlagPre)&upFlagAbs&fFlagPost);
decoder.nChangeAbs = nChangeAbs;
%% draw line

f1 = plotMapIdx(img, rois, neuronEachPlane,jet,colormapIndex,...
   {preTuning.tuningPeak, postTuning.tuningPeak}, ...
    {responsiveFlagPre, responsiveFlagPost});
h = subplot(2,2,3);
disp('waiting to draw line')
lineRoi = drawline(h);
pause();disp('waiting key press')
linePos= lineRoi.Position;
close(f1);
plotMapIdxLine(img, rois, neuronEachPlane,jet,colormapIndex,...
   {preTuning.tuningPeak, postTuning.tuningPeak}, ...
    {responsiveFlagPre, responsiveFlagPost}, linePos);
projPre = projectCoordPlaneTone(preTuning.tuningPeak, responsiveFlagPre, roiCentroid,linePos,peakSignifPre, neuronPlane);
projPost = projectCoordPlaneTone(postTuning.tuningPeak, responsiveFlagPre, roiCentroid,linePos,peakSignifPost, neuronPlane);
%% save data
close all;
save([savePath '\' mouse '\' 'prePostTuning.mat']);

%% functions

function [f1] = plotMapIdx(img, rois, neuronEachPlane,C,colormapIndex,plotIndex, plotFlag)
    f1 = figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.04, 0.6, 0.9]);
    if ~iscell(plotIndex)
        plotIndex = {plotIndex};
    end
    if ~iscell(plotFlag)
        plotFlag = {plotFlag};
    end
    nPlanes = length(neuronEachPlane);
    for i = 1:nPlanes
        neuronPlane = [0 cumsum(neuronEachPlane)];
        for k = 1:length(plotIndex)
            subplot(length(plotIndex),2,(k-1)*2+i);             
            imagesc(img{i});colormap gray;hold on;
            ylim([0 size(img{i},1)]);xlim([0 size(img{1},2)]);
            for j = 1:neuronEachPlane(i)
                cellIndex = neuronPlane(i) + j;
                x = rois{i}{j}(:,1); %freehand rois have the outlines in x-y coordinates
                y = rois{i}{j}(:,2); %matlab matrices are inverted so x values are the 2nd column of mn coordinates, and y is the 1st columna
                if plotFlag{k}(cellIndex)
                    %plot(x,y,'.','color',C(colormapIndex(peakIndex(cellIndex)),:),'MarkerSize',1);
                    patch(x,y,C(colormapIndex(plotIndex{k}(cellIndex)),:),'EdgeColor','none');
                else
                    patch(x,y,[0.8 0.8 0.8],'EdgeColor','none');
                    %plot(x,y,'.','color',[0.8 0.8 0.8],'MarkerSize',1);
                end
            end
            title(['Plane' int2str(i)])
            axis off;
        end
    end
end


function [f1] = plotMapIdxLine(img, rois, neuronEachPlane,C,colormapIndex,plotIndex, plotFlag, linePos)
    f1 = figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.04, 0.6, 0.9]);
    if ~iscell(plotIndex)
        plotIndex = {plotIndex};
    end
    if ~iscell(plotFlag)
        plotFlag = {plotFlag};
    end
    nPlanes = length(neuronEachPlane);
    for i = 1:nPlanes
        neuronPlane = [0 cumsum(neuronEachPlane)];
        for k = 1:length(plotIndex)
            h = subplot(length(plotIndex),2,(k-1)*2+i);             
            imagesc(ones(size(img{i})));colormap gray;hold on;
            ylim([0 size(img{i},1)]);xlim([0 size(img{1},2)]);
            for j = 1:neuronEachPlane(i)
                cellIndex = neuronPlane(i) + j;
                x = rois{i}{j}(:,1); %freehand rois have the outlines in x-y coordinates
                y = rois{i}{j}(:,2); %matlab matrices are inverted so x values are the 2nd column of mn coordinates, and y is the 1st columna
                if plotFlag{k}(cellIndex)
                    %plot(x,y,'.','color',C(colormapIndex(peakIndex(cellIndex)),:),'MarkerSize',1);
                    patch(x,y,C(colormapIndex(plotIndex{k}(cellIndex)),:),'EdgeColor','none');
                else
                    patch(x,y,[0.8 0.8 0.8],'EdgeColor','none');
                    %plot(x,y,'.','color',[0.8 0.8 0.8],'MarkerSize',1);
                end
            end
            plot(linePos(:,1),linePos(:,2),'LineWidth',2,'Color',[0 0 0]); 
            set(h, 'Ydir', 'reverse')
            title(['Plane' int2str(i)])
            axis off;
        end
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

tempLim = max([abs(tempMax) abs(tempMin)]);
index1 = index1/tempLim;
index1 = round(index1*127.5 + 128.5);
index1(index1 > 256) = 256;
index1(index1 < 1) = 1;

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

function [peakCount, peakSignif] = getPeakCountPlane(peak,neuronPlane,responsiveFlag,allToneIdx)
    tempIdx = neuronPlane(1)+1:neuronPlane(2);
    temp = peak(tempIdx);  
    peakCount= sum(temp(responsiveFlag(tempIdx))==allToneIdx',2);
    peakSignif = find(peakCount > sum(responsiveFlag(tempIdx))*0.005);
end

function proj = projectCoord(roiCentroid,linePos)
    nNeuron = size(roiCentroid,1);
    lineSlope = (linePos(2,2)-linePos(1,2)) / (linePos(2,1)-linePos(1,1));
    tempVec = [1; lineSlope]; tempVec = tempVec / norm(tempVec);
    proj = roiCentroid * tempVec;
end

function proj = projectCoordPlaneTone(peak, responsiveFlag, roiCentroid,linePos,allToneIdx, neuronPlane)
    nPlanes = length(roiCentroid);
    for i = 1:nPlanes
        nTones = length(allToneIdx{i});
        tempIdx = neuronPlane(i)+1:neuronPlane(i+1);
        tempPeak = peak(tempIdx); 
        tempCent = roiCentroid{i};
        for j = 1:nTones
            cellIdx = (tempPeak == allToneIdx{i}(j));
            proj{i}{j} = projectCoord( tempCent(cellIdx & responsiveFlag(tempIdx),:), linePos);
        end
    end
end

function [toneStat,allCounts] = plotTonotopyProj(proj,peakSignif,allToneIdx,behavToneIdx,cmap,normalized,smoothWindow )
    if nargin == 6
        smoothWindow = 1;
    end
    nData = length(proj);
    nPlanes = length(proj{1});
    toneStat = [];
    figure;
    temp = [];
    for i = 1:nData     
        for j = 1:nPlanes
            nTones = length(proj{i}{j});
            for k = 1:nTones
                temp = [temp;proj{i}{j}{k}];
            end
        end
    end
    edgeMax = max(temp); edgeMax = (ceil(edgeMax/100)+1)*100;
    edgeMin = min(temp); edgeMin = (floor(edgeMin/100)-1)*100;
    edges = edgeMin:50:edgeMax;
    for i = 1:nData     
        for j = 1:nPlanes
            subplot(nData,nPlanes,j+(i-1)*nPlanes); hold on;
            temp = proj{i}{j};
            nTones = length(temp);
            toneStat{i}{j} = zeros(nTones,3);
            for k = 1:nTones
                toneTemp = temp{k};
                toneID = find(peakSignif{i}{j}(k)==allToneIdx);
                if any(toneID == behavToneIdx)
                    linewid = 3;
                else
                    linewid = 1;
                end
                [counts,~] = histcounts(toneTemp,edges);
                allCounts{i}{j}(:,k) = counts;
                if normalized
                    counts = counts / sum(counts);
                end
                edgeMed = edges(1:end-1)/2 + edges(2:end)/2;
                plot(edgeMed, smoothdata(counts,'gaussian',smoothWindow),'Color',cmap(toneID,:),'LineWidth',linewid);
                toneStat{i}{j}(k,1) = median(toneTemp);
                toneStat{i}{j}(k,2) = prctile(toneTemp,75) - prctile(toneTemp,25);
                toneStat{i}{j}(k,3) = length(toneTemp);
            end
            title(['plane' int2str(j)]);
            xlabel('position (pixels)'); ylabel('frequency')
        end
        
    end
end

function [f1, targWidth, foilWidth] = plotTonotopyProjChange(toneStat,peakSignif,toneFreq,targIdx,foilIdx)
%myFun - Description
%
% Syntax: f1 = plotTonotopyProjChange(input)
%
% Long description
    nData = length(toneStat);
    nPlanes = length(toneStat{1});
    nTones = 17;
    f1 = figure;
    set(f1, 'Units', 'Normalized', 'OuterPosition', [0.025, 0.2, 0.95, 0.6]);
    for j = 1:nPlanes
        subplot(1,nPlanes,j)
        preTemp = toneStat{1}{j}; postTemp = toneStat{2}{j};
        preWidth = zeros(nTones,1); postWidth = zeros(nTones,1);
        preWidth(peakSignif{1}{j}) = preTemp(:,2);
        postWidth(peakSignif{2}{j}) = postTemp(:,2);
        plot(1:nTones, preWidth,'LineWidth',2); hold on;
        plot(1:nTones,postWidth,'LineWidth',2);
        ylimm = ylim;
        plot([targIdx targIdx],ylimm,'--g')
        plot([foilIdx foilIdx],ylimm,'--r')
        legend('preLearning','postLearning');
        title(['Plane' int2str(j)]);
        ticksLoc = 1:4:17;
        xlim([1 17])
        xticks(ticksLoc)
        xticklabels(strsplit(num2str(toneFreq(ticksLoc))));
        xlabel('frequency'); ylabel('width (pixels)')

        targWidth{j,1} = preWidth(targIdx);targWidth{j,2} = postWidth(targIdx);
        foilWidth{j,1} = preWidth(foilIdx);foilWidth{j,2} = postWidth(foilIdx);
    end

end

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