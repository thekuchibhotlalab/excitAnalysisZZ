mouse = 'cd017';
if strcmp(mouse,'cd017')
    baselineFile = [6 11 16 21 26 31 36 41 46 51 56 61 66 71 76 82];
elseif strcmp(mouse,'cd036')
    baselineFile = [8 13 18 22 26 30 34 38 42 46 50 54 58 62 66 70 71 73];
end

sessionList = [length(baselineFile)-3 length(baselineFile)-1];
%sessionList = [1 2];
%% compare dff and spikes

dffBaseline = {};
spikeBaseline = {};
smoothSpikeBaseline = {};
shuffleSpikeBaseline= {};

allCorrPairDff = cell(1,length(sessionList));
allCorrPairSpike = cell(1,length(sessionList));
allCorrPairSmoothSpike = cell(1,length(sessionList));
allCorrPairSpikeShuffle = cell(1,length(sessionList));

allCorrMatDff = cell(1,length(sessionList));
allCorrMatSpike = cell(1,length(sessionList));
allCorrMatSmoothSpike = cell(1,length(sessionList));
allCorrMatSpikeShuffle = cell(1,length(sessionList));

ishereBaseline{1} = sum(ishere{1}(:,sessionList+1),2); %+1 beacuse 1 day of tuning curve
ishereBaseline{2} = sum(ishere{2}(:,sessionList+1),2); %+1 beacuse 1 day of tuning curve
sum(ishereBaseline{1} == length(sessionList))
sum(ishereBaseline{2} == length(sessionList))
ishereBaseline{1} = (ishereBaseline{1}==length(sessionList));
ishereBaseline{2} = (ishereBaseline{2}==length(sessionList));

for j = 1:length(baselineFile)
    tempDff = [];
    tempSpike = [];
    tempSpikeSmooth = [];
    tempSpikeShuffle = [];
    for i = 1:nPlanes
        % use +(2-i) as a dirty way to get rid of the extra frame of plane 1 in the beginning
        tempDff = [tempDff;signals{i}(ishereBaseline{i},nFrames_oneplane(baselineFile(j),i)+1+(2-i):nFrames_oneplane(baselineFile(j)+1,i),2)];
        %dffBaseline{i,j} = signals{i}(ishereBaseline{i},nFrames_oneplane(baselineFile(j),i)+1:nFrames_oneplane(baselineFile(j)+1,i),2);
        tempSpike=[tempSpike;s_deconvolve_by_trial{i}(ishereBaseline{i},nFrames_oneplane(baselineFile(j),i)+1+(2-i):nFrames_oneplane(baselineFile(j)+1,i))];
        % smooth using 200ms moving window
        tempSpikeSmooth = [tempSpikeSmooth;smoothdata(s_deconvolve_by_trial{i}(ishereBaseline{i},nFrames_oneplane(baselineFile(j),i)+1+(2-i):nFrames_oneplane(baselineFile(j)+1,i)),'movmean',3)];
        tempSpikeShuffle = [tempSpikeShuffle;Shuffle(s_deconvolve_by_trial{i}(ishereBaseline{i},nFrames_oneplane(baselineFile(j),i)+1+(2-i):nFrames_oneplane(baselineFile(j)+1,i)))];
    end
    dffBaseline{j} = tempDff;
    spikeBaseline{j} = tempSpike; 
    smoothSpikeBaseline{j} = tempSpikeSmooth;
    shuffleSpikeBaseline{j} = tempSpikeShuffle;
end


tic;
nNeurons = size(dffBaseline{j},1);
for j = 1:length(sessionList) % the first and last session
    a = triu(ones(nNeurons));
    trid = (a==0);
    %for k = 1:size(allPairs,1)
    corrMat = corr(dffBaseline{sessionList(j)}');
    allPair = corrMat(trid);
    allCorrPairDff{j} = allPair(:);
    allCorrMatDff{j} = corrMat;
    
    corrMat = corr(spikeBaseline{sessionList(j)}');
    allPair = corrMat(trid);
    allCorrPairSpike{j} = allPair(:);
    allCorrMatSpike{j} = corrMat;

    corrMat = corr(smoothSpikeBaseline{sessionList(j)}');
    allPair = corrMat(trid);
    allCorrPairSmoothSpike{j} = allPair(:);
    allCorrMatSmoothSpike{j} = corrMat;
    
    corrMat = corr(shuffleSpikeBaseline{sessionList(j)}');
    allPair = corrMat(trid);
    allCorrPairSpikeShuffle{j} = allPair(:);
    allCorrMatSpikeShuffle{j} = corrMat;
    
end
toc;
%%
edges = -0.2:0.005:0.8;
figure;subplot(1,2,1);
h = histogram(allCorrPairDff{1,1},edges);h.EdgeColor = 'none'; hold on;
h = histogram(allCorrPairDff{1,2},edges);h.EdgeColor = 'none';
subplot(1,2,2);
h=histogram(allCorrPairDff{2,1},edges);h.EdgeColor = 'none'; hold on; 
h=histogram(allCorrPairDff{2,2},edges);h.EdgeColor = 'none';

figure;subplot(1,2,1);
h=histogram(allCorrPairSpike{1,1},edges);h.EdgeColor = 'none';  hold on; 
h=histogram(allCorrPairSpike{1,2},edges);h.EdgeColor = 'none'; 
subplot(1,2,2);
h=histogram(allCorrPairSpike{2,1},edges);h.EdgeColor = 'none';  hold on; 
h=histogram(allCorrPairSpike{2,2},edges);h.EdgeColor = 'none'; 

figure;subplot(1,2,1);
h=histogram(allCorrPairSmoothSpike{1,1},edges);h.EdgeColor = 'none';  hold on; 
h=histogram(allCorrPairSmoothSpike{1,2},edges);h.EdgeColor = 'none'; 
subplot(1,2,2);
h=histogram(allCorrPairSmoothSpike{2,1},edges);h.EdgeColor = 'none';  hold on; 
h=histogram(allCorrPairSmoothSpike{2,2},edges);h.EdgeColor = 'none'; 

figure;subplot(1,2,1);
h=histogram(allCorrPairSpikeShuffle{1,1},edges);h.EdgeColor = 'none';  hold on; 
h=histogram(allCorrPairSpikeShuffle{1,2},edges);h.EdgeColor = 'none'; 
subplot(1,2,2);
h=histogram(allCorrPairSpikeShuffle{2,1},edges);h.EdgeColor = 'none';  hold on; 
h=histogram(allCorrPairSpikeShuffle{2,2},edges);h.EdgeColor = 'none'; 
%% scatter plot of changes
limit = [-0.2 1];
dotSize = 8;
figure; subplot(4,2,1);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
scatter(allPairCorr{1,1},allPairCorr{1,2},dotSize,'.');xlim(limit);ylim(limit); 
subplot(4,2,2);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
scatter(allPairCorr{2,1},allPairCorr{2,2},dotSize,'.');xlim(limit);ylim(limit); 

subplot(4,2,3);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
scatter(allCorrPairSpike{1,1},allCorrPairSpike{1,2},dotSize,'.');xlim(limit);ylim(limit); 
subplot(4,2,4);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
scatter(allCorrPairSpike{2,1},allCorrPairSpike{2,2},dotSize,'.');xlim(limit);ylim(limit); 

subplot(4,2,5);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
scatter(allCorrPairSmoothSpike{1,1},allCorrPairSmoothSpike{1,2},dotSize,'.');xlim(limit);ylim(limit); 
subplot(4,2,6);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
scatter(allCorrPairSmoothSpike{2,1},allCorrPairSmoothSpike{2,2},dotSize,'.');xlim(limit);ylim(limit); 

subplot(4,2,7);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
scatter(allCorrPairSpikeShuffle{1,1},allCorrPairSpikeShuffle{1,2},dotSize,'.');xlim(limit);ylim(limit); 
subplot(4,2,8);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
scatter(allCorrPairSpikeShuffle{2,1},allCorrPairSpikeShuffle{2,2},dotSize,'.');xlim(limit);ylim(limit); 


%%  
spikeBaseline = {};
shuffleSpikeBaseline= {};

allCorrPairSpike = cell(1,length(sessionList));
allCorrPairSpikeShuffle = cell(1,length(sessionList));