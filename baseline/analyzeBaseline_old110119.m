mouse = 'cd036';
if strcmp(mouse,'cd017')
    baselineFile = [6 11 16 21 26 31 36 41 46 51 56 61 66 71 76 82];
elseif strcmp(mouse,'cd036')
    baselineFile = [8 13 18 22 26 30 34 38 42 46 50 54 58 62 66 70 71 73];
end
sessionList = 1:length(baselineFile);
nPlanes = 2;
%sessionList = [length(baselineFile)-3 length(baselineFile)-1];
%sessionList = [1 2];

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
%%
for j = 1:length(baselineFile)
    tempDff = [];
    tempSpike = [];
    tempSpikeSmooth = [];
    tempSpikeShuffle = [];
    for i = 1:nPlanes
        % use +(2-i) as a dirty way to get rid of the extra frame of plane 1 in the beginning
        tempDff = [tempDff;signals{i}(ishereBaseline{i},nFrames_oneplane(baselineFile(j),i)+1+(2-i):nFrames_oneplane(baselineFile(j)+1,i))];
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

%%
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
limit = [-0.2 1];
dotSize = 8;
edges = -0.2:0.005:0.8;
climCorr = [-0.3 0.3];
climSpike = [0 0.2];
corrMed = [];
corrStd = [];
varExp = {};
angleEV1 = [];
for i = 1:length(sessionList)-1
    figure(1);subplot_tight(5,4,i);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on;
    scatter(allCorrPairSpike{i},allCorrPairSpike{i+1},dotSize,'.');xlim(limit);ylim(limit); 
    title(['day' int2str(sessionList(i))])
    
    figure(2);subplot_tight(5,4,i);h = histogram(allCorrPairSpike{i},edges);h.EdgeColor = 'none'; hold on;
    h = histogram(allCorrPairSpike{i+1},edges);h.EdgeColor = 'none'; hold on;
    ylim([0 8000])
    title(['day' int2str(sessionList(i))])
    if (i==1)
        legend('day before', 'day after')
    end
    
    corrMed(i) = median(allCorrPairSpike{i});
    corrStd(i) = std(allCorrPairSpike{i});
    
end
figure(3)
plot(corrMed); hold on; plot(corrStd);
ylim([0 0.1])
legend('median pairwise corr','std of pairwise corr')

tempCorrMat = allCorrMatSpike{1};
[~,idx] = max(sum(tempCorrMat,2));
[~,corridx] = sortrows(tempCorrMat,idx,'descend');

[EV, ED] = eig(allCorrMatSpike{1});
ED = diag(ED);
[ED,idx] = sort(ED,'descend');
EV = EV(:,idx);
EV5 = {zeros(size(EV,1),length(sessionList)),zeros(size(EV,1),length(sessionList)),...
    zeros(size(EV,1),length(sessionList)),zeros(size(EV,1),length(sessionList)),...
    zeros(size(EV,1),length(sessionList))};
anglePC = {};
for i = 1:length(sessionList)
    figure(4);
    tempCorrMat = allCorrMatSpike{i};
    tempCorrMat = tempCorrMat(corridx,:);
    tempCorrMat = tempCorrMat(:,corridx);
    subplot_tight(5,4,i,[0.01 0.01]);imagesc(tempCorrMat); caxis(climCorr)
    set(gca,'xticklabels',[]);set(gca,'yticklabels',[])
    
    [EV, ED] = eig(allCorrMatSpike{i});
    ED = diag(ED);
    [ED,idx] = sort(ED,'descend');
    EV = EV(:,idx);
    varExp{i} = [cumsum(ED)/sum(ED)];
    for j = 1:5
        EV5{j}(:,i) = EV(:,j);
    end
    
    figure(5);
    subplot_tight(5,4,i,[0.01 0.01]);imagesc(spikeBaseline{i}); caxis(climSpike); 
    set(gca,'xticklabels',[]);set(gca,'yticklabels',[])
    
    figure(6);
    subplot_tight(5,4,i,[0.05 0.05]);plot(varExp{i});
    xlim([0 15])
    ylim([0 varExp{i}(15)+0.05])
    %set(gca,'xticklabels',[]);set(gca,'yticklabels',[])
    
end
for i = 1:5
    tempAngle = EV5{i}' * EV5{i};
    tempAngle(tempAngle<0) = -tempAngle(tempAngle<0);
    anglePC{i} = tempAngle;
end

figure(7);
for i = 1:5
    subplot(1,5,i)
    imagesc(anglePC{i},[-1 1])
end
%%
edges = -0.2:0.005:0.8;
figure;subplot(1,4,1);
h = histogram(allCorrPairDff{1,1},edges);h.EdgeColor = 'none'; hold on;
h = histogram(allCorrPairDff{1,2},edges);h.EdgeColor = 'none';
%subplot(1,2,2);
%h=histogram(allCorrPairDff{2,1},edges);h.EdgeColor = 'none'; hold on; 
%h=histogram(allCorrPairDff{2,2},edges);h.EdgeColor = 'none';

subplot(1,4,2);
h=histogram(allCorrPairSpike{1,1},edges);h.EdgeColor = 'none';  hold on; 
h=histogram(allCorrPairSpike{1,2},edges);h.EdgeColor = 'none'; 
%subplot(1,2,2);
%h=histogram(allCorrPairSpike{2,1},edges);h.EdgeColor = 'none';  hold on; 
%h=histogram(allCorrPairSpike{2,2},edges);h.EdgeColor = 'none'; 

subplot(1,4,3);
h=histogram(allCorrPairSmoothSpike{1,1},edges);h.EdgeColor = 'none';  hold on; 
h=histogram(allCorrPairSmoothSpike{1,2},edges);h.EdgeColor = 'none'; 
%subplot(1,2,2);
%h=histogram(allCorrPairSmoothSpike{2,1},edges);h.EdgeColor = 'none';  hold on; 
%h=histogram(allCorrPairSmoothSpike{2,2},edges);h.EdgeColor = 'none'; 

subplot(1,4,4);
h=histogram(allCorrPairSpikeShuffle{1,1},edges);h.EdgeColor = 'none';  hold on; 
h=histogram(allCorrPairSpikeShuffle{1,2},edges);h.EdgeColor = 'none'; 
%subplot(1,2,2);
%h=histogram(allCorrPairSpikeShuffle{2,1},edges);h.EdgeColor = 'none';  hold on; 
%h=histogram(allCorrPairSpikeShuffle{2,2},edges);h.EdgeColor = 'none'; 
%% scatter plot of changes
limit = [-0.2 1];
dotSize = 8;
figure; subplot(1,4,1);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
scatter(allCorrPairDff{1,1},allCorrPairDff{1,2},dotSize,'.');xlim(limit);ylim(limit); 
%subplot(4,2,2);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
%scatter(allPairCorr{2,1},allPairCorr{2,2},dotSize,'.');xlim(limit);ylim(limit); 

subplot(1,4,2);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
scatter(allCorrPairSpike{1,1},allCorrPairSpike{1,2},dotSize,'.');xlim(limit);ylim(limit); 
%subplot(4,2,4);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
%scatter(allCorrPairSpike{2,1},allCorrPairSpike{2,2},dotSize,'.');xlim(limit);ylim(limit); 

subplot(1,4,3);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
scatter(allCorrPairSmoothSpike{1,1},allCorrPairSmoothSpike{1,2},dotSize,'.');xlim(limit);ylim(limit); 
%subplot(4,2,6);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
%scatter(allCorrPairSmoothSpike{2,1},allCorrPairSmoothSpike{2,2},dotSize,'.');xlim(limit);ylim(limit); 

subplot(1,4,4);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
scatter(allCorrPairSpikeShuffle{1,1},allCorrPairSpikeShuffle{1,2},dotSize,'.');xlim(limit);ylim(limit); 
%subplot(4,2,8);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
%scatter(allCorrPairSpikeShuffle{2,1},allCorrPairSpikeShuffle{2,2},dotSize,'.');xlim(limit);ylim(limit); 