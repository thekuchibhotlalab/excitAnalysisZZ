%% load data
function getBaselineCorr(mouse)
%mouse = 'cd041';
datapath = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' mouse];
load([datapath '\ar1_allday_fixBadSession\' mouse '_calman_ar1_foo90_pars_allday_s_plane0.mat'],'s');
%load(['C:\Users\zzhu34\Documents\tempdata\deconv_test_' mouse...
%    '\allSessions\' mouse '_calman_ar2_foo95_optb_nosmin_c.mat'],'c');
try
    load(['C:\Users\zzhu34\Documents\tempdata\excitData\' mouse '\roi\ishere_plane0_final.mat'],'ishere');
catch
    load(['C:\Users\zzhu34\Documents\tempdata\excitData\' mouse '\roi\ishere_plane0.mat'],'ishere');  
end
if iscell(ishere)
    ishere = ishere{1};
end
%cellselect = sum(cell2mat(spikeBaseline),2)>1000;
% data contains variables: dffBaseline spikeBaseline smoothSpikeBaseline shuffleSpikeBaseline

%%
tic;
%figure; scatter(mean(s{6},2),mean(c{6},2));
configPath = 'C:\Users\zzhu34\Documents\tempdata\excitData\config\mouse';
filename = [configPath '\' mouse '_config.csv'];
configTable = readtable(filename);
baselineStemp = s(strcmp(configTable.BehavType,'Baseline'));
baselineS = cell2mat(baselineStemp);
baselineS = reshape(baselineS,size(baselineStemp{1},1),size(baselineStemp{1},2),length(baselineStemp));
baselineDay = 2:(size(ishere,2)-2);
selectNeuron = sum(ishere(:,baselineDay),2)==length(baselineDay);
baselineS = baselineS(selectNeuron,:,:); 
nanflag = (sum(sum(isnan(baselineS),2),3)>0); 
baselineS = baselineS(~nanflag,:,:);

nNeuron = size(baselineS,1);
%%
tridIdx = zeros(nNeuron, nNeuron); count = 0;
for i = 1:nNeuron
    tempIdx = (i+1):nNeuron;
    tridIdx(tempIdx,i) = count + (1:length(tempIdx));
    count = count + length(tempIdx);
end

%%

nRep = 2000;
trid = triu(ones(nNeuron));
trid = logical(1-trid);
nPair = sum(trid(:)>0);
allPair = zeros(nPair,size(baselineS,3)); 
allPairZscore =zeros(nPair,size(baselineS,3));
allPairH =zeros(nPair,size(baselineS,3));
for j = 1:size(baselineS,3) % the first and last session
    

    corrMat = corr(baselineS(:,:,j)');
    tempAllPair = corrMat(trid);
    allCorrMatSpike{j} = corrMat;
    [basis{j}, varExp{j}, proj{j}, covMat{j}] = fn_pca(baselineS(:,:,j),'zscore',false);

    allPairShuffle = zeros(length(tempAllPair),nRep);
    for i = 1:nRep
        shuffledS = shuffleS(baselineS(:,:,j));
        corrMat = corr(shuffledS');
        allPairShuffle(:,i) = corrMat(trid);
    end
    zscorePair = (tempAllPair - mean(allPairShuffle,2)) ./std(allPairShuffle,0,2);
    allPair(:,j) = tempAllPair;
    allPairZscore(:,j) = zscorePair;
    allPairH(:,j) = zscorePair > 3;
end
toc;
%% Plot results 

% plot 1 - distribution of significant (any day) correlation coeffecients
%colorScale = repmat(linspace(0,1,size(baselineS,3)),3,1) .* repmat([0;0;0],1,size(baselineS,3)) ...
%    + repmat(linspace(1,0,size(baselineS,3)),3,1) .* repmat([1;0.5;0.5],1,size(baselineS,3)) ;
%respFlag = sum(allPairH,2) > 0;
%figure; hold on; 
%for i = 1:size(baselineS,3)
    %fn_plotHistLine(allPair(respFlag,i));
%    h = cdfplot(allPair(respFlag,i)); set( h, 'Color', colorScale(:,i), 'LineWidth', 1.5);
%end
% plot 1.1 - distribution of significant (any day) correlation coeffecients
%figure; [~,sortIdx] = sort(sum(allPairH,2),'descend'); imagesc(allPairH(sortIdx,:));

save([datapath '\' mouse '_baseline_data_pairwiseCorr.mat'],'selectNeuron','baselineS','tridIdx',...
    'covMat','allPair','allPairZscore','allPairH');
t = toc;
disp([mouse ' is finished. Time = ' num2str(t,'%.2f') ' sec'])
end

%{
limit = [-0.2 1];
dotSize = 8;
edges = -0.2:0.005:0.8;
climCorr = [-0.3 0.3];
climSpike = [0 0.2];
corrMed = [];
corrStd = [];
varExp = {};
angleEV1 = [];
for i = 1:length(spikeBaseline)-1
    figure(1);subplot_tight(5,4,i);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on;
    scatter(allCorrPairSpike{i},allCorrPairSpike{i+1},dotSize,'.');xlim(limit);ylim(limit); 
    title(['day' int2str(i)])
    
    figure(2);subplot_tight(5,4,i);h = histogram(allCorrPairSpike{i},edges);h.EdgeColor = 'none'; hold on;
    h = histogram(allCorrPairSpike{i+1},edges);h.EdgeColor = 'none'; hold on;
    ylim([0 8000])
    title(['day' int2str(i)])
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
%}
function s = shuffleS(s)
for i = 1:size(s,1)
    tempIdx = randperm(size(s,2));
    s(i,:) = s(i,tempIdx);
end

end


%%

% [EV, ED] = eig(allCorrMatSpike{1});
% ED = diag(ED);
% [ED,idx] = sort(ED,'descend');
% EV = EV(:,idx);
% EV5 = {zeros(size(EV,1),length(sessionList)),zeros(size(EV,1),length(sessionList)),...
%     zeros(size(EV,1),length(sessionList)),zeros(size(EV,1),length(sessionList)),...
%     zeros(size(EV,1),length(sessionList))};
% anglePC = {};
% for i = 1:length(sessionList)
%     figure(4);
%     tempCorrMat = allCorrMatSpike{i};
%     tempCorrMat = tempCorrMat(corridx,:);
%     tempCorrMat = tempCorrMat(:,corridx);
%     subplot_tight(5,4,i,[0.01 0.01]);imagesc(tempCorrMat); caxis(climCorr)
%     set(gca,'xticklabels',[]);set(gca,'yticklabels',[])
%     
%     [EV, ED] = eig(allCorrMatSpike{i});
%     ED = diag(ED);
%     [ED,idx] = sort(ED,'descend');
%     EV = EV(:,idx);
%     varExp{i} = [cumsum(ED)/sum(ED)];
%     for j = 1:5
%         EV5{j}(:,i) = EV(:,j);
%     end
%     
%     figure(5);
%     subplot_tight(5,4,i,[0.01 0.01]);imagesc(spikeBaseline{i}); caxis(climSpike); 
%     set(gca,'xticklabels',[]);set(gca,'yticklabels',[])
%     
%     figure(6);
%     subplot_tight(5,4,i,[0.05 0.05]);plot(varExp{i});
%     xlim([0 15])
%     ylim([0 varExp{i}(15)+0.05])
%     %set(gca,'xticklabels',[]);set(gca,'yticklabels',[])
%     
% end
% for i = 1:5
%     tempAngle = EV5{i}' * EV5{i};
%     tempAngle(tempAngle<0) = -tempAngle(tempAngle<0);
%     anglePC{i} = tempAngle;
% end
% 
% figure(7);
% for i = 1:5
%     subplot(1,5,i)
%     imagesc(anglePC{i},[-1 1])
% end
%%
% edges = -0.2:0.005:0.8;
% figure;subplot(1,4,1);
% h = histogram(allCorrPairDff{1,1},edges);h.EdgeColor = 'none'; hold on;
% h = histogram(allCorrPairDff{1,2},edges);h.EdgeColor = 'none';
% %subplot(1,2,2);
% %h=histogram(allCorrPairDff{2,1},edges);h.EdgeColor = 'none'; hold on; 
% %h=histogram(allCorrPairDff{2,2},edges);h.EdgeColor = 'none';
% 
% subplot(1,4,2);
% h=histogram(allCorrPairSpike{1,1},edges);h.EdgeColor = 'none';  hold on; 
% h=histogram(allCorrPairSpike{1,2},edges);h.EdgeColor = 'none'; 
% %subplot(1,2,2);
% %h=histogram(allCorrPairSpike{2,1},edges);h.EdgeColor = 'none';  hold on; 
% %h=histogram(allCorrPairSpike{2,2},edges);h.EdgeColor = 'none'; 
% 
% subplot(1,4,3);
% h=histogram(allCorrPairSmoothSpike{1,1},edges);h.EdgeColor = 'none';  hold on; 
% h=histogram(allCorrPairSmoothSpike{1,2},edges);h.EdgeColor = 'none'; 
% %subplot(1,2,2);
% %h=histogram(allCorrPairSmoothSpike{2,1},edges);h.EdgeColor = 'none';  hold on; 
% %h=histogram(allCorrPairSmoothSpike{2,2},edges);h.EdgeColor = 'none'; 
% 
% subplot(1,4,4);
% h=histogram(allCorrPairSpikeShuffle{1,1},edges);h.EdgeColor = 'none';  hold on; 
% h=histogram(allCorrPairSpikeShuffle{1,2},edges);h.EdgeColor = 'none'; 
%subplot(1,2,2);
%h=histogram(allCorrPairSpikeShuffle{2,1},edges);h.EdgeColor = 'none';  hold on; 
%h=histogram(allCorrPairSpikeShuffle{2,2},edges);h.EdgeColor = 'none'; 
%% scatter plot of changes
% limit = [-0.2 1];
% dotSize = 8;
% figure; subplot(1,4,1);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
% scatter(allCorrPairDff{1,1},allCorrPairDff{1,2},dotSize,'.');xlim(limit);ylim(limit); 
% %subplot(4,2,2);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
% %scatter(allPairCorr{2,1},allPairCorr{2,2},dotSize,'.');xlim(limit);ylim(limit); 
% 
% subplot(1,4,2);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
% scatter(allCorrPairSpike{1,1},allCorrPairSpike{1,2},dotSize,'.');xlim(limit);ylim(limit); 
% %subplot(4,2,4);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
% %scatter(allCorrPairSpike{2,1},allCorrPairSpike{2,2},dotSize,'.');xlim(limit);ylim(limit); 
% 
% subplot(1,4,3);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
% scatter(allCorrPairSmoothSpike{1,1},allCorrPairSmoothSpike{1,2},dotSize,'.');xlim(limit);ylim(limit); 
% %subplot(4,2,6);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
% %scatter(allCorrPairSmoothSpike{2,1},allCorrPairSmoothSpike{2,2},dotSize,'.');xlim(limit);ylim(limit); 
% 
% subplot(1,4,4);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
% scatter(allCorrPairSpikeShuffle{1,1},allCorrPairSpikeShuffle{1,2},dotSize,'.');xlim(limit);ylim(limit); 
%subplot(4,2,8);plot(limit,limit,'Color',[0.8 0.8 0.8]); hold on
%scatter(allCorrPairSpikeShuffle{2,1},allCorrPairSpikeShuffle{2,2},dotSize,'.');xlim(limit);ylim(limit);
%