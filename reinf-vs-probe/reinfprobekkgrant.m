%% load data
clear; 
load('stimDffCorrected.mat')
% the data used here is the concatenation of 3 sessions, dff activity is 
% corrected for pretone baseline on a trial by trial basis
% the data is also concatenated version of both cd017 and cd036

% smooth the activity
smoothWindow = 5;
timeAxis = ((0:60)-15)/15;

targDff{1} = smoothdata(targDff{1},2,'gaussian',smoothWindow);
targDff{2} = smoothdata(targDff{2},2,'gaussian',smoothWindow);
foilDff{1} = smoothdata(foilDff{1},2,'gaussian',smoothWindow);
foilDff{2} = smoothdata(foilDff{2},2,'gaussian',smoothWindow);

% simply plot mean activity for visualization
meanTarg{1} = nanmean(targDff{1},3);
meanTarg{2} = nanmean(targDff{2},3);

meanFoil{1} = nanmean(foilDff{1},3);
meanFoil{2} = nanmean(foilDff{2},3);

figure;
subplot(1,2,1)
plot(timeAxis,nanmean(meanTarg{1},1),'LineWidth',2);hold on;
plot(timeAxis,nanmean(meanTarg{2},1),'LineWidth',2);
plot([0 0],[-0.015 0.045],'Color',[0.8 0.8 0.8])
legend('reinforced','probe'); title('Target Response')
ylim([-0.015 0.045]); ylabel('DF/F'); xlabel('Time After Tone Onset (s)')

subplot(1,2,2)
plot(timeAxis,nanmean(meanFoil{1},1),'LineWidth',2);hold on;
plot(timeAxis,nanmean(meanFoil{2},1),'LineWidth',2);
plot([0 0],[-0.015 0.045],'Color',[0.8 0.8 0.8])
legend('reinforced','probe'); title('Foil Response')
ylim([-0.015 0.045]); ylabel('DF/F'); xlabel('Time After Tone Onset (s)')

figure; 
imagesc(nanmean(targDff{2},3))
caxis([0 0.1])

toneWindow = 16:30;

tempAct = cat(3,targDff{1}(:,toneWindow,:),foilDff{1}(:,toneWindow,:));
tempAct = squeeze(nanmean(tempAct,3));
[~,idx] = max(tempAct,[],2);
idx = idx + 15;

meanTarg{1} = squeeze(nanmean(targDff{1}(:,idx,:),2));
meanFoil{1} = squeeze(nanmean(foilDff{1}(:,idx,:),2));

meanTarg{2} = squeeze(nanmean(targDff{2}(:,idx,:),2));
meanFoil{2} = squeeze(nanmean(foilDff{2}(:,idx,:),2));

%% significance test for responsive neurons for target/foil tones in reinf/probe context
targRespFlagR = ttest(meanTarg{1}',squeeze(nanmean(targDff{1}(:,10:15,:),2))','Tail', 'right' );
foilRespFlagR = ttest(meanFoil{1}',squeeze(nanmean(foilDff{1}(:,10:15,:),2))','Tail', 'right' );

targRespFlagP = ttest(meanTarg{2}',squeeze(nanmean(targDff{2}(:,10:15,:),2))','Tail', 'right' );
foilRespFlagP = ttest(meanFoil{2}',squeeze(nanmean(foilDff{2}(:,10:15,:),2))','Tail', 'right' );

%% significant test for modulation of reinf/probe on target tone
[targModFlag,targP] = ttest2(meanTarg{1}',meanTarg{2}');
[targP,targSort] = sort(targP,'ascend');

% one can check whether this result is similar with computing significance
% by each 10 trials on different day
%[h1,p] = ttest2(meanTarg{1}(:,1:105)',meanTarg{2}(:,1:10)');
%[h2,p] = ttest2(meanTarg{1}(:,106:105+140)',meanTarg{2}(:,11:20)');
%[h3,p] = ttest2(meanTarg{1}(:,246:end)',meanTarg{2}(:,21:30)');

[foilModFlag,foilP] = ttest2(meanFoil{1}',meanFoil{2}');
[foilP,foilSort] = sort(foilP,'ascend');

[tfModFlag,tfP] = ttest2(meanTarg{1}'-meanFoil{1}',meanTarg{2}'-meanFoil{2}','Alpha',0.025);
[tfP,tfSort] = sort(tfP,'ascend');
%% plot neurons with modulation of reinf/probe on target tone
for i = 1:10
    plotStuff(timeAxis,targDff, foilDff, targSort,i)
end
%% plot neurons with modulation of reinf/probe on foil tone
for i = 1:10
    plotStuff(timeAxis,targDff, foilDff, foilSort,i)
end
%% plot neurons with modulation of reinf/probe on selectivity (T-F)
for i = 1:10
    plotStuff(timeAxis,targDff, foilDff, tfSort,i)
end
%% plot bar graphs of reinf/probe modulation on target and foil tone
tempR = nanmean(meanTarg{1},2);
tempP = nanmean(meanTarg{2},2);
targMod = (abs(tempR) - abs(tempP)) ./ (abs(tempR) + abs(tempP));
tempR = nanmean(meanFoil{1},2);
tempP = nanmean(meanFoil{2},2);
foilMod = (abs(tempR) - abs(tempP)) ./ (abs(tempR) + abs(tempP));
tfMod = targMod - foilMod;

targModSignif = targMod( (targRespFlagR | targRespFlagP) & targModFlag);
[targModSort,targModIdx] = sort(targModSignif);
figure; subplot(1,2,1) 
bTarg = bar(targModSort, 'EdgeColor', 'None', 'FaceColor', [0.8500 0.3250 0.0980]); 
bTarg.FaceColor = 'flat'; tempIdx = find(targModSort>0);
bTarg.CData(tempIdx,:) = repmat([0 0.4470 0.7410], length(tempIdx),1);
xlabel('Neuron'); ylabel('Modulation Index'); title('Target'); ylim([-1 1])

foilModSignif = foilMod( (foilRespFlagR | foilRespFlagP) & foilModFlag);
[foilModSort,foilModIdx] = sort(foilModSignif);
subplot(1,2,2)
bFoil = bar(foilModSort, 'EdgeColor', 'None', 'FaceColor', [0.8500 0.3250 0.0980]); 
bFoil.FaceColor = 'flat'; tempIdx = find(foilModSort>0);
bFoil.CData(tempIdx,:) = repmat([0 0.4470 0.7410], length(tempIdx),1);
xlabel('Neuron'); ylabel('Modulation Index'); title('Foil'); ylim([-1 1])

% also plot it for modulation on selectivity
%tfModSignif = tfMod( (targRespFlagR | targRespFlagP | foilRespFlagR | foilRespFlagP) & tfModFlag);
%[tfModSort,tfModIdx] = sort(tfModSignif);
%figure;
%bFoil = bar(tfModSort, 'EdgeColor', 'None', 'FaceColor', [0.8500 0.3250 0.0980]); 
%bFoil.FaceColor = 'flat'; tempIdx = find(tfModSort>0);
%bFoil.CData(tempIdx,:) = repmat([0 0.4470 0.7410], length(tempIdx),1);
%xlabel('Neuron'); ylabel('Modulation Index'); title('Foil'); ylim([-1 1])

%% 
toneRespFlag = sum((foilRespFlagR | foilRespFlagP | targRespFlagR | targRespFlagP));
noMod = sum(~targModFlag & ~foilModFlag & toneRespFlag);
targModOnly = sum((targModFlag & toneRespFlag) - (targModFlag & foilModFlag & toneRespFlag));
foilModOnly = sum((foilModFlag & toneRespFlag) - (targModFlag & foilModFlag & toneRespFlag));
targFoilMod =  sum((targModFlag & foilModFlag & toneRespFlag));
figure;
pie([targModOnly foilModOnly targFoilMod noMod]./toneRespFlag,[1 1 1 0])
legend({'Target Only', 'Foil Only', 'Target and Foil', 'No Modulation'},'Location','Best')
%% function for plotting single cells
function plotStuff(timeaxis, targDff,foilDff,idx, i)
figure; 
subplot(1,2,1)
cont1 = plotStufff(timeaxis, {squeeze(targDff{1}(idx(i),:,:)),...
    squeeze(targDff{2}(idx(i),:,:))} );
title('Target'); hold on;

subplot(1,2,2)
cont2 = plotStufff(timeaxis, {squeeze(foilDff{1}(idx(i),:,:)),...
    squeeze(foilDff{2}(idx(i),:,:))});
title('Foil'); hold on;

minDff = min([cont1(:) ; cont2(:)]);
maxDff = max([cont1(:) ; cont2(:)]);

subplot(1,2,1)
plot([0 0],[minDff*1.2 maxDff*1.2],'Color',[0.8 0.8 0.8])
legend('reinforced','probe'); xlabel('Time after Tone Onset (s)'); ylabel('DF/F')
ylim([minDff*1.2 maxDff*1.2])
subplot(1,2,2)
plot([0 0],[minDff*1.2 maxDff*1.2],'Color',[0.8 0.8 0.8])
legend('reinforced','probe'); xlabel('Time after Tone Onset (s)'); ylabel('DF/F')
ylim([minDff*1.2 maxDff*1.2])

suptitle(['cell ' int2str(i)])
end

function cont = plotStufff(timeaxis, act)
sem{1} = nanstd(act{1},0,2) ./ sqrt(size(act{1},2));
sem{2} = nanstd(act{2},0,2) ./ sqrt(size(act{2},2));

cont{1} = [nanmean(act{1},2)+sem{1};flip(nanmean(act{1},2)-sem{1})];
cont{2} = [nanmean(act{2},2)+sem{2};flip(nanmean(act{2},2)-sem{2})];

plot(timeaxis,nanmean(act{1},2),'LineWidth',2); hold on; plot(timeaxis,nanmean(act{2},2),'LineWidth',2)
fill([timeaxis flip(timeaxis)],cont{1},[0 0.4470 0.7410], 'faceAlpha',0.2, 'edgeColor', 'None')
fill([timeaxis flip(timeaxis)],cont{2},[0.8500 0.3250 0.0980], 'faceAlpha',0.2, 'edgeColor', 'None')
cont = [cont{1} cont{2}];
end


