function CleanToneEvokedActivityFullLearning()

mice = {'cd017','cd036'};
RECAP_TONEEVOKED = cell(length(mice),6);
thismouse = 2;
%%
mouse = mice{thismouse};
results = GetPreprocessedInfo(mouse);
nFrames_oneplane = results{1};
signals = results{2};
matrix = results{3};
ishere = results{4};
ks = results{5};
ks2 = results{6};
ctx = results{7};
acq = results{8};
TONEF = results{9};
REWF = results{10};
acq = results{11};
dayss_behav = results{12};
startfrom = results{13};
% results = {nFrames_oneplane,signals,matrix,ishere,ks,ks2,ctx,acq,TONEF,REWF,acq,dayss_behav,startfrom};
nDays_behav = length(dayss_behav);
DAY = 1;
BLOC = 2;
TRIAL = 3;
CONTEXT = 4;
RESP = 5;
nPlanes = 2;

if strcmp(mouse,'cd017')
    nop = matrix(:,BLOC)==4; % remove bad trials
else
    nop = ~ismember(matrix(:,BLOC),matrix(:,BLOC)); % take everything
end

%% Find set of days with max number of same cells
global_here = [];
for p=1:nPlanes % loop through planes
    here = ishere{p}(:,dayss_behav+1);
    here(isnan(here) | here~=1) = 0;
    nCellsHere = sum(here); 

    global_here = [global_here logical(here)']; % days x neurons
end

nDayPerCell = sum(global_here);
nCellPerDays = sum(global_here,2);
figure;
subplot(2,1,1);hist(nCellPerDays,50);
subplot(2,1,2);hist(nDayPerCell,20);

% 1) remove the worst days
if strcmp(mouse,'cd017')
    nthreshold = 400;
elseif strcmp(mouse,'cd036')
    nthreshold = 800;
end
findgooddays = find(nCellPerDays>nthreshold);
nCellInGoodDays = sum(global_here(findgooddays,:),2);

% Find the best nb of days to have the best nb of cells
n = length(findgooddays); % start with all good days
results = [];
while n>1    
    takenCells = sum(global_here(findgooddays,:))>=n;
    nKeptCells = sum(takenCells);
    results = [results; [n nKeptCells]];
    n = n-1;
end

figure;plot(results(:,1),results(:,2),'.');

% decided to keep consistent cells from all the good days
n = length(findgooddays);
takenCells = sum(global_here(findgooddays,:))>=n;
nKeptCells = sum(takenCells);

% Just to check day identity
sum(global_here(findgooddays,takenCells),2);

% Stats
totcellsalldays = sum(sum(global_here)~=0);
totdays = size(global_here,1);
mcells = mean(sum(global_here(findgooddays,takenCells),2)./sum(global_here(findgooddays,:),2)*100);
semcells = sem(sum(global_here(findgooddays,takenCells),2)./sum(global_here(findgooddays,:),2)*100);
disp([num2str(mcells) '+-' num2str(semcells) ' % of cells kept over ' num2str(n/totdays*100) '% of days']);

%% MEDIAN TRACE PER TONE (ALL PLANES) PER DAY
ndays = length(findgooddays);

pretone = 1; % for PSTH tone, in s
posttone = 4; % for PSTH tone, in s
nframes_psth = pretone*round(acq) + posttone*round(acq);

ms = nan(nframes_psth,ndays,2);
ss = nan(nframes_psth,ndays,2);

msprobe = nan(nframes_psth,ndays,2);
ssprobe = nan(nframes_psth,ndays,2);

sfull = cell(ndays,2);
ssfull = cell(ndays,2);

sfulltrial = cell(ndays,2);
sfulltrialprobe = cell(ndays,2);

sprobefull = cell(ndays,2);
ssprobefull = cell(ndays,2);

figure;hold on;

for j=1:ndays  
    day = findgooddays(j);
    this_day = matrix(:,DAY)==day+1;    
    for k=1:2        
        for p=1:nPlanes % loop through planes
            signals_thisplane = signals{p};
            trace = signals_thisplane(:,:,2);

            if p==2, prevncells = nCellsHere; end % continue matrix

            nCellsHere = size(trace,1);

            ok = this_day & ctx(:,1) & ~nop & ks2(:,k); % all tone (i.e. all trials)
            l = sum(ok); % nb of (these) tones
            matidx = repelem(0:nframes_psth-1,l,1);
            idx = (matidx+repelem(matrix(ok,TONEF(p))- pretone*round(acq),1,nframes_psth))';

            ok_probe = this_day & ctx(:,2) & ~nop & ks2(:,k); % all tone (i.e. all trials)
            l_probe = sum(ok_probe); % nb of (these) tones
            matidx_probe = repelem(0:nframes_psth-1,l_probe,1);
            idx_probe = (matidx_probe+repelem(matrix(ok_probe,TONEF(p))- pretone*round(acq),1,nframes_psth))';

            if p==1
                s = trace(:,idx(:))';    
                s = reshape(s,nframes_psth,l,nCellsHere);     

                s_probe = trace(:,idx_probe(:))';    
                s_probe = reshape(s_probe,nframes_psth,l_probe,nCellsHere);     
            else
                s2 = trace(:,idx(:))'; 
                s2 = reshape(s2,nframes_psth,l,nCellsHere); 
                s(:,:,prevncells+1:prevncells+nCellsHere) = s2;

                s2_probe = trace(:,idx_probe(:))'; 
                s2_probe = reshape(s2_probe,nframes_psth,l_probe,nCellsHere); 
                s_probe(:,:,prevncells+1:prevncells+nCellsHere) = s2_probe;            
            end
        end

        if length(takenCells)~=size(s,3), error('Error size'); end
        
        ms(:,j,k) = mean(squeeze(median(s(:,:,takenCells),2)),2);
        ss(:,j,k) = sem(squeeze(median(s(:,:,takenCells),2))')';
        
        sfull{j,k} = squeeze(median(s(:,:,takenCells),2));
        ssfull{j,k} = squeeze(std(s(:,:,takenCells),0,2));
        
        sfulltrial{j,k} = s(:,:,takenCells);
        sfulltrialprobe{j,k} = s_probe(:,:,takenCells);
        
                
        msprobe(:,j,k) = mean(squeeze(median(s_probe(:,:,takenCells),2)),2);
        ssprobe(:,j,k) = sem(squeeze(median(s_probe(:,:,takenCells),2))')';
        
        sprobefull{j,k} = squeeze(median(s_probe(:,:,takenCells),2));
        ssprobefull{j,k} = squeeze(std(s_probe(:,:,takenCells),0,2));  
    end
end
%%
n = length(findgooddays);
contextNames = {'Reinforced','Probe'};

figure;
colors_target = [0*ones(n,1) linspace(0.2,0.8,n)' 0*ones(n,1)];
colors_foil = [linspace(0.2,0.8,n)' 0*ones(n,1) 0*ones(n,1)];   
xmax = nframes_psth-25;
tonenames = {'Target','Foil'};
for j=1:2
    subplot(2,2,j+(j-1));hold on;
    for i=1:n
        if j==1
            shadedErrorBar(1:xmax,ms(1:xmax,i,j),ss(1:xmax,i,j),{'color',colors_target(i,:)});
        else
            shadedErrorBar(1:xmax,ms(1:xmax,i,j),ss(1:xmax,i,j),{'color',colors_foil(i,:)});
        end
    end
    ylim([-0.04 0.1]);
    PlotHVLines(pretone*round(acq),'v','k','linewidth',1);   
    PlotHVLines(0,'h','k:','linewidth',1);
    title(tonenames{j});    
    
    subplot(2,2,j+(j-1)+1);hold on;
    for i=1:n
        tms = ms(1:xmax,i,j)-mean(ms(1:pretone*round(acq),i,j));
        if j==1
            shadedErrorBar(1:xmax,tms,ss(1:xmax,i,j),{'color',colors_target(i,:)});
        else
            shadedErrorBar(1:xmax,tms,ss(1:xmax,i,j),{'color',colors_foil(i,:)});
        end
    end
    ylim([-0.04 0.1]);
    PlotHVLines(pretone*round(acq),'v','k','linewidth',1);   
    PlotHVLines(0,'h','k:','linewidth',1);
    title([contextNames{1} ' - ' tonenames{j}]);
end

figure;
tonenames = {'Target','Foil'};
for j=1:2
    subplot(2,2,j+(j-1));hold on;
    for i=1:n
        if j==1
            shadedErrorBar(1:xmax,msprobe(1:xmax,i,j),ssprobe(1:xmax,i,j),{'color',colors_target(i,:)});
        else
            shadedErrorBar(1:xmax,msprobe(1:xmax,i,j),ssprobe(1:xmax,i,j),{'color',colors_foil(i,:)});
        end
    end
    ylim([-0.04 0.1]);
    PlotHVLines(pretone*round(acq),'v','k','linewidth',1);   
    PlotHVLines(0,'h','k:','linewidth',1);
    title(tonenames{j});    
    
    subplot(2,2,j+(j-1)+1);hold on;
    for i=1:n
        tms = msprobe(1:xmax,i,j)-mean(msprobe(1:pretone*round(acq),i,j));
        if j==1
            shadedErrorBar(1:xmax,tms,ssprobe(1:xmax,i,j),{'color',colors_target(i,:)});
        else
            shadedErrorBar(1:xmax,tms,ssprobe(1:xmax,i,j),{'color',colors_foil(i,:)});
        end
    end
    ylim([-0.04 0.1]);
    PlotHVLines(pretone*round(acq),'v','k','linewidth',1);   
    PlotHVLines(0,'h','k:','linewidth',1);
    title([contextNames{2} ' - ' tonenames{j}]);
end
%% Same but as color curves
n = length(findgooddays);

xmax = nframes_psth-25;
xrange = [1 xmax];
tonenames = {'Target','Foil'};
cutoffs = [-0.02 0.05];

figure;
for j=1:2
    subplot(2,2,j+(j-1));hold on;
    PlotColorCurves(flip(squeeze(ms(1:xmax,1:n,j))'),xrange,'cutoffs',cutoffs);
%         PlotColorCurves(squeeze(ms(1:xmax,1:n,j))',xrange,'cutoffs',cutoffs);

    PlotHVLines(pretone*round(acq)+1,'v','w','linewidth',1);   
    title(tonenames{j});    
    
    subplot(2,2,j+(j-1)+1);hold on;
    tms = ms(1:xmax,1:n,j)-mean(ms(1:pretone*round(acq),1:n,j));
    PlotColorCurves(flip(squeeze(tms)'),xrange,'cutoffs',cutoffs);
%     PlotColorCurves(squeeze(tms)',xrange,'cutoffs',cutoffs);
    PlotHVLines(pretone*round(acq)+1,'v','w','linewidth',1);   
    title([contextNames{1} ' - ' tonenames{j}]);
end

figure;
for j=1:2
    subplot(2,2,j+(j-1));hold on;
    PlotColorCurves(flip(squeeze(msprobe(1:xmax,1:n,j))'),xrange,'cutoffs',cutoffs);
    PlotHVLines(pretone*round(acq)+1,'v','w','linewidth',1);   
    title(tonenames{j});    
    
    subplot(2,2,j+(j-1)+1);hold on;
    tms = msprobe(1:xmax,1:n,j)-mean(msprobe(1:pretone*round(acq),1:n,j));
    PlotColorCurves(flip(squeeze(tms)'),xrange,'cutoffs',cutoffs);
    PlotHVLines(pretone*round(acq)+1,'v','w','linewidth',1);   
    title([contextNames{2} ' - ' tonenames{j}]);
end
%% Plot individual median psth

% Amplitude and latency
pretone = 1; % for PSTH tone, in s
posttone = 4; % for PSTH tone, in s
xrange = [pretone pretone*round(acq)+posttone*round(acq)];

cutoffs = [-0.03 0.05];


% ORDERED BY TARGET
figure;
smoothfactor = [0 1];
ndays = length(findgooddays);
for i=1:ndays
%     i = findgooddays(ff);

    subplot(ndays,2,i+(i-1));hold on;
    targettraces = sfull{i,1};
    if i==1
        [~,ordertarget] = sort(mean(targettraces(15:25,:)));
    end
    PlotColorCurves(Smooth(targettraces(:,ordertarget)',smoothfactor),xrange,'cutoffs',cutoffs);
    PlotHVLines(pretone*round(acq)+1,'v','w');
    
    subplot(ndays,2,i+(i-1)+1);hold on;
    foiltraces = sfull{i,2};
    PlotColorCurves(Smooth(foiltraces(:,ordertarget)',smoothfactor),xrange,'cutoffs',cutoffs);   
    PlotHVLines(pretone*round(acq)+1,'v','w');
end

figure;
for i=1:ndays
%     i = findgooddays(ff);    
    
    subplot(ndays,2,i+(i-1));hold on;
    targettraces = sprobefull{i,1};  
%     if ff==1
%         [~,ordertarget] = sort(mean(targettraces(15:25,:)));
%     end
    PlotColorCurves(Smooth(targettraces(:,ordertarget)',smoothfactor),xrange,'cutoffs',cutoffs);
    PlotHVLines(pretone*round(acq)+1,'v','w');
    
    subplot(ndays,2,i+(i-1)+1);hold on;
    foiltraces = sprobefull{i,2};
    PlotColorCurves(Smooth(foiltraces(:,ordertarget)',smoothfactor),xrange,'cutoffs',cutoffs);    
    PlotHVLines(pretone*round(acq)+1,'v','w');
end

% ORDERED BY FOIL
figure;
smoothfactor = [0 1];
ndays = 5;
for ff=1:ndays
    i = findgooddays(ff);
   
    subplot(ndays,2,i+(i-1)+1);hold on;
    foiltraces = sfull{i,2};
    [~,orderfoil] = sort(mean(foiltraces(15:25,:)));
    PlotColorCurves(Smooth(foiltraces(:,orderfoil)',smoothfactor),xrange,'cutoffs',cutoffs);   
    PlotHVLines(pretone*round(acq)+1,'v','w');
    
    subplot(ndays,2,i+(i-1));hold on;
    targettraces = sfull{i,1};
    PlotColorCurves(Smooth(targettraces(:,orderfoil)',smoothfactor),xrange,'cutoffs',cutoffs);
    PlotHVLines(pretone*round(acq)+1,'v','w');
end


figure;
for ff=1:ndays
    i = findgooddays(ff);    
    
    subplot(ndays,2,i+(i-1)+1);hold on;
    foiltraces = sprobefull{i,2};
    [~,orderfoil] = sort(mean(foiltraces(15:25,:)));
    PlotColorCurves(Smooth(foiltraces(:,orderfoil)',smoothfactor),xrange,'cutoffs',cutoffs);    
    PlotHVLines(pretone*round(acq)+1,'v','w');
    
    subplot(ndays,2,i+(i-1));hold on;
    targettraces = sprobefull{i,1};    
    PlotColorCurves(Smooth(targettraces(:,orderfoil)',smoothfactor),xrange,'cutoffs',cutoffs);
    PlotHVLines(pretone*round(acq)+1,'v','w');
end
%% Plot for fig SfN
if strcmp(mouse,'cd017')
    daytoplot = findgooddays;
elseif strcmp(mouse,'cd036')
    daytoplot = findgooddays;
end

% Amplitude and latency
pretone = 1; % for PSTH tone, in s
posttone = 4; % for PSTH tone, in s
xrange = [pretone pretone*round(acq)+posttone*round(acq)];

cutoffs = [-0.03 0.05];

% ORDERED BY TARGET
figure;
smoothfactor = [0 1];
ndays = length(daytoplot);
for i=1:ndays

    subplot(ndays,2,i+(i-1));hold on;
    targettraces = sfull{i,1};
    if i==1
        [~,ordertarget] = sort(mean(targettraces(15:25,:)));
    end
    PlotColorCurves(Smooth(targettraces(:,ordertarget)',smoothfactor),xrange,'cutoffs',cutoffs);
    PlotHVLines(pretone*round(acq)+1,'v','w');
    
    subplot(ndays,2,i+(i-1)+1);hold on;
    foiltraces = sfull{i,2};
    PlotColorCurves(Smooth(foiltraces(:,ordertarget)',smoothfactor),xrange,'cutoffs',cutoffs);   
    PlotHVLines(pretone*round(acq)+1,'v','w');
    title([contextNames{1} ' - Day ' num2str(daytoplot(i))]);
end

figure;
for i=1:ndays
    
    subplot(ndays,2,i+(i-1));hold on;
    targettraces = sprobefull{i,1};  
%     if ff==1
%         [~,ordertarget] = sort(mean(targettraces(15:25,:)));
%     end
    PlotColorCurves(Smooth(targettraces(:,ordertarget)',smoothfactor),xrange,'cutoffs',cutoffs);
    PlotHVLines(pretone*round(acq)+1,'v','w');
    
    subplot(ndays,2,i+(i-1)+1);hold on;
    foiltraces = sprobefull{i,2};
    PlotColorCurves(Smooth(foiltraces(:,ordertarget)',smoothfactor),xrange,'cutoffs',cutoffs);    
    PlotHVLines(pretone*round(acq)+1,'v','w');
    title([contextNames{2} ' - Day ' num2str(daytoplot(i))]);
end


figure;
colors_target = [0*ones(ndays,1) linspace(0.2,0.8,ndays)' 0*ones(ndays,1)];
colors_foil = [linspace(0.2,0.8,ndays)' 0*ones(ndays,1) 0*ones(ndays,1)];   
xmax = nframes_psth;
tonenames = {'Target','Foil'};
for ff=1:ndays
    for j=1:2
        subplot(ndays,2,j+(ff-1)*2);hold on;
        i = daytoplot(ff);
        if j==1
            tms = ms(1:xmax,i,j)-mean(ms(1:pretone*round(acq),i,j));
            shadedErrorBar(1:xmax,tms,ss(1:xmax,i,j),{'color',colors_target(ff,:)});
            tmsprobe = msprobe(1:xmax,i,j)-mean(msprobe(1:pretone*round(acq),i,j));
            shadedErrorBar(1:xmax,tmsprobe,ssprobe(1:xmax,i,j),{'color',colors_target(ff,:)});
        elseif j==2
            tms = ms(1:xmax,i,j)-mean(ms(1:pretone*round(acq),i,j));
            shadedErrorBar(1:xmax,tms,ss(1:xmax,i,j),{'color',colors_foil(ff,:)});
            tmsprobe = msprobe(1:xmax,i,j)-mean(msprobe(1:pretone*round(acq),i,j));
            shadedErrorBar(1:xmax,tmsprobe,ssprobe(1:xmax,i,j),{'color',colors_foil(ff,:)});
        end
        ylim([-0.02 0.06]);
        PlotHVLines(pretone*round(acq),'v','k','linewidth',1);   
        PlotHVLines(0,'h','k:','linewidth',1);
        title(['Day ' num2str(i) ' - ' tonenames{j}]);   
    end
end
    
%% Get selectivity for each neurons
peakwind = 17:26;
selectivityidx = nan(nKeptCells,ndays);
selectivityidxprobe = nan(nKeptCells,ndays);
for i=1:ndays
    targettraces = sfull{i,1};
    foiltraces = sfull{i,2};
    
    selectivityidx(:,i) = abs((max(targettraces(peakwind,:))-max(foiltraces(peakwind,:)))...
        ./sum([max(abs(targettraces(peakwind,:)));max(abs(foiltraces(peakwind,:)))]));
    
    targettracesprobe = sprobefull{i,1};
    foiltracesprobe = sprobefull{i,2};    
    
    selectivityidxprobe(:,i) = abs((max(targettracesprobe(peakwind,:))-max(foiltracesprobe(peakwind,:)))...
        ./sum([max(abs(targettracesprobe(peakwind,:)));max(abs(foiltracesprobe(peakwind,:)))]));
end

smoothfactor1 = [0 1];
smoothfactor = 0;
%% Evolution of idx selectivity (deltaSI)
x = findgooddays;
nclusters=4;
smoothselectivity = Smooth(selectivityidx,smoothfactor1);
deltaselectivityidx = smoothselectivity - repmat(smoothselectivity(:,1),1,ndays);
idxreinf = kmeans(Smooth(deltaselectivityidx,smoothfactor),nclusters);

if strcmp(mouse,'cd017')
    smoothselectivityProbe = Smooth(selectivityidxprobe(:,2:end),smoothfactor1);
    deltaselectivityidxprobe = smoothselectivityProbe - repmat(smoothselectivityProbe(:,2:end),1,ndays-1);
else
    smoothselectivityProbe = Smooth(selectivityidxprobe,smoothfactor1);
    deltaselectivityidxprobe = smoothselectivityProbe - repmat(smoothselectivityProbe(:,1),1,ndays);
end
idxprobe = kmeans(Smooth(deltaselectivityidxprobe,smoothfactor),nclusters);


figure;
subplot(2,3,1);hold on;
plot(x,Smooth(deltaselectivityidx(idxreinf==1,:),smoothfactor)','r-');
hold on;plot(x,Smooth(deltaselectivityidx(idxreinf==2,:),smoothfactor)','b-');
hold on;plot(x,Smooth(deltaselectivityidx(idxreinf==3,:),smoothfactor)','g-');
hold on;plot(x,Smooth(deltaselectivityidx(idxreinf==4,:),smoothfactor)','k-');
ylim([-1 1]);

subplot(2,3,2);hold on;
plot(x,Smooth(deltaselectivityidx(idxprobe==1,:),smoothfactor)','r-');
hold on;plot(x,Smooth(deltaselectivityidx(idxprobe==2,:),smoothfactor)','b-');
hold on;plot(x,Smooth(deltaselectivityidx(idxprobe==3,:),smoothfactor)','g-');
hold on;plot(x,Smooth(deltaselectivityidx(idxprobe==4,:),smoothfactor)','k-');
ylim([-1 1]);

mgp = nan(nclusters,ndays);
semgrp = nan(nclusters,ndays);
mgpprobe = nan(nclusters,ndays);
semgrpprobe = nan(nclusters,ndays);
colors = {'r','b','g','k'};
for i=1:nclusters
    subplot(2,3,3);hold on;
    mgp(i,:) = mean(Smooth(deltaselectivityidx(idxreinf==i,:),smoothfactor));
    semgrp(i,:) = sem(Smooth(deltaselectivityidx(idxreinf==i,:),smoothfactor));
    shadedErrorBar(x,mgp(i,:),semgrp(i,:),{'color',colors{i}});
    ylim([-1 1]);
    subplot(2,3,4);hold on;
    mgpprobe(i,:) = mean(Smooth(deltaselectivityidx(idxprobe==i,:),smoothfactor));
    semgrpprobe(i,:) = sem(Smooth(deltaselectivityidx(idxprobe==i,:),smoothfactor));
    shadedErrorBar(x,mgpprobe(i,:),semgrpprobe(i,:),{'color',colors{i}});
    ylim([-1 1]);
end
PlotHVLines(0,'h','k:');

subplot(2,3,5);hold on;
plot(x,Smooth(deltaselectivityidx,smoothfactor)','color',[0.8 0.8 0.8]);
subplot(2,3,6);hold on;
plot(x,Smooth(deltaselectivityidx,smoothfactor)','color',[0.8 0.8 0.8]);
for i=1:nclusters
    subplot(2,3,5);hold on;
    shadedErrorBar(x,mgp(i,:),semgrp(i,:),{'color',colors{i}});
    ylim([-1 1]);
    subplot(2,3,6);hold on;
    shadedErrorBar(x,mgpprobe(i,:),semgrpprobe(i,:),{'color',colors{i}});
    ylim([-1 1]);
end



figure;
if strcmp(mouse,'cd017')
    xprobe=x(2:end);
else
    xprobe=x;
end
subplot(2,3,1);hold on;
plot(xprobe,Smooth(deltaselectivityidxprobe(idxprobe==1,:),smoothfactor)','r-');
hold on;plot(xprobe,Smooth(deltaselectivityidxprobe(idxprobe==2,:),smoothfactor)','b-');
hold on;plot(xprobe,Smooth(deltaselectivityidxprobe(idxprobe==3,:),smoothfactor)','g-');
hold on;plot(xprobe,Smooth(deltaselectivityidxprobe(idxprobe==4,:),smoothfactor)','k-');
ylim([-1 1]);

subplot(2,3,2);hold on;
plot(xprobe,Smooth(deltaselectivityidxprobe(idxreinf==1,:),smoothfactor)','r-');
hold on;plot(xprobe,Smooth(deltaselectivityidxprobe(idxreinf==2,:),smoothfactor)','b-');
hold on;plot(xprobe,Smooth(deltaselectivityidxprobe(idxreinf==3,:),smoothfactor)','g-');
hold on;plot(xprobe,Smooth(deltaselectivityidxprobe(idxreinf==4,:),smoothfactor)','k-');
ylim([-1 1]);

mgp = nan(nclusters,length(xprobe));
semgrp = nan(nclusters,length(xprobe));
mgpprobe = nan(nclusters,ndays);
semgrpprobe = nan(nclusters,ndays);
for i=1:nclusters
    subplot(2,3,3);hold on;
    mgp(i,:) = mean(Smooth(deltaselectivityidxprobe(idxprobe==i,:),smoothfactor));
    semgrp(i,:) = sem(Smooth(deltaselectivityidxprobe(idxprobe==i,:),smoothfactor));
    shadedErrorBar(xprobe,mgp(i,:),semgrp(i,:),{'color',colors{i}});
    ylim([-1 1]);
    subplot(2,3,4);hold on;
    mgpprobe(i,:) = mean(Smooth(deltaselectivityidxprobe(idxreinf==i,:),smoothfactor));
    semgrpprobe(i,:) = sem(Smooth(deltaselectivityidxprobe(idxreinf==i,:),smoothfactor));
    shadedErrorBar(xprobe,mgpprobe(i,:),semgrpprobe(i,:),{'color',colors{i}});
    ylim([-1 1]);
end
PlotHVLines(0,'h','k:');

subplot(2,3,5);hold on;
plot(x,Smooth(deltaselectivityidxprobe,smoothfactor)','color',[0.8 0.8 0.8]);
subplot(2,3,6);hold on;
plot(x,Smooth(deltaselectivityidxprobe,smoothfactor)','color',[0.8 0.8 0.8]);
for i=1:nclusters
    subplot(2,3,5);hold on;
    shadedErrorBar(x,mgp(i,:),semgrp(i,:),{'color',colors{i}});
    subplot(2,3,6);hold on;
    shadedErrorBar(x,mgpprobe(i,:),semgrpprobe(i,:),{'color',colors{i}});
end

%% DEFINE TONE-RESPONSIVE CELLS IN BOTH CONTEXTS

pretone = 1; % for PSTH tone, in s
posttone = 4; % for PSTH tone, in s
nframes_psth = pretone*round(acq) + posttone*round(acq);

showfig = false;
alphas = [0.001,0.05];
toneresponsive = cell(nPlanes,1);
toneselective = cell(nPlanes,1);
for p=1:nPlanes
    signals_thisplane = signals{p};
    nCells = size(signals_thisplane,1);
    
    toneresponsive{p} = nan(nCells,nDays_behav,2);
    toneselective{p} = cell(2,1);
    
    for i=1:nCells
        for j=1:nDays_behav            
            this_day = matrix(:,DAY)==dayss_behav(j)+1; 
            trace = signals_thisplane(i,:,2);
            
            here = ishere{p}(i,dayss_behav(j)+1);
            here(isnan(here) | here~=1) = 0;
            here = logical(here);
            if ~here, continue, end
    
            for ct = 1:2
                ok = this_day & ctx(:,ct);
                l = sum(ok);            
                matidx = repelem(0:nframes_psth-1,l,1);
                idx = (matidx+repelem(matrix(ok,TONEF(p))- pretone*round(acq),1,nframes_psth))';
                s = trace(:,idx(:))';
                s = reshape(s,nframes_psth,l);
                ms = median(s,2);
                bef=median(s(5:14,:));
                aft=median(s(17:26,:));
                
                [h0,pttest] = ttest(bef(:),aft(:),'alpha',alphas(ct));  
                toneresponsive{p}(i,j,ct) = h0;
                
                yep = nan(2,1);
                for k=1:2
                    ok = this_day & ctx(:,ct) & ks2(:,k);
                    l = sum(ok);            
                    matidx = repelem(0:nframes_psth-1,l,1);
                    idx = (matidx+repelem(matrix(ok,TONEF(p))- pretone*round(acq),1,nframes_psth))';
                    s = trace(:,idx(:))';
                    s = reshape(s,nframes_psth,l);
                    ms = median(s,2);
                    bef=median(s(5:14,:));
                    aft=median(s(17:26,:));
                    [h0,pttest] = ttest(bef(:),aft(:),'alpha',alphas(ct));  
                    yep(k) = h0;
                end
                toneselective{p}{ct}(i,j,:) = logical(yep);

                if showfig
                    fig=figure;hold on;
                    plot(s,'color',[0.8 0.8 0.8]);
                    plot(ms,'k','linewidth',2);
                    ylim([-0.2 0.2]);                    
                    PlotHVLines(pretone*round(acq),'v','color','k','linewidth',1);  
                    title(['Cell ' num2str(i) ', n=' num2str(length(bef)) ' ' contextNames{ct} ' trials, p=' ...
                        num2str(pttest) '(' num2str(h0) ')']); 
                    PlotIntervals([[5 14];[17 26]]);
                    pause();
                    close(fig);
                end
            end
        end
    end
end

%%
% Stats
respp = [];selecc = [];
for p=1:nPlanes % loop through planes
    respp = [respp;toneresponsive{p}]; % neurons x days x contexts
    seleccthisplane = [];
    for ct=1:2
        seleccthisplane(:,:,:,ct) = toneselective{p}{ct};
    end
    selecc = [selecc;seleccthisplane]; % neurons x days x tones x ctx
end
thesetakenCells = find(takenCells);

colors_target = [0*ones(nGoodDays,1) linspace(0.2,0.8,nGoodDays)' 0*ones(nGoodDays,1)];
colors_foil = [linspace(0.2,0.8,nGoodDays)' 0*ones(nGoodDays,1) 0*ones(nGoodDays,1)]; 
figure;
xrange = [1 nframes_psth];
cutoffs = [-0.1 0.1];
for i=1:nGoodDays
    subplot(4,4,i);hold on;
    day = findgooddays(i);
    oktarget = selecc(takenCells,day,1,1);
    tracereinftarget = sfull{i,1}(:,logical(oktarget));
%     tracereinftarget = sfull{i,1};
%     tracereinftarget(:,logical(oktarget)) = nan;
    PlotColorCurves(tracereinftarget',xrange,'cutoffs',cutoffs);
end
colormap(newCmap)

figure;
for i=1:nGoodDays
    subplot(2,2,1);hold on;
    day = findgooddays(i);
    oktarget = selecc(takenCells,day,1,1);
    tracereinftarget = sfull{i,1}(:,logical(oktarget));
    mtracereinftarget = mean(tracereinftarget,2);
    plot(1:nframes_psth,mtracereinftarget,'color',colors_target(i,:));
    
    subplot(2,2,2);hold on;
    day = findgooddays(i);
    okfoil = selecc(takenCells,day,2,1);
    tracereinffoil = sfull{i,2}(:,logical(okfoil));
    mtracereinffoil = mean(tracereinffoil,2);
    plot(1:nframes_psth,mtracereinffoil,'color',colors_foil(i,:));
    
    subplot(2,2,3);hold on;
    day = findgooddays(i);
    oktarget = selecc(takenCells,day,1,2);
    traceprobefoil = sprobefull{i,2}(:,logical(oktarget));
    mtraceprobefoil = mean(traceprobefoil,2);
    plot(1:nframes_psth,mtraceprobefoil,'color',colors_target(i,:));
    
    subplot(2,2,4);hold on;
    day = findgooddays(i);
    okfoil = selecc(takenCells,day,1,2);
    traceprobefoil = sprobefull{i,2}(:,logical(okfoil));
    mtraceprobefoil = mean(traceprobefoil,2);
    plot(1:nframes_psth,mtraceprobefoil,'color',colors_foil(i,:));
end



figure;
for i=1:nGoodDays
    subplot(2,2,1);hold on;
    tracereinftarget = sfull{i,1};
    mtracereinftarget = mean(tracereinftarget,2);
    plot(1:nframes_psth,mtracereinftarget,'color',colors_target(i,:));
    
    subplot(2,2,2);hold on;
    day = findgooddays(i);
    tracereinffoil = sfull{i,2};
    mtracereinffoil = mean(tracereinffoil,2);
    plot(1:nframes_psth,mtracereinffoil,'color',colors_foil(i,:));
    
    subplot(2,2,3);hold on;
    traceprobefoil = sprobefull{i,2};
    mtraceprobefoil = mean(traceprobefoil,2);
    plot(1:nframes_psth,mtraceprobefoil,'color',colors_target(i,:));
    
    subplot(2,2,4);hold on;
    traceprobefoil = sprobefull{i,2};
    mtraceprobefoil = mean(traceprobefoil,2);
    plot(1:nframes_psth,mtraceprobefoil,'color',colors_foil(i,:));
end

%% Plot max selectivity index

smoothfactor = 2;

[~,w] = sort(max(selectivityidxprobe,[],2));
ok = w(1:10);
ok = [86;90;91;94;97;100;101];
for c =1:3
    figure;
    for i=1:n
    
        subplot(2,n,i);hold on
        plot(Smooth(sfull{i,1}(:,ok(c)),smoothfactor));
        plot(Smooth(sfull{i,2}(:,ok(c)),smoothfactor));
        ylim([-0.15 0.14]);
        
        subplot(2,n,n+i);hold on
        plot(Smooth(sprobefull{i,1}(:,ok(c)),smoothfactor));
        plot(Smooth(sprobefull{i,2}(:,ok(c)),smoothfactor));
        ylim([-0.15 0.14]);
    end
end

%% Plot av selectivity index
savefig = false;
figure;
groups = 1:nKeptCells;
m = [[selectivityidx groups'];[selectivityidxprobe groups']];
m = sortrows(m,n+1);
m = m(:,1:end-1);
PlotColorCurves(m)

% No sorting
figure;
subplot(1,2,1);hold on;
PlotColorCurves(selectivityidx);
subplot(1,2,2);hold on;
PlotColorCurves(selectivityidxprobe);

% Max sorting
figure;
subplot(1,2,1);hold on;
[~,maxSI] = max(selectivityidx,[],2);
[~,order] = sort(maxSI);
PlotColorCurves(selectivityidx(order,:));
subplot(1,2,2);hold on;
PlotColorCurves(selectivityidxprobe(order,:));

% Mean sorting
figure;
subplot(1,2,1);hold on;
meanSI = mean(selectivityidx,2);
[~,order] = sort(meanSI);
PlotColorCurves(selectivityidx(order,:));
subplot(1,2,2);hold on;
PlotColorCurves(selectivityidxprobe(order,:));

smoothfactor = 0;
nclusters=4;
idx = kmeans(Smooth(selectivityidx,smoothfactor),nclusters);
figure;
colors = {'r','b','g','k','c','y','r','b','g','k'};
for i=1:nclusters
    subplot(nclusters/2,2,i); hold on;
    plot(Smooth(selectivityidx(idx==i,:),smoothfactor)','-','color',colors{i});
    ylim([0 1]);
    title([num2str(sum(idx==i)) ' cells']);
end    

figure;
colors = {'r','b','g','k','c','y','r','b','g','k'};
for i=1:nclusters
    subplot(nclusters/2,2,i); hold on;
    plot(Smooth(selectivityidx(idx==i,:),smoothfactor)','-','color',[0.8 0.8 0.8]);
    plot(mean(Smooth(selectivityidx(idx==i,:),smoothfactor))','-','color',colors{i});
    ylim([0 1]);
    title([num2str(sum(idx==i)) ' cells']);
end  
figure;
colors = {'r','b','g','k','c','y','r','b','g','k'};
for i=1:nclusters
    subplot(nclusters/2,2,i); hold on;
    plot(Smooth(selectivityidxprobe(idx==i,:),smoothfactor)','-','color',[0.8 0.8 0.8]);
    plot(mean(Smooth(selectivityidxprobe(idx==i,:),smoothfactor))','-','color',colors{i});
    ylim([0 1]);
    title([num2str(sum(idx==i)) ' cells']);
end 

smoothfactor = 0;
nclusters=4;
idxprobe = kmeans(Smooth(selectivityidxprobe,smoothfactor),nclusters);
figure;
colors = {'r','b','g','k','c','y','r','b','g','k'};
for i=1:nclusters
    subplot(nclusters/2,2,i); hold on;
    plot(Smooth(selectivityidxprobe(idxprobe==i,:),smoothfactor)','-','color',colors{i});
    ylim([0 1]);
    title([num2str(sum(idxprobe==i)) ' cells']);
end    
figure;
colors = {'r','b','g','k','c','y','r','b','g','k'};
for i=1:nclusters
    subplot(nclusters/2,2,i); hold on;
    plot(Smooth(selectivityidxprobe(idxprobe==i,:),smoothfactor)','-','color',[0.8 0.8 0.8]);
    plot(mean(Smooth(selectivityidxprobe(idxprobe==i,:),smoothfactor))','-','color',colors{i});
    ylim([0 1]);
    title([num2str(sum(idxprobe==i)) ' cells']);
end 


%% 
SLOW = idx==1;
FAST = idx==4;

%% Get behavioral results
if strcmp(mouse,'cd017')
    disk = 4;
elseif strcmp(mouse,'cd036')
    disk = 5;
end
results = Behavior_2P_Auditory_GNG(mouse,disk);
xdayreinf = results{1,1};
dprimereinf = results{1,2};
xdayprobe = results{2,1};
xdayprobefull = results{2,1};
dprimeprobe = results{2,2};
% figure;plot(xdayprobe,dprime,'.-');
% if strcmp(mouse,'cd017')
%     dprimeprobe(2) = [];
%     xdayprobe(2) = [];
% elseif strcmp(mouse,'cd036')
%     dprimeprobe(3:5) = [];
%     xdayprobe(3:5) = [];
% end
% figure;plot(xdayprobe,dprime,'.-');
% dprimeinterp = interp1(xdayprobe,dprimeprobe,xdayprobefull);
% figure;plot(xdayprobefull,dprimeinterp,'.-');
% xdayprobe = xdayprobefull;
% dprimeprobe = dprimeinterp;
dprimeprobe([3;4;5;7;8;9;10;11;12;13;14;15;16]) = nan;
dprimeprobe(7:end) = dprimeprobe(6);
figure;hold on;
plot(xdayreinf,dprimereinf,'k-');
plot(xdayprobe,dprimeprobe,'b--');
xdayprobe([3;4;5]) = [];
dprimeprobe([3;4;5]) = [];
dprimeinterp = interp1(xdayprobe,dprimeprobe,xdayprobefull);

%% Normalized
dprimeinterpNORM = dprimeinterp-min(dprimeinterp);
dprimeinterpNORM = dprimeinterpNORM/max(dprimeinterpNORM);
dprimereinfNORM = dprimereinf-min(dprimereinf);
dprimereinfNORM = dprimereinfNORM/max(dprimereinfNORM);
% figure;hold on;plot(dprimeinterpNORM);plot(dprimereinfNORM);
zselectivityidx = selectivityidx-repmat(min(selectivityidx,[],2),1,n);
zselectivityidx = zselectivityidx./repmat(max(zselectivityidx,[],2),1,n);
zselectivityidxprobe = selectivityidxprobe-repmat(min(selectivityidxprobe,[],2),1,n);
zselectivityidxprobe = zselectivityidxprobe./repmat(max(zselectivityidxprobe,[],2),1,n);
% c = nan(nKeptCells,2);
% p = nan(nKeptCells,2);
% for i=1:nKeptCells
%     [c(i,1),p(i,1)] = corr(zselectivityidx(i,:)',dprimeinterpNORM(findgooddays),'type','spearman');
%     [c(i,2),p(i,2)] = corr(zselectivityidx(i,:)',dprimereinfNORM(findgooddays),'type','spearman');
% end

%% All clusters
figure;hold on;
subplot(2,2,1);hold on;
plot(Smooth(selectivityidx,smoothfactor)','-','color',[0.8 0.8 0.8]);
plot(mean(Smooth(selectivityidx,smoothfactor))','k-');
ylim([0 1]);  

subplot(2,2,2);hold on;
mm = mean(Smooth(selectivityidx,smoothfactor));
normMM = (mm-min(mm))/max(mm-min(mm));
plot(normMM,'k-');
plot(dprimereinfNORM(findgooddays),'k-','linewidth',3);

call = nan(nKeptCells,2);
pall = nan(nKeptCells,2);
% for i=1:nKeptCells
%     [call(i,1),pall(i,1)] = corr(zselectivityidx(i,:)',dprimereinfNORM(findgooddays),'type','pearson');
%     [call(i,2),pall(i,2)] = corr(zselectivityidx(i,:)',dprimeinterpNORM(findgooddays),'type','pearson');
% end
for i=1:nKeptCells
    [call(i,1),pall(i,1)] = corr(selectivityidx(i,:)',dprimereinf(findgooddays),'type','pearson');
    [call(i,2),pall(i,2)] = corr(selectivityidx(i,:)',dprimeinterp(findgooddays),'type','pearson');
end
subplot(2,2,3);hold on;
g1 = linspace(0.75,1.25,nKeptCells);
plot(g1,call(:,1),'.')
xlim([0.5 1.5]);
h1 = jbtest(call(:,1));
if h1==0
    bar(1,mean(call(:,1)));
    errorbar(1,mean(call(:,1)),sem(call(:,1)));
    [~,p0] = ttest(call(:,1));
    title(['ttest, p0=' num2str(p0)]);
else
    bar(1,median(call(:,1)));
    errorbar(1,median(call(:,1)),std(call(:,1)));
    p0 = signrank(call(:,1),0,'tail','right');
    title(['W, p0=' num2str(p0)]);
    
    subplot(2,2,4);hold on;
    plot(g1(pall(:,1)<=0.05)',call(pall(:,1)<=0.05,1),'.');
    bar(1,median(call(pall(:,1)<=0.05,1)));
    errorbar(1,median(call(pall(:,1)<=0.05,1)),std(call(pall(:,1)<=0.05,1)));
end


%% Only FAST cluster

figure;hold on;
subplot(2,2,1);hold on;
plot(Smooth(selectivityidx(FAST,:),smoothfactor)','-','color',[0.8 0.8 0.8]);
plot(mean(Smooth(selectivityidx(FAST,:),smoothfactor))','b-');
ylim([0 1]);  

subplot(2,2,2);hold on;
mm = mean(Smooth(selectivityidx(FAST,:),smoothfactor));
normMM = (mm-min(mm))/max(mm-min(mm));
plot(normMM,'b-');
plot(dprimereinfNORM(findgooddays),'k-','linewidth',3);
plot(dprimeinterpNORM(findgooddays),'-','linewidth',3,'color',[0.8 0.8 0.8]);

nFAST = sum(FAST);
subplot(2,2,3);hold on;
g1 = linspace(0.75,1.25,nFAST);
g2 = linspace(1.75,2.25,nFAST);
plot([g1' g2'],[call(FAST,1) call(FAST,2)],'.')
xlim([0.5 2.5]);

h1 = jbtest(call(FAST,1));
h2 = jbtest(call(FAST,2));
if h1+h2==0
    bar([1 2],[mean(call(FAST,1)) mean(call(FAST,2))]);
    errorbar(1,mean(call(FAST,1)),sem(call(FAST,1)));
    errorbar(2,mean(call(FAST,2)),sem(call(FAST,2)));
    [~,pttest] = ttest2(call(FAST,1),call(FAST,2));
    title(['ttest, p=' num2str(pttest)]);
else
    bar([1 2],[median(call(FAST,1)) median(call(FAST,2))]);
    errorbar(1,median(call(FAST,1)),std(call(FAST,1)));
    errorbar(2,median(call(FAST,2)),std(call(FAST,2)));
    pW = ranksum(call(FAST,1),call(FAST,2));
    title(['W, p=' num2str(pW)]);
    
    subplot(2,2,4);hold on;
    plot(g1(pall(FAST,1)<=0.05)',call(FAST & (pall(:,1)<=0.05),1),'.');
    bar(1,median(call(FAST & (pall(:,1)<=0.05),1)));
    errorbar(1,median(call(FAST & (pall(:,1)<=0.05),1)),std(call(FAST & (pall(:,1)<=0.05),1)));
    plot(g2(pall(FAST,2)<=0.05)',call(FAST & (pall(:,2)<=0.05),2),'.');
    bar(2,median(call(FAST & (pall(:,2)<=0.05),2)));
    errorbar(2,median(call(FAST & (pall(:,2)<=0.05),2)),std(call(FAST & (pall(:,2)<=0.05),2)));
    
end

%% Only SLOW cluster

figure;hold on;
subplot(2,2,1);hold on;
plot(Smooth(selectivityidx(SLOW,:),smoothfactor)','-','color',[0.8 0.8 0.8]);
plot(mean(Smooth(selectivityidx(SLOW,:),smoothfactor))','b-');
ylim([0 1]);  

subplot(2,2,2);hold on;
mm = mean(Smooth(selectivityidx(SLOW,:),smoothfactor));
normMM = (mm-min(mm))/max(mm-min(mm));
plot(normMM,'b-');
plot(dprimereinfNORM(findgooddays),'k-','linewidth',3);
plot(dprimeinterpNORM(findgooddays),'-','linewidth',3,'color',[0.8 0.8 0.8]);

nSLOW = sum(SLOW);
subplot(2,2,3);hold on;
g1 = linspace(0.75,1.25,nSLOW);
g2 = linspace(1.75,2.25,nSLOW);
plot([g1' g2'],[call(SLOW,1) call(SLOW,2)],'.')
xlim([0.5 2.5]);

h1 = jbtest(call(SLOW,1));
h2 = jbtest(call(SLOW,2));

% if h1+h2==0
    bar([1 2],[median(call(SLOW,1)) median(call(SLOW,2))]);
    errorbar(1,median(call(SLOW,1)),std(call(SLOW,1)));
    errorbar(2,median(call(SLOW,2)),std(call(SLOW,2)));
    pW = ranksum(call(SLOW,1),call(SLOW,2));
    title(['W, p=' num2str(pW)]);
    ylim([-1 1]);
    
    subplot(2,2,4);hold on;
    plot(g1(pall(SLOW,1)<=0.05)',call(SLOW & (pall(:,1)<=0.05),1),'.');
    bar(1,median(call(SLOW & (pall(:,1)<=0.05),1)));
    errorbar(1,median(call(SLOW & (pall(:,1)<=0.05),1)),std(call(SLOW & (pall(:,1)<=0.05),1)));
    plot(g2(pall(SLOW,2)<=0.05)',call(SLOW & (pall(:,2)<=0.05),2),'.');
    bar(2,median(call(SLOW & (pall(:,2)<=0.05),2)));
    errorbar(2,median(call(SLOW & (pall(:,2)<=0.05),2)),std(call(SLOW & (pall(:,2)<=0.05),2)));
    
% else
%     bar([1 2],[mean(call(SLOW,1)) mean(call(SLOW,2))]);
%     errorbar(1,mean(call(SLOW,1)),sem(call(SLOW,1)));
%     errorbar(2,mean(call(SLOW,2)),sem(call(SLOW,2)));
%     [~,pttest] = ttest2(call(SLOW,1),call(SLOW,2));
%     title(['ttest, p=' num2str(pttest)]);
%     
%     subplot(2,2,4);hold on;
%     plot(g1(pall(SLOW,1)<=0.05)',call(SLOW & (pall(:,1)<=0.05),1),'.');
%     bar(1,mean(call(SLOW & (pall(:,1)<=0.05),1)));
%     errorbar(1,mean(call(SLOW & (pall(:,1)<=0.05),1)),sem(call(SLOW & (pall(:,1)<=0.05),1)));
%     plot(g2(pall(SLOW,2)<=0.05)',call(SLOW & (pall(:,2)<=0.05),2),'.');
%     bar(2,mean(call(SLOW & (pall(:,2)<=0.05),2)));
%     errorbar(2,mean(call(SLOW & (pall(:,2)<=0.05),2)),sem(call(SLOW & (pall(:,2)<=0.05),2)));
% end

%% Look at PRE tuning curves of neurons in SLOW and FAST cluster
cd(sbxpath)
names = {'TuningCurve_cd036_000_003';'TuningCurve_cd036_017_001';'TuningCurve_cd036_018_001'};
fignames = {'Pre learning';'Post learning D1';'Post learning D7'};
nameidx = 2;
name = names{nameidx};
table = load([sbxpath '/' name '/population/neuronRespTable.mat']);
table = table.neuronRespTable;
tuningMedian = load([sbxpath '/' name '/population/tuningMedian.mat']);
tuningMedian = tuningMedian.tuningMedian;
neuronIDX = load([sbxpath '/' name '/population/neuronIndex.mat']);
neuronIDX = neuronIDX.neuronIDX;
responsiveCellFlag = load([sbxpath '/' name '/population/responsiveCellFlag.mat'],'responsiveCellFlag');
responsiveCellFlag = responsiveCellFlag.responsiveCellFlag';

takenCellIDX = find(takenCells)';

totakeTC = ismember(neuronIDX,takenCellIDX);
neuronIDX_restricted = neuronIDX(totakeTC);

totakeSI = ismember(takenCellIDX,neuronIDX);
takenCellIDX_restricted = takenCellIDX(totakeSI);

% sum(totakeTC)
% sum(totakeSI)
ntotake = sum(totakeTC);




%% Plot SLOW and FAST cluster in FOV

neuronPerPlane = [0;size(signals{1},1);size(signals{1},1)];
idxSLOW = takenCellIDX(totakeSI & SLOW);
idxFAST = takenCellIDX(totakeSI & FAST);

responsiveCellFlagSelected = responsiveCellFlag(totakeTC);

[~,peakIndex] = max(tuningMedian);
peakIndex = peakIndex(totakeTC');
colormapIndex = round(linspace(1,64,17));
C = colormap('jet');
figure;
for p = 1:nPlanes
    roiName = [suite2ppath 'plane' num2str(p-1) '/' mouse '_roi' int2str(p-1) '.zip'];% CD 7/10/2019
    rois = ReadImageJROI(roiName);
    tcells = takenCellIDX_restricted((takenCellIDX_restricted > neuronPerPlane(p)) & (takenCellIDX_restricted <= neuronPerPlane(p+1)+neuronPerPlane(p)));
    tresponsiveCellFlagSelected = responsiveCellFlagSelected((takenCellIDX_restricted > neuronPerPlane(p)) & (takenCellIDX_restricted <= neuronPerPlane(p+1)+neuronPerPlane(p)));
    tidxSLOW = idxSLOW(ismember(idxSLOW,tcells));
    tidxFAST = idxFAST(ismember(idxFAST,tcells));
    tpeakIndex = peakIndex((takenCellIDX_restricted > neuronPerPlane(p)) & (takenCellIDX_restricted <= neuronPerPlane(p+1)+neuronPerPlane(p)));
    if p==2        
        tcells = tcells-neuronPerPlane(p); 
        tidxSLOW = tidxSLOW-neuronPerPlane(p); 
        tidxFAST = tidxFAST-neuronPerPlane(p); 
    end

    load([suite2ppath 'plane' num2str(p-1) '/Fall.mat'],'ops');
    meanImg = ops.meanImgE;       
    subplot(3,3,1+(p-1)*3);imagesc(meanImg);
    colormap gray;    
    for jj=1:length(tcells)
        j = tcells(jj);
        x=rois{1,j}.mnCoordinates(:,1); %freehand rois have the outlines in x-y coordinates
        y=rois{1,j}.mnCoordinates(:,2); %matlab matrices are inverted so x values are the 2nd column of mn coordinates, and y is the 1st columna
        if ismember(tcells(jj),tidxSLOW)
            patch(x,y,'g','EdgeColor','none'); 
        elseif ismember(tcells(jj),tidxFAST)
            patch(x,y,'b','EdgeColor','none'); 
        end
    end
    yStart = ops.xrange(1);
    xStart = ops.yrange(1);
    ylim([xStart+1 xStart+size(ops.max_proj,1)]);
    xlim([yStart+1 yStart+size(ops.max_proj,2)]);
    title(['Cluster index - Plane ' num2str(p)]); 
        
    subplot(3,3,2+(p-1)*3);imagesc(meanImg);
    colormap gray;
    subplot(3,3,3+(p-1)*3);imagesc(meanImg);
    colormap gray;    
    for jj=1:length(tcells)
        j = tcells(jj);
        x=rois{1,j}.mnCoordinates(:,1); %freehand rois have the outlines in x-y coordinates
        y=rois{1,j}.mnCoordinates(:,2); %matlab matrices are inverted so x values are the 2nd column of mn coordinates, and y is the 1st columna
        if ismember(tcells(jj),tidxSLOW)
            subplot(3,3,2+(p-1)*3);hold on;
            if tresponsiveCellFlagSelected(jj)
                 patch(x,y,C(colormapIndex(tpeakIndex(jj)),:),'EdgeColor','none'); 
            else 
                patch(x,y,[0.8 0.8 0.8],'EdgeColor','none');
            end        
        elseif ismember(tcells(jj),tidxFAST)
            subplot(3,3,3+(p-1)*3);hold on;
            if tresponsiveCellFlagSelected(jj)
                 patch(x,y,C(colormapIndex(tpeakIndex(jj)),:),'EdgeColor','none'); 
            else 
                patch(x,y,[0.8 0.8 0.8],'EdgeColor','none');
            end        
        end
    end
    subplot(3,3,2+(p-1)*3);hold on;
    yStart = ops.xrange(1);
    xStart = ops.yrange(1);
    ylim([xStart+1 xStart+size(ops.max_proj,1)]);
    xlim([yStart+1 yStart+size(ops.max_proj,2)]);
    title([fignames{nameidx} ' - SLOW cluster']); 
    
    subplot(3,3,3+(p-1)*3);hold on;
    yStart = ops.xrange(1);
    xStart = ops.yrange(1);
    ylim([xStart+1 xStart+size(ops.max_proj,1)]);
    xlim([yStart+1 yStart+size(ops.max_proj,2)]);
    title([fignames{nameidx} ' - FAST cluster']); 
    
end



%%
figure;hold on;
plot(xdayreinf,dprimereinf,'k-');
plot(xdayprobefull,dprimeinterp,'b--');

c = nan(nKeptCells,2);
p = nan(nKeptCells,2);
for i=1:nKeptCells
    [c(i,1),p(i,1)] = corr(selectivityidx(i,:)',dprimeinterp(findgooddays),'type','pearson');
    [c(i,2),p(i,2)] = corr(selectivityidx(i,:)',dprimereinf(findgooddays),'type','pearson');
end


%%
[cm1,pm1] = corr(mean(selectivityidx(SLOW,:))',dprimeinterp(findgooddays),'type','pearson');
[cm2,pm2] = corr(mean(selectivityidx(SLOW,:))',dprimereinf(findgooddays),'type','pearson');
cm1
pm1
cm2
pm2
[cm1,pm1] = corr(mean(selectivityidx(FAST,:))',dprimeinterp(findgooddays),'type','pearson');
[cm2,pm2] = corr(mean(selectivityidx(FAST,:))',dprimereinf(findgooddays),'type','pearson');
cm1
pm1
cm2
pm2

figure;hold on;
plot(mean(selectivityidx(FAST,:))',dprimeinterp(findgooddays),'.')
plot(mean(selectivityidx(FAST,:))',dprimereinf(findgooddays),'k.')

figure;
subplot(2,2,1);hold on;
plot(findgooddays,dprimeinterp(findgooddays),'color',[0.8 0.8 0.8]);
plot(findgooddays,dprimereinf(findgooddays),'color','k');
plot(findgooddays,mean(selectivityidx(FAST,:)),'color',colors{4});

subplot(2,2,2);hold on;
plot(findgooddays,dprimeinterp(findgooddays),'color',[0.8 0.8 0.8]);
plot(findgooddays,dprimereinf(findgooddays),'color','k');
plot(findgooddays,mean(selectivityidx(SLOW,:)),'color',colors{4});

figure;
subplot(2,2,1);hold on;
plot(findgooddays,dprimeinterpNORM(findgooddays),'color',[0.8 0.8 0.8]);
plot(findgooddays,dprimereinfNORM(findgooddays),'color','k');
mm = mean(selectivityidx(FAST,:));
plot(findgooddays,(mm-min(mm))/max(mm-min(mm)),'g');
mm = mean(selectivityidx(SLOW,:));
plot(findgooddays,(mm-min(mm))/max(mm-min(mm)),'b');

subplot(2,2,2);hold on;
plot(findgooddays,dprimeinterpNORM(findgooddays),'color',[0.8 0.8 0.8]);
plot(findgooddays,dprimereinfNORM(findgooddays),'color','k');
plot(findgooddays,mean(zselectivityidx(SLOW,:)),'r');

    
%%
figure;hold on;
subplot(2,2,1);hold on;
hist(c(:,1),-1:0.07:1);xlim([-1 1]);
PlotHVLines(0,'v','k--');
subplot(2,2,2);hold on;
hist(c(:,2),-1:0.07:1);xlim([-1 1]);
PlotHVLines(0,'v','k--');
% subplot(2,2,4);hold on;
% plot(c(:,1),c(:,2),'.');plot([-1 1],[-1 1],'k');
% xlim([-1 1]);
% ylim([-1 1]);
subplot(2,2,3);hold on;
hist(c(p(:,1)<=0.05,1),-1:0.07:1);xlim([-1 1]);
PlotHVLines(0,'v','k--');
subplot(2,2,4);hold on;
hist(c(p(:,2)<=0.05,2),-1:0.07:1);xlim([-1 1]);
PlotHVLines(0,'v','k--');

%%
figure;hold on;
subplot(2,2,1);hold on;
hist(c(SLOW,1),-1:0.07:1);xlim([-1 1]);
PlotHVLines(0,'v','k--');
title(['SLOW Corr w/ probe - med = ' num2str(median(c(SLOW,1))) ', ' num2str(sum(SLOW)) ' cells'])
subplot(2,2,2);hold on;
hist(c(SLOW,2),-1:0.07:1);xlim([-1 1]);
PlotHVLines(0,'v','k--');
title(['SLOW Corr w/ reinf - med = ' num2str(median(c(SLOW,2))) ', ' num2str(sum(SLOW)) ' cells'])
% subplot(2,2,4);hold on;
% plot(c(:,1),c(:,2),'.');plot([-1 1],[-1 1],'k');
% xlim([-1 1]);
% ylim([-1 1]);
subplot(2,2,3);hold on;
hist(c(p(:,1)<=0.05 & SLOW,1),-1:0.07:1);xlim([-1 1]);
PlotHVLines(0,'v','k--');
title(['med = ' num2str(median(c(p(:,1)<=0.05 & SLOW,1))) ', ' num2str(sum(p(:,1)<=0.05 & SLOW)) ' signif cells'])
subplot(2,2,4);hold on;
hist(c(p(:,2)<=0.05 & SLOW,2),-1:0.07:1);xlim([-1 1]);
PlotHVLines(0,'v','k--');
title(['med = ' num2str(median(c(p(:,2)<=0.05 & SLOW,2))) ', ' num2str(sum(p(:,2)<=0.05 & SLOW)) ' signif cells'])

figure;hold on;
subplot(2,2,1);hold on;
hist(c(FAST,1),-1:0.07:1);xlim([-1 1]);
PlotHVLines(0,'v','k--');
title(['FAST Corr w/ probe - med = ' num2str(median(c(FAST,1))) ', ' num2str(sum(SLOW)) ' cells'])
ylim([0 15])
subplot(2,2,2);hold on;
hist(c(FAST,2),-1:0.07:1);xlim([-1 1]);
PlotHVLines(0,'v','k--');
ylim([0 15])
title(['FAST Corr w/ reinf - med = ' num2str(median(c(FAST,2))) ', ' num2str(sum(SLOW)) ' cells'])
subplot(2,2,3);hold on;
hist(c(p(:,1)<=0.05 & FAST,1),-1:0.07:1);xlim([-1 1]);
PlotHVLines(0,'v','k--');
title(['med = ' num2str(median(c(p(:,1)<=0.05 & FAST,1))) ', ' num2str(sum(p(:,1)<=0.05 & FAST)) ' signif cells'])
subplot(2,2,4);hold on;
hist(c(p(:,2)<=0.05 & FAST,2),-1:0.07:1);xlim([-1 1]);
PlotHVLines(0,'v','k--');
title(['med = ' num2str(median(c(p(:,2)<=0.05 & FAST,2))) ', ' num2str(sum(p(:,2)<=0.05 & FAST)) ' signif cells'])


%%
figure;hold on;
subplot(2,2,1);hold on;
hist(c(SLOW | FAST,1),-1:0.07:1);xlim([-1 1]);
PlotHVLines(0,'v','k--');
title(['SLOW+FAST Corr w/ probe - med = ' num2str(median(c(SLOW | FAST,1))) ', ' num2str(sum(SLOW | FAST)) ' cells'])
subplot(2,2,2);hold on;
hist(c(SLOW | FAST,2),-1:0.07:1);xlim([-1 1]);
PlotHVLines(0,'v','k--');
title(['SLOW+FAST Corr w/ reinf - med = ' num2str(median(c(SLOW | FAST,2))) ', ' num2str(sum(SLOW | FAST)) ' cells'])
% subplot(2,2,4);hold on;
% plot(c(:,1),c(:,2),'.');plot([-1 1],[-1 1],'k');
% xlim([-1 1]);
% ylim([-1 1]);
subplot(2,2,3);hold on;
hist(c(p(:,1)<=0.05 & (SLOW | FAST),1),-1:0.07:1);xlim([-1 1]);
PlotHVLines(0,'v','k--');
title(['med = ' num2str(median(c(p(:,1)<=0.05 & (SLOW | FAST),1))) ', ' num2str(sum(p(:,1)<=0.05 & (SLOW | FAST))) ' signif cells'])
subplot(2,2,4);hold on;
hist(c(p(:,2)<=0.05 & (SLOW | FAST),2),-1:0.07:1);xlim([-1 1]);
PlotHVLines(0,'v','k--');
title(['med = ' num2str(median(c(p(:,2)<=0.05 & (SLOW | FAST),2))) ', ' num2str(sum(p(:,2)<=0.05 & (SLOW | FAST))) ' signif cells'])

%% Clean hist fig

figure;hold on;
subplot(2,2,1);hold on;
hist(c(SLOW,2),-1:0.07:1);xlim([-1 1]);
PlotHVLines(0,'v','k--');
title(['SLOW Corr w/ Reinf - med = ' num2str(median(c(SLOW,2))) ', ' num2str(sum(SLOW)) ' cells'])
PlotHVLines(median(c(SLOW,2)),'v','r-');
subplot(2,2,2);hold on;
hist(c(FAST,1),-1:0.07:1);xlim([-1 1]);
PlotHVLines(0,'v','k--');
PlotHVLines(median(c(FAST,1)),'v','r-');
title(['FAST Corr w/ Probe - med = ' num2str(median(c(FAST,1))) ', ' num2str(sum(FAST)) ' cells']);
subplot(2,2,3);hold on;
hist(c(FAST | SLOW,2),-1:0.07:1);xlim([-1 1]);
PlotHVLines(0,'v','k--');
PlotHVLines(median(c(FAST | SLOW,2)),'v','r-');
title(['FAST Corr w/ Reinf - med = ' num2str(median(c(FAST | SLOW,2))) ', ' num2str(sum(FAST | SLOW)) ' cells']);




%%


% Plot cdf corr fast/reinf & probe AND slow/reinf & probe
figure;
bins = -1:0.07:1;
subplot(2,2,1);hold on;
[h,ht] = hist(c(SLOW,1),bins);
s = CumSum(h/sum(h));
% toplot = ht>=-0.05 & ht<=0.05;
% plot(ht(toplot),s(toplot),'linewidth',2);
plot(ht,s,'linewidth',2);
PlotHVLines(0,'v','k:');
[h,ht] = hist(c(SLOW,2),bins);
s = CumSum(h/sum(h));
plot(ht,s,'linewidth',2,'color','k');
[~,pKS] = kstest2(c(SLOW,1),c(SLOW,2));
title(['SLOW, pKS=' num2str(pKS)]);

subplot(2,2,2);hold on;
[h,ht] = hist(c(FAST,1),bins);
s = CumSum(h/sum(h));
plot(ht,s,'linewidth',2);
PlotHVLines(0,'v','k:');
[h,ht] = hist(c(FAST,2),bins);
s = CumSum(h/sum(h));
plot(ht,s,'linewidth',2,'color','k');
[~,pKS] = kstest2(c(FAST,1),c(FAST,2));
title(['FAST, pKS=' num2str(pKS)]);


subplot(2,2,3);hold on;
[h,ht] = hist(c(SLOW,1),bins);
s = CumSum(h/sum(h));
plot(ht,s,'b-','linewidth',2);
PlotHVLines(0,'v','k:');
[h,ht] = hist(c(FAST,1),bins);
s = CumSum(h/sum(h));
plot(ht,s,'b--','linewidth',2);
[~,pKS] = kstest2(c(SLOW,1),c(FAST,1));
title(['SLOW & FAST, Probe, pKS=' num2str(pKS)]);

subplot(2,2,4);hold on;
[h,ht] = hist(c(SLOW,2),bins);
s = CumSum(h/sum(h));
plot(ht,s,'k-','linewidth',2);
PlotHVLines(0,'v','k:');
[h,ht] = hist(c(FAST,2),bins);
s = CumSum(h/sum(h));
plot(ht,s,'k--','linewidth',2);
[~,pKS] = kstest2(c(SLOW,2),c(FAST,2));
title(['SLOW & FAST, Reinf, pKS=' num2str(pKS)]);


%% Plot cdf corr FAST+SLOW with Reinf and probe AND separated
figure;
bins = -1:0.07:1;
subplot(2,2,1);hold on;
[h,ht] = hist(c(SLOW,1),bins);
s = CumSum(h/sum(h));
% toplot = ht>=-0.05 & ht<=0.05;
% plot(ht(toplot),s(toplot),'linewidth',2);
plot(ht,s,'linewidth',2);
PlotHVLines(0,'v','k:');
[h,ht] = hist(c(SLOW,2),bins);
s = CumSum(h/sum(h));
plot(ht,s,'linewidth',2,'color','k');
[~,pKS] = kstest2(c(SLOW,1),c(SLOW,2));
title(['SLOW, reinf and probe, pKS=' num2str(pKS)]);
[h,ht] = hist(c(SLOW | FAST,2),bins);
s = CumSum(h/sum(h));
plot(ht,s,'linewidth',2,'color','r');

subplot(2,2,2);hold on;
[h,ht] = hist(c(FAST,1),bins);
s = CumSum(h/sum(h));
plot(ht,s,'linewidth',2);
PlotHVLines(0,'v','k:');
[h,ht] = hist(c(FAST,2),bins);
s = CumSum(h/sum(h));
plot(ht,s,'linewidth',2,'color','k');
[~,pKS] = kstest2(c(FAST,1),c(FAST,2));
title(['FAST, reinf and probe, pKS=' num2str(pKS)]);
[h,ht] = hist(c(SLOW | FAST,2),bins);
s = CumSum(h/sum(h));
plot(ht,s,'linewidth',2,'color','r');

subplot(2,2,3);hold on;
[h,ht] = hist(c(SLOW,1),bins);
s = CumSum(h/sum(h));
plot(ht,s,'b-','linewidth',2);
PlotHVLines(0,'v','k:');
[h,ht] = hist(c(FAST,1),bins);
s = CumSum(h/sum(h));
plot(ht,s,'b--','linewidth',2);
[~,pKS] = kstest2(c(SLOW,1),c(FAST,1));
title(['SLOW & FAST, Probe, pKS=' num2str(pKS)]);
[h,ht] = hist(c(SLOW | FAST,2),bins);
s = CumSum(h/sum(h));
plot(ht,s,'linewidth',2,'color','r');

subplot(2,2,4);hold on;
[h,ht] = hist(c(SLOW,2),bins);
s = CumSum(h/sum(h));
plot(ht,s,'k-','linewidth',2);
PlotHVLines(0,'v','k:');
[h,ht] = hist(c(FAST,2),bins);
s = CumSum(h/sum(h));
plot(ht,s,'k--','linewidth',2);
[~,pKS] = kstest2(c(SLOW,2),c(FAST,2));
title(['SLOW & FAST, Reinf, pKS=' num2str(pKS)]);
[h,ht] = hist(c(SLOW | FAST,2),bins);
s = CumSum(h/sum(h));
plot(ht,s,'linewidth',2,'color','r');

[~,pKS] = kstest2(c(SLOW,2),c(SLOW | FAST,2));
%%


sum([p(:,1)<=0.05 p(:,2)<=0.05],2)

cellIncluster3 = find(SLOW);
c(cellIncluster3(1),1)
c(cellIncluster3(1),2)

cellIncluster4 = find(FAST);
c(cellIncluster4(1),1)
c(cellIncluster4(1),2)

figure;
subplot(2,2,1);hold on;
plot([1 2],[c(SLOW,1) c(SLOW,2)],'.-');
xlim([0.5 2.5]);ylim([-1 1]);
mean(c(SLOW,1))
mean(c(SLOW,2))
subplot(2,2,2);hold on;
plot([1 2],[c(FAST,1) c(FAST,2)],'.-');
xlim([0.5 2.5]);ylim([-1 1]);
mean(c(FAST,1))
mean(c(FAST,2))


sum(SLOW)
sum(idx==2)
sum(idx==3)
sum(FAST)

figure;hold on;
for i=1:nclusters
    subplot(2,2,i);hold on;
    plot(c(idx==i,1),c(idx==i,2),'.','color',colors{i});
    plot([-1 1],[-1 1],'k');
    xlim([-1 1]);
    ylim([-1 1]);
end

figure;
subplot(2,1,1);hold on;
toplot = selectivityidx;
plot(repmat(dprimereinf(daytoplot)',size(selectivityidx,1),1)',toplot','.-');
xlim([0 3]); %ylim([-0.4 0.4]);
subplot(2,1,2);hold on;
toplotprobe = selectivityidxprobe;
plot(repmat(dprimeprobe(daytoplot)',size(selectivityidxprobe,1),1)',toplotprobe','.-');
xlim([0 3]); %ylim([-0.4 0.4]);





% z-score selectivity
figure;
subplot(1,2,1);hold on;
zselectivityidx = selectivityidx-repmat(min(selectivityidx,[],2),1,n);
zselectivityidx = zselectivityidx./repmat(max(zselectivityidx,[],2),1,n);
zselectivityidxprobe = selectivityidxprobe-repmat(min(selectivityidxprobe,[],2),1,n);
zselectivityidxprobe = zselectivityidxprobe./repmat(max(zselectivityidxprobe,[],2),1,n);
[~,maxSI] = max(zselectivityidx,[],2);
[~,order] = sort(maxSI);
PlotColorCurves(zselectivityidx(order,:));
subplot(1,2,2);hold on;
PlotColorCurves(zselectivityidxprobe(order,:));


for i=118:nKeptCells
%     m = [selectivityidx(i,:);selectivityidxprobe(i,:)];
%     PlotColorCurves(m);
    f=figure;hold on;
    plot(selectivityidx(i,:),'k.-');
    plot(selectivityidxprobe(i,:),'b.--');
    title(['cell ' num2str(i)]);
    ylim([0 1]);
    pause();
    close(f);
end

for c=[65 68 93 119 120 147 123 269 394 438 551 554 586]
%     m = [selectivityidx(i,:);selectivityidxprobe(i,:)];
%     PlotColorCurves(m);
    f=figure;
    subplot(3,1,1);hold on;
    plot(selectivityidx(c,:),'k.-');
    plot(selectivityidxprobe(c,:),'b.--');
    title(['cell ' num2str(c)]);
    ylim([0 1]);
    ylabel('Selectivity Index');
    xlabel('Days')
    
    for i=1:n    
        subplot(3,n,n+i);hold on
        plot(Smooth(sfull{i,1}(:,c),smoothfactor),'g');
        plot(Smooth(sfull{i,2}(:,c),smoothfactor),'r');
        ylim([-0.15 2]);
        xlabel('Days')
        ylabel('DF/F');
%         axis off
        
        subplot(3,n,2*n+i);hold on
        plot(Smooth(sprobefull{i,1}(:,c),smoothfactor),'g');
        plot(Smooth(sprobefull{i,2}(:,c),smoothfactor),'r');
        ylim([-0.15 2]);
        xlabel('Days');
        ylabel('DF/F');
%         axis off
    end
    if savefig
        if strcmp(mouse,'cd036')
            saveas(f,['T:\LabData5\cd036\analysis\' mouse '_SelectivityTrajectory_FAST_Cell' num2str(c) '.pdf']);
        elseif strcmp(mouse,'cd017')
            saveas(f,['U:\LabData4\celine\cd017\analysis\' mouse '_SelectivityTrajectory_FAST_Cell' num2str(c) '.pdf']);
        end
    end
    close(f);
end


for c=[27 29 30 31 32 33 35 37 71 72 46 82 89 122]
%     m = [selectivityidx(i,:);selectivityidxprobe(i,:)];
%     PlotColorCurves(m);
    f=figure;
    subplot(3,1,1);hold on;
    plot(selectivityidx(c,:),'k.-');
    plot(selectivityidxprobe(c,:),'b.--');
    title(['cell ' num2str(c)]);
    ylim([0 1]);
    ylabel('Selectivity Index');
    xlabel('Days')
       
    for i=1:n    
        subplot(3,n,n+i);hold on
        plot(Smooth(sfull{i,1}(:,c),smoothfactor),'g');
        plot(Smooth(sfull{i,2}(:,c),smoothfactor),'r');
        ylim([-0.15 0.5]);
        xlabel('Days')
        ylabel('DF/F');
%         axis off
        
        subplot(3,n,2*n+i);hold on
        plot(Smooth(sprobefull{i,1}(:,c),smoothfactor),'g');
        plot(Smooth(sprobefull{i,2}(:,c),smoothfactor),'r');
        ylim([-0.15 0.5]);
        xlabel('Days');
        ylabel('DF/F');
%         axis off
    end
    if savefig
        if strcmp(mouse,'cd036')
            saveas(f,['T:\LabData5\cd036\analysis\' mouse '_SelectivityTrajectory_SLOW_Cell' num2str(c) '.pdf']);
        elseif strcmp(mouse,'cd017')
            saveas(f,['U:\LabData4\celine\cd017\analysis\' mouse '_SelectivityTrajectory_SLOW_Cell' num2str(c) '.pdf']);
        end
    end
%     pause();
    
    close(f);
end



for c=[68 147 123 394 438 551 554 586]
%     figure;
    for i=1:n    
        subplot(3,n,n+i);hold on
        plot(Smooth(sfull{i,1}(:,c),smoothfactor));
        plot(Smooth(sfull{i,2}(:,c),smoothfactor));
        ylim([-0.15 2]);
        
        subplot(3,n,2*n+i);hold on
        plot(Smooth(sprobefull{i,1}(:,c),smoothfactor));
        plot(Smooth(sprobefull{i,2}(:,c),smoothfactor));
        ylim([-0.15 2]);
    end
end



for i=22:nKeptCells
    if selectivityidx(i,2)-selectivityidxprobe(i,2)>-0.2, continue, end
%     m = [selectivityidx(i,:);selectivityidxprobe(i,:)];
%     PlotColorCurves(m);
    f=figure;hold on;
    plot(selectivityidx(i,:),'k.-');
    plot(selectivityidxprobe(i,:),'b.--');
    plot(selectivityidx(i,:)-selectivityidxprobe(i,:),'r-');
    title(['cell ' num2str(i)]);
    ylim([-1 1]);
    pause();
    close(f);
end

f=figure;
for i=1:n
    subplot(4,5,i); hold on;
    c = corr(selectivityidx(:,i),selectivityidxprobe(:,i));
    PlotColorMap(c);
end




