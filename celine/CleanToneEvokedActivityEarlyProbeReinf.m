function CleanToneEvokedActivityEarlyProbeReinf()

mice = {'cd017','cd036'};
RECAP_TONEEVOKED = cell(length(mice),6);
thismouse = 1;
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

%% MEDIAN TRACE PER TONE (ALL PLANES) PER DAY
ndays = 5;

pretone = 1; % for PSTH tone, in s
posttone = 4; % for PSTH tone, in s
nframes_psth = pretone*round(acq) + posttone*round(acq);


ms = nan(nframes_psth,nDays_behav,2);
ss = nan(nframes_psth,nDays_behav,2);

msprobe = nan(nframes_psth,nDays_behav,2);
ssprobe = nan(nframes_psth,nDays_behav,2);

sfull = cell(nDays_behav,2);
ssfull = cell(nDays_behav,2);

sprobefull = cell(nDays_behav,2);
ssprobefull = cell(nDays_behav,2);

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
                
        msprobe(:,j,k) = mean(squeeze(median(s_probe(:,:,takenCells),2)),2);
        ssprobe(:,j,k) = sem(squeeze(median(s_probe(:,:,takenCells),2))')';
        
        sprobefull{j,k} = squeeze(median(s_probe(:,:,takenCells),2));
        ssprobefull{j,k} = squeeze(std(s_probe(:,:,takenCells),0,2));  
    end
end
%%
n = 5;

figure;
colors_target = [0*ones(n,1) linspace(0.2,0.8,n)' 0*ones(n,1)];
colors_foil = [linspace(0.2,0.8,n)' 0*ones(n,1) 0*ones(n,1)];   
xmax = nframes_psth-25;
tonenames = {'Target','Foil'};
for j=1:2
    subplot(2,2,j+(j-1));hold on;
    for ff=1:n
        i = findgooddays(ff);
        if j==1
            shadedErrorBar(1:xmax,ms(1:xmax,i,j),ss(1:xmax,i,j),{'color',colors_target(ff,:)});
        else
            shadedErrorBar(1:xmax,ms(1:xmax,i,j),ss(1:xmax,i,j),{'color',colors_foil(ff,:)});
        end
    end
    ylim([-0.04 0.1]);
    PlotHVLines(pretone*round(acq),'v','k','linewidth',1);   
    PlotHVLines(0,'h','k:','linewidth',1);
    title(tonenames{j});    
    
    subplot(2,2,j+(j-1)+1);hold on;
    for ff=1:n
        i = findgooddays(ff);
        tms = ms(1:xmax,i,j)-mean(ms(1:pretone*round(acq),i,j));
        if j==1
            shadedErrorBar(1:xmax,tms,ss(1:xmax,i,j),{'color',colors_target(ff,:)});
        else
            shadedErrorBar(1:xmax,tms,ss(1:xmax,i,j),{'color',colors_foil(ff,:)});
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
    for ff=1:n
        i = findgooddays(ff);
        if j==1
            shadedErrorBar(1:xmax,msprobe(1:xmax,i,j),ssprobe(1:xmax,i,j),{'color',colors_target(ff,:)});
        else
            shadedErrorBar(1:xmax,msprobe(1:xmax,i,j),ssprobe(1:xmax,i,j),{'color',colors_foil(ff,:)});
        end
    end
    ylim([-0.04 0.1]);
    PlotHVLines(pretone*round(acq),'v','k','linewidth',1);   
    PlotHVLines(0,'h','k:','linewidth',1);
    title(tonenames{j});    
    
    subplot(2,2,j+(j-1)+1);hold on;
    for ff=1:n
        i = findgooddays(ff);
        tms = msprobe(1:xmax,i,j)-mean(msprobe(1:pretone*round(acq),i,j));
        if j==1
            shadedErrorBar(1:xmax,tms,ssprobe(1:xmax,i,j),{'color',colors_target(ff,:)});
        else
            shadedErrorBar(1:xmax,tms,ssprobe(1:xmax,i,j),{'color',colors_foil(ff,:)});
        end
    end
    ylim([-0.04 0.1]);
    PlotHVLines(pretone*round(acq),'v','k','linewidth',1);   
    PlotHVLines(0,'h','k:','linewidth',1);
    title([contextNames{2} ' - ' tonenames{j}]);
end
%% Same but as color curves
n = 5;

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
ndays = 5;
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
    daytoplot = [2,4,5];
elseif strcmp(mouse,'cd036')
    daytoplot = [1,2,5];
end
savefig = false;

% Amplitude and latency
pretone = 1; % for PSTH tone, in s
posttone = 4; % for PSTH tone, in s
xrange = [pretone pretone*round(acq)+posttone*round(acq)];

cutoffs = [-0.03 0.05];

% ORDERED BY TARGET for Target and by FOIL for foil
fig=figure;
smoothfactor = [0 1];
ndays = length(daytoplot);
for ff=1:ndays
    i = daytoplot(ff);

    subplot(ndays,2,ff+(ff-1));hold on;
    targettraces = sfull{i,1};
    if ff==1
        [~,ordertarget] = sort(mean(targettraces(15:25,:)));
    end
    PlotColorCurves(Smooth(targettraces(:,ordertarget)',smoothfactor),xrange,'cutoffs',cutoffs);
    PlotHVLines(pretone*round(acq)+1,'v','w');
    
    subplot(ndays,2,ff+(ff-1)+1);hold on;
    foiltraces = sfull{i,2};
    if ff==1
        [~,orderfoil] = sort(mean(foiltraces(15:25,:)));
    end
    
    PlotColorCurves(Smooth(foiltraces(:,orderfoil)',smoothfactor),xrange,'cutoffs',cutoffs);   
    PlotHVLines(pretone*round(acq)+1,'v','w');
    title([contextNames{1} ' - Day ' num2str(i)]);
end
if savefig
    cd 'E:\KishoreLab\celine\matlab\Analysis_Exci_Project\AnalysisSfN'
    print(fig,'-dpdf',[mouse '_IndivTargetFoilReinfSelectedDays_OrderedByTargetAndFoil.pdf']);
    close(fig);
end

fig = figure;
for ff=1:ndays
    i = daytoplot(ff);    
    
    subplot(ndays,2,ff+(ff-1));hold on;
    targettraces = sprobefull{i,1};  
    PlotColorCurves(Smooth(targettraces(:,ordertarget)',smoothfactor),xrange,'cutoffs',cutoffs);
    PlotHVLines(pretone*round(acq)+1,'v','w');
    
    subplot(ndays,2,ff+(ff-1)+1);hold on;
    foiltraces = sprobefull{i,2};
    PlotColorCurves(Smooth(foiltraces(:,orderfoil)',smoothfactor),xrange,'cutoffs',cutoffs);    
    PlotHVLines(pretone*round(acq)+1,'v','w');
    title([contextNames{2} ' - Day ' num2str(i)]);
end
if savefig
    cd 'E:\KishoreLab\celine\matlab\Analysis_Exci_Project\AnalysisSfN'
    export_fig(fig,[mouse '_IndivTargetFoilProbeSelectedDays_OrderedByTargetAndFoil.pdf']);
    close(fig);
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
    
% Get selectivity for each neurons
peakwind = 17:26;
selectivityidx = nan(nKeptCells,ndays);
selectivityidxprobe = nan(nKeptCells,ndays);
for ff=1:ndays
    i = daytoplot(ff);
    targettraces = sfull{i,1};
    foiltraces = sfull{i,2};
    
    selectivityidx(:,ff) = abs((max(targettraces(peakwind,:))-max(foiltraces(peakwind,:)))...
        ./sum([max(abs(targettraces(peakwind,:)));max(abs(foiltraces(peakwind,:)))]));
    
    targettracesprobe = sprobefull{i,1};
    foiltracesprobe = sprobefull{i,2};    
    
    selectivityidxprobe(:,ff) = abs((max(targettracesprobe(peakwind,:))-max(foiltracesprobe(peakwind,:)))...
        ./sum([max(abs(targettracesprobe(peakwind,:)));max(abs(foiltracesprobe(peakwind,:)))]));
end

% Get behavioral results
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
if strcmp(mouse,'cd017')
    dprimeprobe(2) = [];
    xdayprobe(2) = [];
elseif strcmp(mouse,'cd036')
    dprimeprobe(3:5) = [];
    xdayprobe(3:5) = [];
end
% figure;plot(xdayprobe,dprime,'.-');
dprimeinterp = interp1(xdayprobe,dprimeprobe,xdayprobefull);
% figure;plot(xdayprobefull,dprimeinterp,'.-');
xdayprobe = xdayprobefull;
dprimeprobe = dprimeinterp;


figure;
subplot(2,1,1);hold on;
toplot = selectivityidx;
plot(repmat(dprimereinf(daytoplot)',size(selectivityidx,1),1)',toplot','.-');
xlim([0 3]); %ylim([-0.4 0.4]);
subplot(2,1,2);hold on;
toplotprobe = selectivityidxprobe;
plot(repmat(dprimeprobe(daytoplot)',size(selectivityidxprobe,1),1)',toplotprobe','.-');
xlim([0 3]); %ylim([-0.4 0.4]);

figure;
subplot(2,1,1);hold on;
toplot = selectivityidx;
plot(repmat([1 2 3],size(selectivityidx,1),1)',toplot','.-');
toplotprobe = selectivityidxprobe;
plot(repmat([4 5 6],size(selectivityidxprobe,1),1)',toplotprobe','.-');
xlim([0 7]);
PlotHVLines(0,'h','k:');

subplot(2,1,2);hold on;
groups = [1:6];
boxplot([toplot toplotprobe],groups)
% ylim([-0.05 0.05]);
PlotHVLines(0,'h','k:');

% plot reinforced vs probe for each day
figure;hold on;
for i=1:ndays
    x = [i+2*(i-1) i+2*(i-1)+1];
    toplot = selectivityidx(:,i);
    toplotprobe = selectivityidxprobe(:,i);
    plot(repmat(x,size(toplot,1),1)',[toplot toplotprobe]','.-');
end
xlim([0 9]);ylabel('Selectivity Idx');

figure;
toplot = selectivityidx(:,2);
toplotprobe = selectivityidxprobe(:,2);
barwitherrnn([toplot;toplotprobe],[ones(size(toplot,1),1);2*ones(size(toplot,1),1)]);
[h0,p] = lillietest(toplot);
[~,p] = lillietest(toplotprobe);
p = signrank(toplot,toplotprobe);
figure;plot(toplot,toplotprobe,'.');

smoothfactor = 0;
toplot = selectivityidx(:,2);
toplotprobe = selectivityidxprobe(:,2);
deltaselectivityidx = [toplot toplotprobe];
nclusters=4;
idx = kmeans(Smooth(deltaselectivityidx,smoothfactor),nclusters);
x = [1 2];
figure;
subplot(1,2,1);hold on;
plot(x,Smooth(deltaselectivityidx(idx==1,:),smoothfactor)','r.-');
hold on;plot(x,Smooth(deltaselectivityidx(idx==2,:),smoothfactor)','b.-');
hold on;plot(x,Smooth(deltaselectivityidx(idx==3,:),smoothfactor)','g.-');
hold on;plot(x,Smooth(deltaselectivityidx(idx==4,:),smoothfactor)','k.-');
xlim([0 3]);

mgp = nan(nclusters,2);
semgrp = nan(nclusters,2);
colors = {'r','b','g','k'};
subplot(1,2,2);hold on;
for i=1:nclusters
    mgp(i,:) = mean(Smooth(deltaselectivityidx(idx==i,:),smoothfactor));
    semgrp(i,:) = sem(Smooth(deltaselectivityidx(idx==i,:),smoothfactor));
    shadedErrorBar(x,mgp(i,:),semgrp(i,:),{'color',colors{i}});
    ylim([0 1]);
end
PlotHVLines(0,'h','k:');



