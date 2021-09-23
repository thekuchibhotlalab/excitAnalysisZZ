function CleanToneEvokedActivity()

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
% nb of cells
H = [];
for p=1:nPlanes
    here = ishere{p};
    here(isnan(here) | here~=1) = 0;
    H = [H;here]; % neurons x days
end
nTotCells = sum(sum(H==1,2)~=0);

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

nCellSignReinf = sum(any(respp(:,:,1),2));
nCellSignProbe = sum(any(respp(:,:,2),2));

nCellSignInBothCtxs = sum(any(respp(:,:,1),2) & any(respp(:,:,2),2));
nCellSignReinfOnly = sum(any(respp(:,:,1),2) & ~any(respp(:,:,2),2));
nCellSignProbeOnly = sum(~any(respp(:,:,1),2) & any(respp(:,:,2),2));

if nCellSignProbe~=nCellSignInBothCtxs+nCellSignProbeOnly, error('Miscalculation nb cells Probe'); end
if nCellSignReinf~=nCellSignInBothCtxs+nCellSignReinfOnly, error('Miscalculation nb cells Reinfoced'); end

figure;pie([nCellSignInBothCtxs;nCellSignReinfOnly;nCellSignProbeOnly;nTotCells-nCellSignReinf],{'Both','Reinf only','Probe only','None'});

nCellSignAllDays = sum(sum(respp(:,:,1),2)==nDays_behav);
signifdays = nan(nDays_behav,2);
for i=1:nDays_behav
    for ct=1:2
        signifdays(i,ct) = sum(sum(respp(:,:,ct),2)==i);
    end
end

signifdayscum = nan(nDays_behav,2);
for i=1:nDays_behav
    for ct=1:2
        signifdayscum(i,ct) = sum(sum(respp(:,:,ct),2)>=i);
    end
end

% Evolution of signif cells across days
figure;
subplot(3,1,1); hold on;
plot(1:nDays_behav,nansum(respp(:,:,1))./sum(~isnan(respp(:,:,1))),'k.-');
plot(1:nDays_behav,nansum(respp(:,:,2))./sum(~isnan(respp(:,:,2))),'b.-');
ylim([0 1]);
PlotHVLines(0.5,'h','k:');
xlabel('Days');ylabel('Cell proportion');
title('Reinf and Probe - Tone responsive');

subplot(3,1,2); hold on;
target_reinf = squeeze(selecc(:,:,1,1)); % neurons x days
foil_reinf = squeeze(selecc(:,:,2,1));
tone_reinf = logical(squeeze(sum(squeeze(selecc(:,:,:,1)),3))); % neurons x days

target_probe = squeeze(selecc(:,:,1,2)); 
foil_probe = squeeze(selecc(:,:,2,2));
tone_probe = logical(squeeze(sum(squeeze(selecc(:,:,:,2)),3))); % neurons x days

nCellPerDay = sum(~isnan(respp(:,:,1)));

plot(1:nDays_behav,nansum(target_reinf)./nCellPerDay,'g.-');
plot(1:nDays_behav,nansum(foil_reinf)./nCellPerDay,'r.-');
plot(1:nDays_behav,nansum(respp(:,:,1))./nCellPerDay,'k.-');
ylim([0 1]);
PlotHVLines(0.5,'h','k:');
xlabel('Days');ylabel('Cell proportion');
title('Reinf - Tone, Target, Foil responsive');

subplot(3,1,3); hold on;
plot(1:nDays_behav,nansum(target_probe)./nCellPerDay,'g.-');
plot(1:nDays_behav,nansum(foil_probe)./nCellPerDay,'r.-');
plot(1:nDays_behav,nansum(respp(:,:,2))./nCellPerDay,'b.-');
ylim([0 1]);
PlotHVLines(0.5,'h','k:');
xlabel('Days');ylabel('Cell proportion');
title('Probe - Tone, Target, Foil responsive');

%% Selectivity, signif cells only
pretone = 1; % for PSTH tone, in s
posttone = 4; % for PSTH tone, in s
nframes_psth = pretone*round(acq) + posttone*round(acq);
smoothfactor = 1;
peakwindinterval = [17 26];
peakwind = 17:26;

showfig = false;
selectivity = cell(nPlanes,1);
signselectivity = cell(nPlanes,1);
for p=1:nPlanes
    signals_thisplane = signals{p};
    nCells = size(signals_thisplane,1);
    
    selectivity{p} = nan(nCells,nDays_behav,2);
    signselectivity{p} = nan(nCells,nDays_behav,2);
    for i=1:nCells
        for j=1:nDays_behav            
            this_day = matrix(:,DAY)==dayss_behav(j)+1; 
            trace = signals_thisplane(i,:,2);
            
            for ct = 1:2
                
                here = ishere{p}(i,dayss_behav(j)+1);
                here(isnan(here) | here~=1) = 0;
                here = logical(here);
                
                yep = toneresponsive{p}(i,j,ct);
                yep(isnan(yep)) = 0;
                here = logical(here) & logical(yep);
                if ~here, continue, end
                          
                ms = nan(nframes_psth,2);
                for k=1:2
                    ok = this_day & ctx(:,ct) & ks2(:,k);
                    l = sum(ok);            
                    matidx = repelem(0:nframes_psth-1,l,1);
                    idx = (matidx+repelem(matrix(ok,TONEF(p))- pretone*round(acq),1,nframes_psth))';
                    s = trace(:,idx(:))';
                    s = reshape(s,nframes_psth,l);
                    ms(:,k) = median(s,2);
                end
                [~,wm] = max([max(abs(ms(peakwind,1)));max(abs(ms(peakwind,2)))]);
                [~,wval] = max(abs(ms(peakwind,wm)));
                mms = ms(peakwind,wm);
                val = sign(mms(wval));
                selectivity{p}(i,j,ct) = max(abs(ms(peakwind,1)))-max(abs(ms(peakwind,2)));
                signselectivity{p}(i,j,ct) = val;
                
                if showfig
                    figure;hold on;
                    plot(Smooth(ms(:,1),smoothfactor),'g');plot(Smooth(ms(:,2),smoothfactor),'r');
                    PlotHVLines(pretone*round(acq),'v','k');
                    PlotIntervals(peakwindinterval);
                    title(['Cell ' num2str(i) ', ctx = ' num2str(ct) ', sign = ' num2str(val) ', idx = ' num2str(selectivity{p}(i,j,ct))]);
                    pause();
                    close();
                end
            end
        end
    end
end
%%
selecc = [];signselecc = [];
for p=1:nPlanes
    selecc = [selecc; selectivity{p}];
    signselecc = [signselecc; signselectivity{p}];
end

figure;
subplot(2,1,1); hold on;
[~,order] =  sortrows(squeeze(selecc(:,:,1)),1:nDays_behav);
PlotColorCurves(squeeze(selecc(order,:,1)),[1 nDays_behav],'cutoffs',[-0.05 0.05]);

subplot(2,1,2); hold on;
PlotColorCurves(squeeze(selecc(:,:,2)),[1 nDays_behav],'cutoffs',[-0.05 0.05]);

figure;
subplot(2,1,1); hold on;
transformsingselec = signselecc;
transformsingselec(isnan(transformsingselec)) = 0;
[~,order] = sort(sum(squeeze(transformsingselec(:,:,1)),2));
PlotColorCurves(squeeze(transformsingselec(order,:,1)));

subplot(2,1,2); hold on;
PlotColorCurves(squeeze(transformsingselec(order,:,2)));

figure;
mseleccperday = squeeze(nanmean(selecc)); % days x contexts
subplot(3,1,1); hold on;
plot(mseleccperday(:,1),'k.-');
plot(mseleccperday(:,2),'b.-');

subplot(3,1,2); hold on;
plot(squeeze(selecc(:,:,1))','.-');

subplot(3,1,3); hold on;
plot(squeeze(selecc(:,:,2))','.-');
%% MEDIAN TRACE PER TONE (ALL PLANES) PER DAY

pretone = 1; % for PSTH tone, in s
posttone = 4; % for PSTH tone, in s
nframes_psth = pretone*round(acq) + posttone*round(acq);

responsive = true;

ms = nan(nframes_psth,nDays_behav,2);
ss = nan(nframes_psth,nDays_behav,2);

msprobe = nan(nframes_psth,nDays_behav,2);
ssprobe = nan(nframes_psth,nDays_behav,2);

sfull = cell(nDays_behav,2);
ssfull = cell(nDays_behav,2);

colors = {'g','r'};
figure;hold on;

for j=1:nDays_behav            
    this_day = matrix(:,DAY)==dayss_behav(j)+1;    
    for k=1:2        
        for p=1:nPlanes % loop through planes
            signals_thisplane = signals{p};
            trace = signals_thisplane(:,:,2);

            if p==2, prevncells = nCellsHere; end % continue matrix

            here = ishere{p}(:,dayss_behav(j)+1);
            here(isnan(here) | here~=1) = 0;  
            
            if responsive
%                 yep = sum(toneresponsive{p}(:,j,:),3); % keep according to BOTH contexts
                yep = toneresponsive{p}(:,j,1); % keep according to REINFORCED contexts
                yep(isnan(yep)) = 0;
                here = logical(here) & logical(yep);
            end
            
            nCellsHere = sum(here);

            ok = this_day & ctx(:,1) & ~nop & ks2(:,k); % all tone (i.e. all trials)
            l = sum(ok); % nb of (these) tones
            matidx = repelem(0:nframes_psth-1,l,1);
            idx = (matidx+repelem(matrix(ok,TONEF(p))- pretone*round(acq),1,nframes_psth))';

            ok_probe = this_day & ctx(:,2) & ~nop & ks2(:,k); % all tone (i.e. all trials)
            l_probe = sum(ok_probe); % nb of (these) tones
            matidx_probe = repelem(0:nframes_psth-1,l_probe,1);
            idx_probe = (matidx_probe+repelem(matrix(ok_probe,TONEF(p))- pretone*round(acq),1,nframes_psth))';

            if p==1
                s = trace(logical(here),idx(:))';    
                s = reshape(s,nframes_psth,l,nCellsHere);     

                s_probe = trace(logical(here),idx_probe(:))';    
                s_probe = reshape(s_probe,nframes_psth,l_probe,nCellsHere);     
            else
                s2 = trace(logical(here),idx(:))'; 
                s2 = reshape(s2,nframes_psth,l,nCellsHere); 
                s(:,:,prevncells+1:prevncells+nCellsHere) = s2;

                s2_probe = trace(logical(here),idx_probe(:))'; 
                s2_probe = reshape(s2_probe,nframes_psth,l_probe,nCellsHere); 
                s_probe(:,:,prevncells+1:prevncells+nCellsHere) = s2_probe;            
            end
        end

        ms(:,j,k) = mean(squeeze(median(s,2)),2);
        ss(:,j,k) = sem(squeeze(median(s,2))')';
        
        sfull{j,k} = squeeze(median(s,2));
        ssfull{j,k} = squeeze(std(s,0,2));
                
        msprobe(:,j,k) = mean(squeeze(median(s_probe,2)),2);
        ssprobe(:,j,k) = sem(squeeze(median(s_probe,2))')';
        
        sprobefull{j,k} = squeeze(median(s_probe,2));
        ssprobefull{j,k} = squeeze(std(s_probe,0,2));
        
%         subplot(nDays_behav,2,k+(2*(j-1)));hold on;
%         plot(ms(:,j,k),'k','linewidth',2);
%         PlotHVLines(0,'h','k:')
%         ylim([-0.03 0.04]);
%         PlotHVLines(pretone*round(acq),'v','color',colors{k},'linewidth',1);     
    end
end
%%
figure;
colors_target = [0*ones(nDays_behav,1) linspace(0.2,0.8,nDays_behav)' 0*ones(nDays_behav,1)];
colors_foil = [linspace(0.2,0.8,nDays_behav)' 0*ones(nDays_behav,1) 0*ones(nDays_behav,1)];   
xmax = nframes_psth-25;
tonenames = {'Target','Foil'};
for j=1:2
    subplot(2,2,j+(j-1));hold on;
    for i=1:nDays_behav
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
    for i=1:nDays_behav
        tms = ms(1:xmax,i,j)-ms(1,i,j);
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
    for i=1:nDays_behav
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
    for i=1:nDays_behav
        tms = msprobe(1:xmax,i,j)-msprobe(1,i,j);
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
%% Look at indiv cell variation

% Amplitude and latency
pretone = 1; % for PSTH tone, in s
posttone = 4; % for PSTH tone, in s
nframes_psth = pretone*round(acq) + posttone*round(acq);
xrange = [-pretone posttone];
savefig = true;
binpeak = -0.1:0.001:0.1;
binlatency = 1:1:15;

if strcmp(mouse,'cd036')
    early = find(ismember(findgooddays,1));
    acquisition = find(ismember(findgooddays,[2;3;4;5;6;7;8;9;10]));
    expression = find(ismember(findgooddays,11));
    expert = find(ismember(findgooddays,[12;13;14;15;16]));
end
if strcmp(mouse,'cd017')
%     early = findgooddays(ismember(findgooddays,[1;2;3]));
%     acquisition = findgooddays(ismember(findgooddays,[4;5]));
%     expression = findgooddays(ismember(findgooddays,[6;7;8;9;10;11]));
%     expert = findgooddays(ismember(findgooddays,[12;13;14]));
    % if findgooddays
    early = find(ismember(findgooddays,[1;2;3]));
    acquisition = find(ismember(findgooddays,[4;5]));
    expression = find(ismember(findgooddays,[6;7;8;9;10;11]));
    expert = find(ismember(findgooddays,[12;13;14]));
end
groups = {early,acquisition,expression,expert};

f=figure;
for i=1:length(groups)
    group = groups{i};
    tt = [];tf = [];
    for j=1:length(group)
        g = group(j);
        tt(:,:,j) = sfull{g,1};
        disp(['groups = ' num2str(i) ', nCells = ' num2str(size(sfull{g,1},2))])
        tf(:,:,j) = sfull{g,2};
    end
    tt = nanmean(tt,3);
    tf = nanmean(tf,3);
    subplot(4,2,i+(i-1));hold on;
    [~,ordertarget] = sort(mean(tt(15:25,:)));
    PlotColorCurves(tt(:,ordertarget)',xrange,'cutoffs',[-0.05 0.05]);
    PlotHVLines(0,'w');
    
    subplot(4,2,i+(i-1)+1);hold on;
    [~,orderfoil] = sort(mean(tf(15:25,:)));
    PlotColorCurves(tf(:,orderfoil)',xrange,'cutoffs',[-0.05 0.05]);
    PlotHVLines(0,'w');
end
if savefig
    if strcmp(mouse,'cd036')
        saveas(f,['T:\LabData5\cd036\analysis\' mouse '_PsthTF_Reinforced_LearningStages.pdf']);
    elseif strcmp(mouse,'cd017')
        saveas(f,['U:\LabData4\celine\cd017\analysis\' mouse '_PsthTF_Reinforced_LearningStages.pdf']);
    end
    close(f);
end
smoothfactor = [0 3];
f=figure;
for i=1:length(groups)
    group = groups{i};
    tt = [];tf = [];
    for j=1:length(group)
        g = group(j);
        tt(:,:,j) = sprobefull{g,1};
        disp(['groups = ' num2str(i) ', nCells = ' num2str(size(sprobefull{g,1},2))])
        tf(:,:,j) = sprobefull{g,2};
    end
    tt = Smooth(nanmean(tt,3),smoothfactor);
    tf = Smooth(nanmean(tf,3),smoothfactor);
    subplot(4,2,i+(i-1));hold on;
    [~,ordertarget] = sort(nanmean(tt(15:25,:)));
    PlotColorCurves(tt(:,ordertarget)',xrange,'cutoffs',[-0.05 0.05]);
    PlotHVLines(0,'w');
    
    subplot(4,2,i+(i-1)+1);hold on;
    [~,orderfoil] = sort(mean(tf(15:25,:)));
    PlotColorCurves(tf(:,orderfoil)',xrange,'cutoffs',[-0.05 0.05]);
    PlotHVLines(0,'w');
end
if savefig
    if strcmp(mouse,'cd036')
        saveas(f,['T:\LabData5\cd036\analysis\' mouse '_PsthTF_Probe_LearningStages.pdf']);
    elseif strcmp(mouse,'cd017')
        saveas(f,['U:\LabData4\celine\cd017\analysis\' mouse '_PsthTF_Probe_LearningStages.pdf']);
    end
    close(f);
end

f=figure;
for i=1:length(groups)
    group = groups{i};
    tt = [];tf = [];
    for j=1:length(group)
        g = group(j);
        tt(:,:,j) = sfull{g,1};
        disp(['groups = ' num2str(i) ', nCells = ' num2str(size(sfull{g,1},2))])
        tf(:,:,j) = sfull{g,2};
    end
    tt = nanmean(tt,3);
    tf = nanmean(tf,3);
    subplot(4,2,i+(i-1));hold on;
    dtf = tt-tf;
    [~,orderdtf] = sort(mean(dtf(15:25,:)));
    PlotColorCurves(dtf(:,orderdtf)',xrange,'cutoffs',[-0.05 0.05]);
    PlotHVLines(0,'w');
    
    subplot(4,2,i+(i-1)+1);hold on;
    tpref = mean(dtf(15:25,:))>0;
    fpref = mean(dtf(15:25,:))<0;
    plot(mean(dtf(:,tpref),2),'g');plot(mean(dtf(:,fpref),2),'r');
    ylim([-0.05 0.05]);
    PlotHVLines(0,'h','k:');
    PlotHVLines(pretone*round(acq),'v','k');    
end
if savefig
    if strcmp(mouse,'cd036')
        saveas(f,['T:\LabData5\cd036\analysis\' mouse '_DiffTF_Reinforced_LearningStages.pdf']);
    elseif strcmp(mouse,'cd017')
        saveas(f,['U:\LabData4\celine\cd017\analysis\' mouse '_DiffTF_Reinforced_LearningStages.pdf']);
    end
    close(f);
end
smoothfactor = [0 3];
f=figure;
for i=1:length(groups)
    group = groups{i};
    tt = [];tf = [];
    for j=1:length(group)
        g = group(j);
        tt(:,:,j) = sprobefull{g,1};
        disp(['groups = ' num2str(i) ', nCells = ' num2str(size(sprobefull{g,1},2))])
        tf(:,:,j) = sprobefull{g,2};
    end
    tt = Smooth(nanmean(tt,3),smoothfactor);
    tf = Smooth(nanmean(tf,3),smoothfactor);
    subplot(4,2,i+(i-1));hold on;
    dtf = tt-tf;
    [~,orderdtf] = sort(nanmean(dtf(15:25,:)));
    PlotColorCurves(dtf(:,orderdtf)',xrange,'cutoffs',[-0.05 0.05]);
    PlotHVLines(0,'w');
    
    subplot(4,2,i+(i-1)+1);hold on;
    tpref = mean(dtf(15:25,:))>0;
    fpref = mean(dtf(15:25,:))<0;
    plot(nanmean(dtf(:,tpref),2),'g');plot(nanmean(dtf(:,fpref),2),'r');
    ylim([-0.05 0.05]);
    PlotHVLines(0,'h','k:');
    PlotHVLines(pretone*round(acq),'v','k');    
end
if savefig
    if strcmp(mouse,'cd036')
        saveas(f,['T:\LabData5\cd036\analysis\' mouse '_DiffTF_Probe_LearningStages.pdf']);
    elseif strcmp(mouse,'cd017')
        saveas(f,['U:\LabData4\celine\cd017\analysis\' mouse '_DiffTF_Probe_LearningStages.pdf']);
    end
    close(f);
end

%% 
namestones = {'Target','Foil'};
namesctxs = {'Reinf','Probe'};
colors_target = [0*ones(n,1) linspace(0.2,0.8,n)' 0*ones(n,1)];
colors_foil = [linspace(0.2,0.8,n)' 0*ones(n,1) 0*ones(n,1)];   
xmax = nframes_psth;
xrange = -pretone*round(acq):posttone*round(acq)-1;
toneframe = 0;
savefig = true;

% f=figure;
% for j=1:2
%     subplot(1,4,j);hold on;
%     for i=1:n
%         tms = ms(1:xmax,i,j)-mean(ms(1:pretone*round(acq),i,j));
%         if ismember(j,1)
%             shadedErrorBar(xrange,tms,ss(1:xmax,i,j),{'color',colors_target(i,:)});
%         else
%             shadedErrorBar(xrange,tms,ss(1:xmax,i,j),{'color',colors_foil(i,:)});
%         end
%     end
%     ylim([-0.07 0.07]);
%     PlotHVLines(toneframe,'v','k','linewidth',1);   
%     PlotHVLines(0,'h','k:','linewidth',1);
%     xlabel('Frames');ylabel('df/f');
%     title([namesctxs{1} ' - ' namestones{j}]);
% end
% for j=1:2
%     subplot(1,4,2+j);hold on;
%     for i=1:n
%         tms = msprobe(1:xmax,i,j)-mean(msprobe(1:pretone*round(acq),i,j));
%         if ismember(j,1)
%             shadedErrorBar(xrange,tms,ssprobe(1:xmax,i,j),{'color',colors_target(i,:)});
%         else
%             shadedErrorBar(xrange,tms,ssprobe(1:xmax,i,j),{'color',colors_foil(i,:)});
%         end
%     end
%     ylim([-0.07 0.07]);
%     xlabel('Frames');ylabel('df/f');
%     PlotHVLines(toneframe,'v','k','linewidth',1);   
%     PlotHVLines(0,'h','k:','linewidth',1);
%     title([namesctxs{2} ' - ' namestones{j}]);
% end

% f=figure;
% for j=1:2
%     subplot(1,4,j);hold on;
%     for i=1:n
%         tms = ms(1:xmax,i,j)-mean(ms(1:pretone*round(acq),i,j));
%         if ismember(j,1)
%             plot(xrange,tms,'color',colors_target(i,:));
%         else
%             plot(xrange,tms,'color',colors_foil(i,:));
%         end
%     end
%     ylim([-0.07 0.07]);
%     PlotHVLines(toneframe,'v','k','linewidth',1);   
%     PlotHVLines(0,'h','k:','linewidth',1);
%     xlabel('Frames');ylabel('df/f');
%     title([namesctxs{1} ' - ' namestones{j}]);
% end
% for j=1:2
%     subplot(1,4,2+j);hold on;
%     for i=1:n
%         tms = msprobe(1:xmax,i,j)-mean(msprobe(1:pretone*round(acq),i,j));
%         if ismember(j,1)
%             plot(xrange,tms,'color',colors_target(i,:));
%         else
%             plot(xrange,tms,'color',colors_foil(i,:));
%         end
%     end
%     ylim([-0.07 0.07]);
%     xlabel('Frames');ylabel('df/f');
%     PlotHVLines(toneframe,'v','k','linewidth',1);   
%     PlotHVLines(0,'h','k:','linewidth',1);
%     title([namesctxs{2} ' - ' namestones{j}]);
% end

f=figure;
smoothfactor = 0;
colors_target = [0*ones(length(groups),1) linspace(0.2,0.8,length(groups))' 0*ones(length(groups),1)];
colors_foil = [linspace(0.2,0.8,length(groups))' 0*ones(length(groups),1) 0*ones(length(groups),1)];   
for i=1:length(groups)
    group = groups{i};
    tt = [];tf = [];
    for j=1:length(group)
        g = group(j);
        tt(:,j) =  ms(1:xmax,g,1);
        disp(['groups = ' num2str(i) ', nCells = ' num2str(size(ms(1:xmax,g,1),2))])
        tf(:,j) = ms(1:xmax,g,2);
    end
    tt = Smooth(nanmedian(tt,2),smoothfactor);
    tf = Smooth(nanmedian(tf,2),smoothfactor);
    
    subplot(1,4,1);hold on;
    tms = tt(1:xmax,:)-mean(tt(1:pretone*round(acq),:));
    plot(xrange,tms,'color',colors_target(i,:));
    ylim([-0.02 0.06]);
    xlabel('Frames');ylabel('df/f');
    PlotHVLines(toneframe,'v','k','linewidth',1);   
    PlotHVLines(0,'h','k:','linewidth',1);
    
    subplot(1,4,2);hold on;
    tms = tf(1:xmax,:)-mean(tf(1:pretone*round(acq),:));
    plot(xrange,tms,'color',colors_foil(i,:));    
    ylim([-0.02 0.06]);
    xlabel('Frames');ylabel('df/f');
    PlotHVLines(toneframe,'v','k','linewidth',1);   
    PlotHVLines(0,'h','k:','linewidth',1);
end
smoothfactor = 0; 
for i=1:length(groups)
    group = groups{i};
    tt = [];tf = [];
    for j=1:length(group)
        g = group(j);
        tt(:,j) =  msprobe(1:xmax,g,1);
        disp(['groups = ' num2str(i) ', nCells = ' num2str(size(msprobe(1:xmax,g,1),2))])
        tf(:,j) = msprobe(1:xmax,g,2);
    end
    tt = Smooth(nanmedian(tt,2),smoothfactor);
    tf = Smooth(nanmedian(tf,2),smoothfactor);
    
    subplot(1,4,3);hold on;
    tms = tt(1:xmax,:)-mean(tt(1:pretone*round(acq),:));
    plot(xrange,tms,'color',colors_target(i,:));    
    ylim([-0.02 0.06]);
    xlabel('Frames');ylabel('df/f');
    PlotHVLines(toneframe,'v','k','linewidth',1);   
    PlotHVLines(0,'h','k:','linewidth',1);
    
    subplot(1,4,4);hold on;
    tms = tf(1:xmax,:)-mean(tf(1:pretone*round(acq),:));
    plot(xrange,tms,'color',colors_foil(i,:));    
    ylim([-0.02 0.06]);
    xlabel('Frames');ylabel('df/f');
    PlotHVLines(toneframe,'v','k','linewidth',1);   
    PlotHVLines(0,'h','k:','linewidth',1);
end
if savefig
    if strcmp(mouse,'cd036')
        saveas(f,['T:\LabData5\cd036\analysis\' mouse '_MeanAmplitudeTF_ReinfProbe_LearningStages.pdf']);
    elseif strcmp(mouse,'cd017')
        saveas(f,['U:\LabData4\celine\cd017\analysis\' mouse '_MeanAmplitudeTF_ReinfProbe_LearningStages.pdf']);
    end
    close(f);
end


f=figure;
colors_target = [0*ones(length(groups),1) linspace(0.2,0.8,length(groups))' 0*ones(length(groups),1)];
colors_foil = [linspace(0.2,0.8,length(groups))' 0*ones(length(groups),1) 0*ones(length(groups),1)];  
peakwin = 15:27;
for i=1:length(groups)
    group = groups{i};
    tt = [];tf = [];
    for j=1:length(group)
        g = group(j);
        tt(:,:,j) =  sfull{g,1};
        disp(['groups = ' num2str(i) ', nCells = ' num2str(size(sfull{g,1},2))])
        tf(:,:,j) = sfull{g,2};
    end
    tt = nanmean(tt,3);
    tf = nanmean(tf,3);
    
    [~,peaktt] = max(tt(peakwin,:));
    htt= hist(peaktt,0:15);
    
    [~,peaktf] = max(tf(peakwin,:));
    htf= hist(peaktf,0:15);
    
    subplot(1,4,1);hold on;
    plot(0:15,htt,'color',colors_target(i,:));
    PlotHVLines(median(peaktt),'v','color',colors_target(i,:));
   
    subplot(1,4,2);hold on;
    plot(0:15,htf,'color',colors_foil(i,:));    
    PlotHVLines(median(peaktf),'v','color',colors_foil(i,:));
end
for i=1:length(groups)
    group = groups{i};
    tt = [];tf = [];
    for j=1:length(group)
        g = group(j);
        tt(:,:,j) =  sprobefull{g,1};
        disp(['groups = ' num2str(i) ', nCells = ' num2str(size(sprobefull{g,1},2))])
        tf(:,:,j) = sprobefull{g,2};
    end
    tt = nanmean(tt,3);
    tf = nanmean(tf,3);
    
    [~,peaktt] = max(tt(peakwin,:));
    htt= hist(peaktt,0:15);
    
    [~,peaktf] = max(tf(peakwin,:));
    htf= hist(peaktf,0:15);
    
    subplot(1,4,3);hold on;
    plot(0:15,htt,'color',colors_target(i,:));
    PlotHVLines(median(peaktt),'v','color',colors_target(i,:));
   
    subplot(1,4,4);hold on;
    plot(0:15,htf,'color',colors_foil(i,:));    
    PlotHVLines(median(peaktf),'v','color',colors_foil(i,:));
end
if savefig
    if strcmp(mouse,'cd036')
        saveas(f,['T:\LabData5\cd036\analysis\' mouse '_PeakLatencyTF_ReinfProbe_LearningStages.pdf']);
    elseif strcmp(mouse,'cd017')
        saveas(f,['U:\LabData4\celine\cd017\analysis\' mouse '_PeakLatencyTF_ReinfProbe_LearningStages.pdf']);
    end
    close(f);
end

%% plot peak latency across days
colors_actions = {'g','c','y','r'};
figure;hold on;
posttoneframes = round(acq)*pretone:round(acq)*pretone+9;
subplot(2,1,1);hold on;
for j=1:2
    [~,maxframes] = max(ms(posttoneframes,:,j));
    plot(maxframes,'.-','markersize',20,'color',colors_actions{j});
end
ylim([0 10]);
PlotHVLines(0,'h','k:','linewidth',1);
xlabel('Peak latency (frame)');ylabel('df/f peak in tone-evocked window');
title(namesctxs{1});

subplot(2,1,2);hold on;
for j=1:2
    [~,maxframes] = max(ms(posttoneframes,:,j));
    plot(maxframes,'.-','markersize',20,'color',colors_actions{j});
end
ylim([0 10]);
PlotHVLines(0,'h','k:','linewidth',1);
xlabel('Peak latency (frame)');ylabel('df/f peak in tone-evocked window');
title(namesctxs{2});

%%

figure;
nd = 6;
for i=1:nd
%     subplot(nDays_behav,2,i+(i-1));hold on;
    subplot(nd,2,i+(i-1));hold on;
    targettraces = sfull{i,1};
    
%     [peak,latency] = max(abs(targettraces(15:25,:)));
%     hist(peak,binpeak);
%     PlotHVLines(median(peak),'v','r');
%     xlim([-0.1 0.1])
%     
%     subplot(nDays_behav,2,i+(i-1)+1);hold on;
%     hist(latency,binlatency);
%     PlotHVLines(median(latency),'v','r');
    
    
%     figure;
    [~,ordertarget] = sort(mean(targettraces(15:25,:)));
    PlotColorCurves(targettraces(:,ordertarget)',xrange,'cutoffs',[-0.05 0.05]);
    
%     ztargettraces = zscore(targettraces);
%     figure;
%     PlotColorCurves(ztargettraces',xrange,'cutoffs',[-2 3])
    
    foiltraces = sfull{i,2};
% 	figure;
%     subplot(nDays_behav,2,i+(i-1)+1);hold on;
    subplot(nd,2,i+(i-1)+1);hold on;
%     [~,orderfoil] = sort(mean(foiltraces(15:25,:)));
    PlotColorCurves(foiltraces(:,ordertarget)',xrange,'cutoffs',[-0.05 0.05]);    
end


figure;
for i=1:nDays_behav
%     subplot(nDays_behav,2,i+(i-1));hold on;
%     foiltraces = sfull{i,2};
%     [peak,latency] = max(abs(foiltraces(15:25,:)));
%     hist(peak,binpeak);
%     PlotHVLines(median(peak),'v','r');
%     xlim([-0.1 0.1])
%     
%     subplot(nDays_behav,2,i+(i-1)+1);hold on;
%     hist(latency,binlatency);
%     PlotHVLines(median(latency),'v','r');     
    subplot(5,2,i+(i-1));hold on;
    targettraces = sfull{i,2};    
    [~,ordertarget] = sort(mean(targettraces(15:25,:)));
    PlotColorCurves(targettraces(:,ordertarget)',xrange,'cutoffs',[-0.05 0.05]);
    
%     ztargettraces = zscore(targettraces);
%     figure;
%     PlotColorCurves(ztargettraces',xrange,'cutoffs',[-2 3])
    
    foiltraces = sfull{i,2};
% 	figure;
%     subplot(nDays_behav,2,i+(i-1)+1);hold on;
    subplot(5,2,i+(i-1)+1);hold on;
%     [~,orderfoil] = sort(mean(foiltraces(15:25,:)));
    PlotColorCurves(foiltraces(:,ordertarget)',xrange,'cutoffs',[-0.05 0.05]);    


end



%% MEDIAN TRACE (ALL PLANES) ACTIONS REINFORCED AND PROBE PER DAY
pretone = 1; % for PSTH tone, in s
posttone = 4; % for PSTH tone, in s
nframes_psth = pretone*round(acq) + posttone*round(acq);

msaction = nan(nframes_psth,nDays_behav,2);
ssaction = nan(nframes_psth,nDays_behav,2);
msactionprobe = nan(nframes_psth,nDays_behav,2);
ssactionprobe = nan(nframes_psth,nDays_behav,2);

for j=1:nDays_behav            
    this_day = matrix(:,DAY)==dayss_behav(j)+1;    
    for k=1:4        
        for p=1:nPlanes % loop through planes
            signals_thisplane = signals{p};
            trace = signals_thisplane(:,:,2);

            if p==2, prevncells = nCellsHere; end % continue matrix

            here = ishere{p}(:,dayss_behav(j)+1);
            here(isnan(here) | here~=1) = 0;           
            nCellsHere = sum(here);

            ok = this_day & ctx(:,1) & ~nop & ks(:,k); % all tone (i.e. all trials)
            l = sum(ok); % nb of (these) tones
            matidx = repelem(0:nframes_psth-1,l,1);
            idx = (matidx+repelem(matrix(ok,TONEF(p))- pretone*round(acq),1,nframes_psth))';

            ok_probe = this_day & ctx(:,2) & ~nop & ks(:,k); % all tone (i.e. all trials)
            l_probe = sum(ok_probe); % nb of (these) tones
            matidx_probe = repelem(0:nframes_psth-1,l_probe,1);
            idx_probe = (matidx_probe+repelem(matrix(ok_probe,TONEF(p))- pretone*round(acq),1,nframes_psth))';

            if p==1
                s = trace(logical(here),idx(:))';    
                s = reshape(s,nframes_psth,l,nCellsHere);     

                s_probe = trace(logical(here),idx_probe(:))';    
                s_probe = reshape(s_probe,nframes_psth,l_probe,nCellsHere);     
            else
                s2 = trace(logical(here),idx(:))'; 
                s2 = reshape(s2,nframes_psth,l,nCellsHere); 
                s(:,:,prevncells+1:prevncells+nCellsHere) = s2;

                s2_probe = trace(logical(here),idx_probe(:))'; 
                s2_probe = reshape(s2_probe,nframes_psth,l_probe,nCellsHere); 
                s_probe(:,:,prevncells+1:prevncells+nCellsHere) = s2_probe;            
            end
        end

        msaction(:,j,k) = mean(squeeze(median(s,2)),2);
        ssaction(:,j,k) = sem(squeeze(median(s,2))')';
        
        msactionprobe(:,j,k) = mean(squeeze(median(s_probe,2)),2);
        ssactionprobe(:,j,k) = sem(squeeze(median(s_probe,2))')';      
    end
end

%%
figure;
namesactions = {'H','M','FA','CR'};
namesctxs = {'Reinf','Probe'};
colors_target = [0*ones(nDays_behav,1) linspace(0.2,0.8,nDays_behav)' 0*ones(nDays_behav,1)];
colors_foil = [linspace(0.2,0.8,nDays_behav)' 0*ones(nDays_behav,1) 0*ones(nDays_behav,1)];   
xmax = nframes_psth;
xrange = -pretone*round(acq):posttone*round(acq)-1;
toneframe = 0;
for j=1:4
    subplot(1,4,j);hold on;
    for i=1:nDays_behav
        if ismember(j,[1;2])
            shadedErrorBar(xrange,msaction(1:xmax,i,j),ssaction(1:xmax,i,j),{'color',colors_target(i,:)});
        else
            shadedErrorBar(xrange,msaction(1:xmax,i,j),ssaction(1:xmax,i,j),{'color',colors_foil(i,:)});
        end
    end
    ylim([-0.07 0.07]);
    PlotHVLines(toneframe,'v','k','linewidth',1);   
    PlotHVLines(0,'h','k:','linewidth',1);
    xlabel('Frames');ylabel('df/f');
    title([namesctxs{1} ' - ' namesactions{j}]);
end

figure;
for j=1:4
    subplot(1,4,j);hold on;
    for i=1:nDays_behav
        tms = msaction(1:xmax,i,j)-msaction(1,i,j);
        if ismember(j,[1;2])
            shadedErrorBar(xrange,tms,ssaction(1:xmax,i,j),{'color',colors_target(i,:)});
        else
            shadedErrorBar(xrange,tms,ssaction(1:xmax,i,j),{'color',colors_foil(i,:)});
        end
    end
    ylim([-0.07 0.07]);
    PlotHVLines(toneframe,'v','k','linewidth',1);   
    PlotHVLines(0,'h','k:','linewidth',1);
    xlabel('Frames');ylabel('df/f');
    title([namesctxs{1} ' - ' namesactions{j}]);
end


figure;
for j=1:4
    subplot(1,4,j);hold on;
    for i=1:nDays_behav
        if ismember(j,[1;2])
            shadedErrorBar(xrange,msactionprobe(1:xmax,i,j),ssactionprobe(1:xmax,i,j),{'color',colors_target(i,:)});
        else
            shadedErrorBar(xrange,msactionprobe(1:xmax,i,j),ssactionprobe(1:xmax,i,j),{'color',colors_foil(i,:)});
        end
    end
    ylim([-0.07 0.07]);
    xlabel('Frames');ylabel('df/f');
    PlotHVLines(toneframe,'v','k','linewidth',1);   
    PlotHVLines(0,'h','k:','linewidth',1);
    title([namesctxs{2} ' - ' namesactions{j}]);
end

figure;
for j=1:4
    subplot(1,4,j);hold on;
    for i=1:nDays_behav
        tms = msactionprobe(1:xmax,i,j)-msactionprobe(1,i,j);
        if ismember(j,[1;2])
            shadedErrorBar(xrange,tms,ssaction(1:xmax,i,j),{'color',colors_target(i,:)});
        else
            shadedErrorBar(xrange,tms,ssaction(1:xmax,i,j),{'color',colors_foil(i,:)});
        end
    end
    ylim([-0.07 0.07]);
    PlotHVLines(toneframe,'v','k','linewidth',1);   
    PlotHVLines(0,'h','k:','linewidth',1);
    xlabel('Frames');ylabel('df/f');
    title([namesctxs{2} ' - ' namesactions{j}]);
end

%% plot peak across days
colors_actions = {'g','c','y','r'};
figure;hold on;
peakframes = round(acq)*pretone+2:round(acq)*pretone+8;
subplot(2,1,1);hold on;
for j=1:4
    plot(max(msaction(peakframes,:,j)),'.-','markersize',20,'color',colors_actions{j});
end
ylim([-0.03 0.08]);
PlotHVLines(0,'h','k:','linewidth',1);
xlabel('Days');ylabel('df/f peak in tone-evocked window');
title(namesctxs{1});

subplot(2,1,2);hold on;
for j=1:4
    plot(max(msactionprobe(peakframes,:,j)),'.-','markersize',20,'color',colors_actions{j});
end
ylim([-0.03 0.08]);
PlotHVLines(0,'h','k:','linewidth',1);
xlabel('Days');ylabel('df/f peak in tone-evocked window');
title(namesctxs{2});

%% plot peak latency across days
colors_actions = {'g','c','y','r'};
figure;hold on;
posttoneframes = round(acq)*pretone:round(acq)*pretone+9;
subplot(2,1,1);hold on;
for j=1:4
    [~,maxframes] = max(msaction(posttoneframes,:,j));
    plot(maxframes,'.-','markersize',20,'color',colors_actions{j});
end
ylim([0 10]);
PlotHVLines(0,'h','k:','linewidth',1);
xlabel('Peak latency (frame)');ylabel('df/f peak in tone-evocked window');
title(namesctxs{1});

subplot(2,1,2);hold on;
for j=1:4
    [~,maxframes] = max(msactionprobe(posttoneframes,:,j));
    plot(maxframes,'.-','markersize',20,'color',colors_actions{j});
end
ylim([0 10]);
PlotHVLines(0,'h','k:','linewidth',1);
xlabel('Peak latency (frame)');ylabel('df/f peak in tone-evocked window');
title(namesctxs{2});
