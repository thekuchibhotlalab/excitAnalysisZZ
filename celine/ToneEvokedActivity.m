function ToneEvokedActivity(mouse)

global info;

if strcmp(mouse,'cd017')
    suite2ppath = 'I:\celine\cd017\suite2p\';
    h5path = 'V:\LabData4\celine\cd017\h5\'; % h5path = 'H:\celine\cd017\h5\';
    sbxpath = 'V:\LabData4\celine\cd017\';
    behavpath = 'V:\LabData4\celine\cd017\behavior\';
end

if strcmp(mouse,'cd036')
    suite2ppath = 'I:\celine\cd036\suite2p\';
    h5path = 'U:\LabData5\cd036'; 
    sbxpath = 'U:\LabData5\cd036\';
    behavpath = 'U:\LabData5\cd036\behavior\';
end

% Default values
acq = 30.98;
pretone = 1; % for PSTH tone, in s
posttone = 3; % for PSTH tone, in s
% Column names
SESSION = 1;
TRIAL = 2;
TONE = 3;
RESPONSE = 4;
H = 1; M = 2; FA = 3; CR = 4;
NOLICK_PERIOD = 5;
RESP_TIME = 6;
DELAY_AFTER_RESP = 7;
TOTAL_TRIAL_DUR_MINUS_RESP_TIME = 8;
LICKF = 9;
REWARDF = 10;
TOTAL_TRIAL_DUR = 11; % toc 
TONEF = 12;
CONTEXT=13;

% Get the number of frames of all the files process in suite2p
cd(h5path);
files = dir('*.h5');
nFiles = length(files);
names = cell(nFiles,1);
for i=1:nFiles, names{i} = files(i).name(1:end-3); end

cd(sbxpath);
nFrames = nan(nFiles,1);nFrames_add = nan(nFiles,1);
sbxread(names{1},1,1);
nPlanes = info.otparam(3);
nFrames_oneplane = nan(nFiles,nPlanes);
for i=1:nFiles
    sbxread(names{i},1,1);
    nFrames(i) = info.max_idx;    
    if i>1
        nFrames_add(i) = nFrames(i)+nFrames_add(i-1); 
    else        
        nFrames_add(i) = nFrames(i); 
    end
    
    if mod(nFrames(i),2)
        if i>1
            nFrames_oneplane(i,:) = [round(nFrames(i)/nPlanes)+nFrames_oneplane(i-1,1) round(nFrames(i)/nPlanes)-1+nFrames_oneplane(i-1,2)];
        else
            nFrames_oneplane(i,:) = [round(nFrames(i)/nPlanes) round(nFrames(i)/nPlanes)-1];
        end
    else
        nFrames_oneplane(i,:) = [nFrames(i)/nPlanes+nFrames_oneplane(i-1,1) nFrames(i)/nPlanes+nFrames_oneplane(i-1,2)];
    end    
end

% nFrames_oneplane = nan(nFiles,nPlanes);
% nFrames_oneplane(logical(mod(nFrames_add,2)),:) = [round(nFrames_add(logical(mod(nFrames_add,2)))/nPlanes) round(nFrames_add(logical(mod(nFrames_add,2)))/nPlanes)-1];
% nFrames_oneplane(~mod(nFrames_add,2),:) = [nFrames_add(~mod(nFrames_add,2))/nPlanes nFrames_add(~mod(nFrames_add,2))/nPlanes];
nFrames_oneplane = [zeros(1,nPlanes);nFrames_oneplane];

% EXCEPTIONS 'NAMES' FOR CD017
if strcmp(mouse,'cd017'), names{end} = 'cd017_018_000'; end

% Load behavioral files and place them in the right order
cd(behavpath);
behav_files = dir('*_*v*.txt');
f=1;
while ~strcmp(behav_files(f).name,[mouse '_1v1.txt'])
    f=f+1;
end
behav_files = circshift(behav_files,-(f-1));
nSessions = size(behav_files,1);
names_behav = cell(nSessions,1);
for i=1:nSessions, names_behav{i} = behav_files(i).name; end

% define day 1 and nb od days of recording, and find target and foil frames

% find which recording session correspond to the first day of behavior
% idea: check file date
cd(sbxpath);
sbxfiles = dir('*.sbx');
day1_behavior = behav_files(1).date(1:end-8);
good1day='';a=1;
while ~strcmp(good1day,day1_behavior)
    good1day = sbxfiles(a).date(1:end-8);
    a=a+1;
end
day1 = str2double(names{a-1}(7:9));
tcdays = day1-1;
% day1 = str2double(names{1}(7:9)); % of behavior
lastday_behavior = behav_files(end).date(1:end-8);
goodlastday='';a=0;
while ~strcmp(goodlastday,lastday_behavior)
    goodlastday = sbxfiles(end-a).date(1:end-8);
    a=a+1;
end
lastday = str2double(names{end-a-1}(7:9));
lastdaytc = str2double(names{end}(7:9));
if lastdaytc>lastday
    tcdays = [tcdays;(lastday+1:lastdaytc)'];
end
%%
% lastday = str2double(names{end}(7:9));
dayss_behav = day1:lastday; 
dayss = sort([tcdays;dayss_behav']);
nDays_behav = length(dayss_behav);
nDays = nDays_behav + length(tcdays);
matrix= [];
list_b = cell(nDays,1);list_f = cell(nDays,1);
previousTrial = 0;
cd(suite2ppath);
for i=1:nDays
    this_day = dayss(i);
    if length(num2str(this_day))==2, toCompare = {[mouse '_0' num2str(this_day)]};
    else, toCompare = {[mouse '_00' num2str(this_day)]}; end
    toCmpr = {};
    for k=1:nFiles, toCmpr = [toCmpr;toCompare]; end       
    ok = cellfun(@strfind,names,toCmpr,'UniformOutput',false); 
    ok2 = nan(nFiles,1);
    for k=1:nFiles, if isempty(ok{k}), continue, end; ok2(k) = ok{k}; end
    n_rec = nansum(ok2);
    list_rec = find(ok2==1);
    disp(['D' num2str(i) ', files ' num2str(list_rec')]);        

    if ~ismember(this_day,dayss_behav), continue, end

    toCompare = {[mouse '_' num2str(this_day) 'v']};toCmpr = {};
    for k=1:nSessions, toCmpr = [toCmpr;toCompare]; end      
    ok = cellfun(@strfind,names_behav,toCmpr,'UniformOutput',false); 
    ok2 = nan(nFiles,1);
    for k=1:nSessions, if isempty(ok{k}), continue, end; ok2(k) = ok{k}; end
    n_behav = nansum(ok2);
    list_behav = find(ok2==1); 

    
    for k=1:n_behav
        d = importdata([behav_files(list_behav(k)).folder '\' behav_files(list_behav(k)).name]);                

        if (i==1 && k==1)
            foil = unique([d(d(:,RESPONSE)==FA,TONE);d(d(:,RESPONSE)==CR,TONE)]);
            target = unique([d(d(:,RESPONSE)==H,TONE);d(d(:,RESPONSE)==M,TONE)]);
        end
        
        trialnb = d(:,TRIAL) + previousTrial;
        nTrials = length(trialnb);
        
        contxt = d(:,CONTEXT);
        resp = d(:,RESPONSE);
        toneFrame = d(:,TONEF);
        rewardFrame = d(:,LICKF);
        
        toneFramePlane = [toneFrame toneFrame];
        toneFramePlane(logical(mod(toneFrame,2)),:) = [round(toneFrame(logical(mod(toneFrame,2)))/nPlanes) round(toneFrame(logical(mod(toneFrame,2)))/nPlanes)-1];
        toneFramePlane(~mod(toneFrame,2),:) = [toneFrame(~mod(toneFrame,2))/nPlanes toneFrame(~mod(toneFrame,2))/nPlanes];

        rewardFramePlane = [rewardFrame rewardFrame];
        rewardFramePlane(logical(mod(rewardFrame,2)),:) = [round(rewardFrame(logical(mod(rewardFrame,2)))/nPlanes) round(rewardFrame(logical(mod(rewardFrame,2)))/nPlanes)-1];
        rewardFramePlane(~mod(rewardFrame,2),:) = [rewardFrame(~mod(rewardFrame,2))/nPlanes rewardFrame(~mod(rewardFrame,2))/nPlanes];
        
        rewardFramePlane(rewardFrame==1000000,:) = nan;
        
        if strcmp(mouse,'cd036') && (i==2 && k==4)
%             toneFrame = toneFrame+nFrames_add(list_rec(k));            
            toneFramePlane = toneFramePlane+repelem(nFrames_oneplane(list_rec(k+1),:),nTrials,1);
            rewardFramePlane = rewardFramePlane+repelem(nFrames_oneplane(list_rec(k+1),:),nTrials,1);
        else
            if ~(k==1 && i==1)
                toneFramePlane = toneFramePlane+repelem(nFrames_oneplane(list_rec(k),:),nTrials,1);
                rewardFramePlane = rewardFramePlane+repelem(nFrames_oneplane(list_rec(k),:),nTrials,1);
            end
        end            
%         
%         if (i==1 && k==4)
% %             toneFrame = toneFrame+nFrames_add(list_rec(k));            
%             toneFramePlane = toneFramePlane+repelem(nFrames_oneplane(list_rec(k),:),nTrials,1);
%             rewardFramePlane = rewardFramePlane+repelem(nFrames_oneplane(list_rec(k),:),nTrials,1);
%         else
%             if ~(k==1 && i==1)
%                 toneFramePlane = toneFramePlane+repelem(nFrames_oneplane(list_rec(k),:),nTrials,1);
%                 rewardFramePlane = rewardFramePlane+repelem(nFrames_oneplane(list_rec(k),:),nTrials,1);
%             end
%         end            
        matrix = [matrix;i*ones(nTrials,1) list_behav(k)*ones(nTrials,1) trialnb contxt resp toneFramePlane rewardFramePlane];
        
        list_b{i} = list_behav;
        list_f{i} = list_rec;
        previousTrial = previousTrial+nTrials;
    end    
end


%% Column names of the new matrix

DAY = 1;
BLOC = 2;
TRIAL = 3;
CONTEXT = 4;
RESP = 5;
% TONEF = 6;

TONEF1 = 6; % tone frame plane 1
TONEF2 = 7; % tone frame plane 2
REW1 = 8;
REW2 = 9;

% nPlanes = 2;
% tone_fr = matrix(:,TONEF);
% nTones = size(matrix,1);
% if sum(~mod(tone_fr,2))~=nTones
%     start_tone = nan(nTones,nPlanes);
% %     start_tone = [round(tone_fr/nPlanes) round(tone_fr/nPlanes)];
%     start_tone(logical(mod(tone_fr,2)),:) = [round(tone_fr(logical(mod(tone_fr,2)))/nPlanes) round(tone_fr(logical(mod(tone_fr,2)))/nPlanes)-1];
%     start_tone(~mod(tone_fr,2),:) = [tone_fr(~mod(tone_fr,2))/nPlanes tone_fr(~mod(tone_fr,2))/nPlanes];
%     matrix(:,TONEF1) = start_tone(:,1);
%     matrix(:,TONEF2) = start_tone(:,2);
% end

% matrix(:,TONEF1) = matrix(:,TONEF);
% matrix(:,TONEF2) = matrix(:,TONEF);


%% Load cell activity + tracking matrix
cd(suite2ppath);
% ishere = load('ishere_combined.mat');
% ishere = ishere.ishere;
if strcmp(mouse,'cd017')
    cd([suite2ppath 'plane0'])
    ishere1 = load('ishere_plane0.mat');
    ishere1 = ishere1.ishere;
    cd([suite2ppath 'plane1'])
    ishere2 = load('ishere.mat');
    ishere2 = ishere2.ishere;
    ishere = ishere1;
    ishere{2} = ishere2{2};
elseif strcmp(mouse,'cd036')
    cd([suite2ppath 'plane0'])
    ishere1 = load('ishere_plane0.mat');
    ishere1 = ishere1.ishere;
    cd([suite2ppath 'plane1'])
    ishere2 = load('ishere_plane1.mat');
    ishere2 = ishere2.ishere;
    ishere = ishere1;
    ishere{2} = ishere2{2};
end

nbtrialperdays = Accumulate(matrix(:,DAY));
startfrom = dayss_behav(1)+1;
nbtrialperdays = nbtrialperdays(startfrom:end);

signals = cell(nPlanes,1);
for i=1:nPlanes
    cd([suite2ppath 'plane' num2str(i-1)]);
    tc = load([mouse '_TC_plane' num2str(i-1) '.mat']);
    signals{i} = tc.tempTC; 
    for j=1:nFiles
        submat = signals{i}(:,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i));
        signals{i}(:,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i),2) = ...
            (submat-median(submat,2))./median(submat,2);%./repmat(max(submat,[],2)-min(submat,[],2),1,size(submat,2));
    end 
end

%% Define some variables
% figure;plot(signals{1}(2,:,2))
h = matrix(:,RESP)==H;
m = matrix(:,RESP)==M;
fa = matrix(:,RESP)==FA;
cr = matrix(:,RESP)==CR;
ks = [h m fa cr];
ks2 = [(h|m) (fa|cr)];
ctx = [matrix(:,CONTEXT)==1 matrix(:,CONTEXT)==0];

TONEF = [TONEF1 TONEF2];
REWF = [REW1 REW2];


acq = 30.98;
acq = acq/nPlanes;

if strcmp(mouse,'cd017')
    nop = matrix(:,BLOC)==4; % remove bad trials
else
    nop = ~ismember(matrix(:,BLOC),matrix(:,BLOC)); % take everything
end

%% MEDIAN TRACE PER PLANE
for p=1:nPlanes
    signals_thisplane = signals{p};
    trace = signals_thisplane(:,:,2);
    nCells = size(signals_thisplane,1);     
    l = size(matrix,1);
    matidx = repelem(0:nframes_psth-1,l,1);
    idx = (matidx+repelem(matrix(:,TONEF(p))- pretone*round(acq),1,nframes_psth))';
    s = trace(:,idx(:))';
    s = reshape(s,nframes_psth,l*nCells);
    ms = median(s,2);
    figure;hold on;
    plot(ms,'k','linewidth',2);
    ylim([-0.2 1]);
    PlotHVLines(pretone*round(acq),'v','color','k','linewidth',1);    
end


%% MEDIAN TRACE (ALL PLANES) PER DAY

pretone = 1; % for PSTH tone, in s
posttone = 4; % for PSTH tone, in s
nframes_psth = pretone*round(acq) + posttone*round(acq);

ms = nan(nframes_psth,nDays_behav,2);
ss = nan(nframes_psth,nDays_behav,2);
colors = {'g','r'};
figure;hold on;
%%
for j=10:nDays_behav            
    this_day = matrix(:,DAY)==dayss_behav(j)+1;    
    for k=1:2        
        for p=1:nPlanes % loop through planes
            signals_thisplane = signals{p};
            trace = signals_thisplane(:,:,2);

            if p==2, prevncells = nCellsHere; end % continue matrix

            here = ishere{p}(:,dayss_behav(i)+1);
            here(isnan(here) | here~=1) = 0;           
            nCellsHere = sum(here);

            ok = this_day & ctx(:,1) & ~nop & ks2(:,k); % all tone (i.e. all trials)
            l = sum(ok); % nb of (these) tones
            matidx = repelem(0:nframes_psth-1,l,1);
            idx = (matidx+repelem(matrix(ok,TONEF(p))- pretone*round(acq),1,nframes_psth))';

            ok_probe = this_day & ctx(:,2) & ~nop; % all tone (i.e. all trials)
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
        
        subplot(nDays_behav,2,k+(2*(j-1)));hold on;
        plot(ms(:,j,k),'k','linewidth',2);
        PlotHVLines(0,'h','k:')
        ylim([-0.03 0.04]);
        PlotHVLines(pretone*round(acq),'v','color',colors{k},'linewidth',1);     
        
    end
end
%%
figure;
% colors_days = repelem(linspace(0,0.8,nDays_behav)',1,3);
colors_target = [0*ones(nDays_behav,1) linspace(0.2,0.8,nDays_behav)' 0*ones(nDays_behav,1)];
colors_foil = [linspace(0.2,0.8,nDays_behav)' 0*ones(nDays_behav,1) 0*ones(nDays_behav,1)];   
xmax = nframes_psth-25;
for j=1:2
    subplot(2,2,j+(j-1));hold on;
    for i=1:nDays_behav
        if j==1
%             plot(ms(:,i,j),'color',colors_target(i,:));
            shadedErrorBar(1:xmax,ms(1:xmax,i,j),ss(1:xmax,i,j),{'color',colors_target(i,:)});
%             plot(ms(:,i,j)+ss(:,i,j),'--','color',colors_target(i,:));
%             plot(ms(:,i,j)-ss(:,i,j),'--','color',colors_target(i,:));
        else
            shadedErrorBar(1:xmax,ms(1:xmax,i,j),ss(1:xmax,i,j),{'color',colors_foil(i,:)});
%             plot(ms(:,i,j),'color',colors_foil(i,:));
%             plot(ms(:,i,j)+ss(:,i,j),'--','color',colors_foil(i,:));
%             plot(ms(:,i,j)-ss(:,i,j),'--','color',colors_foil(i,:));
        end
    end
    ylim([-0.04 0.07]);
    PlotHVLines(pretone*round(acq),'v','k','linewidth',1);   
    PlotHVLines(0,'h','k:','linewidth',1);
    
    subplot(2,2,j+(j-1)+1);hold on;
    for i=1:nDays_behav
        tms = ms(1:xmax,i,j)-ms(1,i,j);
        if j==1
            shadedErrorBar(1:xmax,tms,ss(1:xmax,i,j),{'color',colors_target(i,:)});
        else
            shadedErrorBar(1:xmax,tms,ss(1:xmax,i,j),{'color',colors_foil(i,:)});
        end
    end
    ylim([-0.04 0.07]);
    PlotHVLines(pretone*round(acq),'v','k','linewidth',1);   
    PlotHVLines(0,'h','k:','linewidth',1);
end
  

%% ALL NEURONS, ALL DAYS, ALL CONTEXTS, FOIL/TARGET
% for p=1:nPlanes
%     signals_thisplane = signals{p};
%     nCells = size(signals_thisplane,1);       
%     for i=1:nCells
%         fig=figure;hold on;
%         trace = signals_thisplane(i,:,2);
%         for j=1:nDays_behav            
%             this_day = matrix(:,DAY)==dayss_behav(j)+1;
%             if sum(this_day)==0, continue, end
%             for c = 1:2
%                 cntxt = matrix(:,CONTEXT)==(c-1);
%                 for k=1:2
%                     ok = this_day & cntxt & ks2(:,k);                    
%                     l = sum(ok);
%                     matidx = repelem(0:nframes_psth-1,l,1);
%                     idx = (matidx+repelem(matrix(ok,TONEF(p))- pretone*round(acq),1,nframes_psth))';
%                     s = trace(:,idx(:));
%                     s = reshape(s,nframes_psth,l);
%                     ms = median(s,2);
%                     subplot(nDays_behav,4,k+(2*(c-1))+(4*(j-1)));hold on;
%                     plot(s,'color',[0.8 0.8 0.8]);hold on;
%                     plot(ms,'k','linewidth',2);
%                     ylim([-0.2 1]);
%                     PlotHVLines(pretone*round(acq),'v','color','k','linewidth',1);                                        
%                 end
%             end
%         end
%         pause();
%         close(fig);
%     end
% end

%% ALL NEURONS, ALL DAYS, ALL CONTEXTS, ALL RESPONSES
% for p=1:nPlanes
%     signals_thisplane = signals{p};
%     nCells = size(signals_thisplane,1);       
%     for i=1:nCells
%         fig=figure;hold on;
%         trace = signals_thisplane(i,:,2);
%         for j=1:nDays_behav            
%             this_day = matrix(:,DAY)==dayss_behav(j)+1;
%             if sum(this_day)==0, continue, end
%             for c = 1:2
%                 cntxt = matrix(:,CONTEXT)==(c-1);
%                 for k=1:4
%                     ok = this_day & cntxt & ks(:,k);                    
%                     l = sum(ok);
%                     matidx = repelem(0:nframes_psth-1,l,1);
%                     idx = (matidx+repelem(matrix(ok,TONEF(p))- pretone*round(acq),1,nframes_psth))';
%                     s = trace(:,idx(:));
%                     s = reshape(s,nframes_psth,l);
%                     ms = mean(s,2);
%                     subplot(nDays_behav,8,k+(4*(c-1))+(8*(j-1)));hold on;
%                     plot(s,'color',[0.8 0.8 0.8]);hold on;
%                     plot(ms,'k','linewidth',2);
%                     ylim([-0.2 1]);
%                     PlotHVLines(pretone*round(acq),'v','color','k','linewidth',1);                                        
%                 end
%             end
%         end
%         pause();
%         close(fig);
%     end
% end
%% MEDIAN TRACE (ALL PLANES) ACTIONS REINFORCED AND PROBE PER DAY
pretone = 1; % for PSTH tone, in s
posttone = 4; % for PSTH tone, in s
nframes_psth = pretone*round(acq) + posttone*round(acq);

msaction = nan(nframes_psth,nDays_behav,2);
ssaction = nan(nframes_psth,nDays_behav,2);
msactionprobe = nan(nframes_psth,nDays_behav,2);
ssactionprobe = nan(nframes_psth,nDays_behav,2);
% colors_actions = {'g','g','r','r'};
figure;hold on;

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
        
%         subplot(nDays_behav,4,k+(4*(j-1)));hold on;
%         plot(msaction(:,j,k),'k','linewidth',2);
%         PlotHVLines(0,'h','k:')
%         ylim([-0.03 0.04]);
%         PlotHVLines(pretone*round(acq),'v','color',colors_actions{k},'linewidth',1);         
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

%% PSTH LICKS, MEDIAN TRACE (ALL PLANES) ACTIONS REINFORCED AND PROBE PER DAY
prerew = 1; % for PSTH tone, in s
postrew = 4; % for PSTH tone, in s
nframes_psth = prerew*round(acq) + postrew*round(acq);

msactionRew = nan(nframes_psth,nDays_behav,2);
ssactionRew = nan(nframes_psth,nDays_behav,2);
msactionprobeRew = nan(nframes_psth,nDays_behav,2);
ssactionprobeRew = nan(nframes_psth,nDays_behav,2);
% colors_actions = {'g','g','r','r'};
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
            nCellsHere = sum(here);

            ok = this_day & ctx(:,1) & ~nop & ks2(:,k); % all tone (i.e. all trials)
            l = sum(ok); % nb of (these) tones
            matidx = repelem(0:nframes_psth-1,l,1);
            idx = (matidx+repelem(matrix(ok,REWF(p))- prerew*round(acq),1,nframes_psth))';
            idx = idx(:);
            idx(isnan(idx)) = [];
            l = length(idx)/nframes_psth;
            
            ok_probe = this_day & ctx(:,2) & ~nop & ks2(:,k); % all tone (i.e. all trials)
            l_probe = sum(ok_probe); % nb of (these) tones
            matidx_probe = repelem(0:nframes_psth-1,l_probe,1);
            idx_probe = (matidx_probe+repelem(matrix(ok_probe,REWF(p))- prerew*round(acq),1,nframes_psth))';
            idx_probe = idx_probe(:);
            idx_probe(isnan(idx_probe)) = [];
            l_probe = length(idx_probe)/nframes_psth;
            
            if p==1
                s = trace(logical(here),idx)';    
                s = reshape(s,nframes_psth,l,nCellsHere);     

                s_probe = trace(logical(here),idx_probe)';    
                s_probe = reshape(s_probe,nframes_psth,l_probe,nCellsHere);     
            else
                s2 = trace(logical(here),idx)'; 
                s2 = reshape(s2,nframes_psth,l,nCellsHere); 
                s(:,:,prevncells+1:prevncells+nCellsHere) = s2;

                s2_probe = trace(logical(here),idx_probe)'; 
                s2_probe = reshape(s2_probe,nframes_psth,l_probe,nCellsHere); 
                s_probe(:,:,prevncells+1:prevncells+nCellsHere) = s2_probe;            
            end
        end

        msactionRew(:,j,k) = mean(squeeze(median(s,2)),2);
        ssactionRew(:,j,k) = sem(squeeze(median(s,2))')';
        
        msactionprobeRew(:,j,k) = mean(squeeze(median(s_probe,2)),2);
        ssactionprobeRew(:,j,k) = sem(squeeze(median(s_probe,2))')';
        
%         subplot(nDays_behav,4,k+(4*(j-1)));hold on;
%         plot(msaction(:,j,k),'k','linewidth',2);
%         PlotHVLines(0,'h','k:')
%         ylim([-0.03 0.04]);
%         PlotHVLines(pretone*round(acq),'v','color',colors_actions{k},'linewidth',1);         
    end
end
%%
figure;
namesactions = {'H','FA'};
namesctxs = {'Reinf','Probe'};
colors_target = [0*ones(nDays_behav,1) linspace(0.2,0.8,nDays_behav)' 0*ones(nDays_behav,1)];
colors_foil = [linspace(0.2,0.8,nDays_behav)' 0*ones(nDays_behav,1) 0*ones(nDays_behav,1)];   
xmax = nframes_psth;
xrange = -prerew*round(acq):postrew*round(acq)-1;
rewframe = 0;
for j=1:2
    subplot(1,2,j);hold on;
    for i=1:nDays_behav
        if ismember(j,1)
            shadedErrorBar(xrange,msactionRew(1:xmax,i,j),ssactionRew(1:xmax,i,j),{'color',colors_target(i,:)});
        else
            shadedErrorBar(xrange,msactionRew(1:xmax,i,j),ssactionRew(1:xmax,i,j),{'color',colors_foil(i,:)});
        end
    end
    ylim([-0.07 0.07]);
    PlotHVLines(rewframe,'v','k','linewidth',1);   
    PlotHVLines(0,'h','k:','linewidth',1);
    xlabel('Frames');ylabel('df/f');
    title([namesctxs{1} ' - ' namesactions{j}]);
end

figure;
for j=1:2
    subplot(1,2,j);hold on;
    for i=1:nDays_behav
        if ismember(j,1)
            shadedErrorBar(xrange,msactionprobeRew(1:xmax,i,j),ssactionprobeRew(1:xmax,i,j),{'color',colors_target(i,:)});
        else
            shadedErrorBar(xrange,msactionprobeRew(1:xmax,i,j),ssactionprobeRew(1:xmax,i,j),{'color',colors_foil(i,:)});
        end
    end
    ylim([-0.07 0.07]);
    xlabel('Frames');ylabel('df/f');
    PlotHVLines(rewframe,'v','k','linewidth',1);   
    PlotHVLines(0,'h','k:','linewidth',1);
    title([namesctxs{2} ' - ' namesactions{j}]);
end
%% PSTH LICKS, MEDIAN TRACE (ALL PLANES) ACTIONS REINFORCED AND PROBE PER DAY - DOWNSAMPLING AND COMPARISON
prerew = 1; % for PSTH tone, in s
postrew = 4; % for PSTH tone, in s
nframes_psth = prerew*round(acq) + postrew*round(acq);

signif = nan(nDays_behav,nframes_psth);
S = cell(nDays_behav,2);
SP = cell(nDays_behav,2);
    
for j=1:nDays_behav            
    this_day = matrix(:,DAY)==dayss_behav(j)+1; 
    ntrials = nan(1,2);
    ntrialsprobe = nan(1,2);
    
    for k=1:2       
        for p=1:nPlanes % loop through planes
            signals_thisplane = signals{p};
            trace = signals_thisplane(:,:,2);

            if p==2, prevncells = nCellsHere; end % continue matrix

            here = ishere{p}(:,dayss_behav(j)+1);
            here(isnan(here) | here~=1) = 0;           
            nCellsHere = sum(here);

            ok = this_day & ctx(:,1) & ~nop & ks2(:,k); % all tone (i.e. all trials)
            l = sum(ok); % nb of (these) tones
            matidx = repelem(0:nframes_psth-1,l,1);
            idx = (matidx+repelem(matrix(ok,REWF(p))- prerew*round(acq),1,nframes_psth))';
            idx = idx(:);
            idx(isnan(idx)) = [];
            l = length(idx)/nframes_psth;
            
            ok_probe = this_day & ctx(:,2) & ~nop & ks2(:,k); % all tone (i.e. all trials)
            l_probe = sum(ok_probe); % nb of (these) tones
            matidx_probe = repelem(0:nframes_psth-1,l_probe,1);
            idx_probe = (matidx_probe+repelem(matrix(ok_probe,REWF(p))- prerew*round(acq),1,nframes_psth))';
            idx_probe = idx_probe(:);
            idx_probe(isnan(idx_probe)) = [];
            l_probe = length(idx_probe)/nframes_psth;
            
            if p==1
                s = trace(logical(here),idx)';    
                s = reshape(s,nframes_psth,l,nCellsHere);     

                s_probe = trace(logical(here),idx_probe)';    
                s_probe = reshape(s_probe,nframes_psth,l_probe,nCellsHere);     
            else
                s2 = trace(logical(here),idx)'; 
                s2 = reshape(s2,nframes_psth,l,nCellsHere); 
                s(:,:,prevncells+1:prevncells+nCellsHere) = s2;

                s2_probe = trace(logical(here),idx_probe)'; 
                s2_probe = reshape(s2_probe,nframes_psth,l_probe,nCellsHere); 
                s_probe(:,:,prevncells+1:prevncells+nCellsHere) = s2_probe;            
            end
        end
        ntrials(:,k) = size(s,2);
        ntrialsprobe(:,k) = size(s_probe,2);
        S{j,k} = s;
        SP{j,k} = s_probe;
    end
        
%     % Check nb of trials in H vs FA (separated between taget and foil but take only trials with licks)
%     [minn,wmin] = min(ntrials);
%     
%     % for this day, compare trace with minn nb of trial to bootstrap
%     % distribution from the larger trial set
%     larger = [1;2];larger = larger(~ismember(larger,wmin));
%     
%     figure;plot(squeeze(mean(mean(S{wmin},2),3)));
%     
%     nbootstrap = 500;
%     bootstrap_traces = nan(nbootstrap,nframes_psth);
%     randomselected = randn(ntrials(larger),nbootstrap);
%     [~,randomselected] = sort(randomselected,1);
%     randomselected = randomselected(1:minn,:);
%     for b = 1:nbootstrap        
%         bootstrap_traces(b,:) = squeeze(mean(mean(S{larger}(:,randomselected(:,b),:),2),3));
%     end
%     
% %     figure;hold on;plot(bootstrap_traces','k');
%     
%     mbootstrap = mean(bootstrap_traces);
% %     stdbootstrap = std(bootstrap_traces);
%     ybootstrap = quantile(bootstrap_traces,.99);
%     figure;hold on;
%     shadedErrorBar(1:nframes_psth,mbootstrap,ybootstrap);
%     plot(squeeze(mean(mean(S{wmin},2),3)));
%     
%     signif(j,:) = squeeze(mean(mean(S{wmin},2),3))>ybootstrap';
%     sw = find(signif);
%     PlotHVLines(sw,'v','g');
end

%% Ratio baseline / reward

baseline = prerew*round(acq)-9:prerew*round(acq);
peakrew = prerew*round(acq)+2:prerew*round(acq)+11;
rewframe = prerew*round(acq);
ratiopercell = cell(nDays_behav,1);
namesactions = {'H','FA'};
figk2 = figure;
figk1 = figure;
% figure;
for i=1:nDays_behav
    for k=1:2
        bef = squeeze(mean(SP{i,k}(baseline,:,:))); % trials x neurons
        aft = squeeze(mean(SP{i,k}(peakrew,:,:))); % trials x neurons
        pttestpos = nan(size(bef,2),1);
        pttestneg = nan(size(bef,2),1);
        for c = 1:size(bef,2) % loop through neurons
            [~,pttestneg(c)] = ttest2(bef(:,c),aft(:,c),'tail','right'); % is it significantly responding to lick
            [~,pttestpos(c)] = ttest2(aft(:,c),bef(:,c),'tail','right');
        end
        % plot traces of signif action modulated neurons
        mtracespos = squeeze(mean(SP{i,k}(:,:,pttestpos<=0.05),2)); % psth x signif neurons 
        mtracesneg = squeeze(mean(SP{i,k}(:,:,pttestneg<=0.05),2)); % psth x signif neurons 
        
        figure;
        subplot(2,2,1); hold on;
        plot(mtracespos);
        ylim([-0.1 1]);
        PlotHVLines(rewframe,'v','k','linewidth',1);
        PlotIntervals([peakrew(1) peakrew(end)]);
        title(['Reward - ' namesactions{k} ' - Day ' num2str(i)]);
        subplot(2,2,2); hold on;
        [~,order] = sort(max(mtracespos(peakrew,:)));
        PlotColorCurves(mtracespos(:,order)');
        PlotHVLines(rewframe,'v','w','linewidth',1);
        PlotIntervals([peakrew(1) peakrew(end)]);
        clim([0 0.8]);
        
        subplot(2,2,3); hold on;
        plot(mtracesneg);
        ylim([-0.1 1]);
        PlotHVLines(rewframe,'v','k');
        PlotIntervals([peakrew(1) peakrew(end)]);
        subplot(2,2,4); hold on;
        [~,order] = sort(max(abs(mtracesneg(peakrew,:))));
        PlotColorCurves(mtracesneg(:,order)');
        PlotHVLines(rewframe,'v','w','linewidth',1);
        PlotIntervals([peakrew(1) peakrew(end)]);
        clim([-0.02 0.1]);
        
        ratiopercell{i} = mean(bef)./mean(aft);
        
%         mratio = mean(ratiopercell{i});
%         semratio = sem(ratiopercell{i});
%         if k==1, figure(figk1); else, figure(figk2); end
%         hold on;
%         plot(i,mratio,'.','markersize',20);
%         errorbar(i,mratio,semratio);
    end
    
end

%% Find responsive cells
showfig = false;
toneresponsive = cell(nPlanes,1);
for p=1:nPlanes
    signals_thisplane = signals{p};
    nCells = size(signals_thisplane,1); 
    toneresponsive{p} = nan(nCells,1);
    for i=1:nCells
%     for i=[86;90;91;94;97;100;101]'
        trace = signals_thisplane(i,:,2);
        l = size(matrix,1);
        matidx = repelem(0:nframes_psth-1,l,1);
        idx = (matidx+repelem(matrix(:,TONEF(p))- pretone*round(acq),1,nframes_psth))';
        s = trace(:,idx(:))';
        s = reshape(s,nframes_psth,l);
        ms = median(s,2);
        bef=median(s(5:15,:));
        aft=median(s(17:26,:));
        [h0,pttest] = ttest(bef(:),aft(:),'alpha',0.001);  
        toneresponsive{p}(i) = h0;
                        
        if showfig
            fig=figure;hold on;
            plot(ms,'k','linewidth',2);
            PlotHVLines(pretone*round(acq),'v','color','k','linewidth',1);  
            title(['Cell ' num2str(i) ', ' num2str(length(bef)) ', h0 = ' num2str(h0) ', p=' num2str(pttest)]); 
            PlotIntervals([[5 15];[17 26]]);
            pause();
            close(fig);
        end
    end
end

%% LEARNING GROUPS - Target VS Foil
if strcmp(mouse,'cd017')
    early = 1:3;
    acquisition = 4:5;
    expression = 6:11;
    expert = 12:15;
elseif strcmp(mouse,'cd036')
    early = 1;
    acquisition = 2:10;
    expression = 11;
    expert = 12:16;
end
groups = {early,acquisition,expression,expert};
groupNames = {'Early','Acquisition','Expression','Expert'};
nGroups = length(groups);
colors = {'g','r'};
showfig = true;
toneNames = {'Target','Foil'};

RES = cell(nGroups,2);
for p=1:nPlanes
    signals_thisplane = signals{p};
        
    trace = signals_thisplane(:,:,2);        
    here = ishere{p};
    here(isnan(here) | here~=1) = 0;
        
    responsive = toneresponsive{p};
    
    for j=1:nGroups         
        groupsok = groups{j}+1;
        cellishere = here(:,groupsok);
        cellhere_group = any(cellishere,2);  
        okcells = cellhere_group & responsive;
        
        these_days = ismember(matrix(:,DAY),groupsok);            
        
        for k=1:2
            ok = these_days & ks2(:,k);                    
            l = sum(ok);
            matidx = repelem(0:nframes_psth-1,l,1);
            idx = (matidx+repelem(matrix(ok,TONEF(p))- pretone*round(acq),1,nframes_psth))';
            
            nCells = sum(okcells);
            s = trace(logical(okcells),idx(:))';
            
            s = reshape(s,nframes_psth,l,nCells);
            ms = median(s,2);
            ms = permute(ms,[1 3 2])';
            
%             RES{j+(p-1)*4,k} = ms;
            if p==1
                RES{j,k} = ms;
            else
                RES{j,k} = [RES{j,k};ms];
            end
        end
    end
end 

% Plot
if showfig,
    fig=figure;hold on;
    for i=1:nGroups
        for j=1:2 % tones
            subplot(nGroups,2,j+(2*(i-1)));hold on;
            ms = RES{i,j};
            if i==1 && j==1
                [~,order] = sort(median(ms(:,30:60),2));
            end
            PlotColorCurves(ms(order,:),[-1 3]);
            PlotHVLines(0,'v','color','w','linewidth',1);
            title([groupNames{i} ' - ' toneNames{j}]);
            clim([-0.04 0.15]);
            ylabel('Cells');xlabel('Time (s)');
        end
    end
end

% Plot mean curves
if showfig,
    fig=figure;hold on;
    for i=1:nGroups
        for j=1:2 % tones
            subplot(nGroups,2,j+(2*(i-1)));hold on;
            ms = RES{i,j};
            if i==1 && j==1
                [~,order] = sort(median(ms(:,30:60),2));
            end
            plot(mean(ms),colors{j},'linewidth',1);
            ylim([-0.02 0.045]);
            PlotHVLines(pretone*round(acq),'v','k','linewidth',1);
            PlotHVLines([0;0.02],'h','k:','linewidth',1);
            title([groupNames{i} ' - ' toneNames{j}]);            
            ylabel('Cells');xlabel('Time (s)');
        end
    end
end

% Plot mean curves (on the same plot for each group)
colors_foil = {[0.2 0 0];[0.5 0 0];[0.8 0 0];[1 0 0]};
colors_target = {[0 0.2 0];[0 0.5 0];[0 0.8 0];[0 1 0]};
colors_groups = [colors_target colors_foil];

if showfig,
    fig=figure;hold on;
    for i=1:nGroups
        for j=1:2 % tones
            subplot(1,2,j);hold on;
            ms = RES{i,j};
            plot(mean(ms),'color',colors_groups{i,j},'linewidth',1);
            ylim([-0.02 0.045]);
            PlotHVLines(pretone*round(acq),'v','k','linewidth',1);
            PlotHVLines([0;0.02],'h','k:','linewidth',1);
            title([groupNames{i} ' - ' toneNames{j}]);            
            ylabel('Cells');xlabel('Time (s)');
        end
    end
end

% Same, indiv neurons
if showfig,
    fig=figure;hold on;
    for i=1:nGroups
        for j=1:2 % tones
            subplot(1,2,j);hold on;
            ms = RES{i,j};
            [y_amplitude,x_deltaPeak] = max(ms(:,15:30),[],2);
            plot(x_deltaPeak,y_amplitude,'.','color',colors_groups{i,j},'markersize',10);
            ylim([-0.1 1.2]);
            PlotHVLines(0,'h','k:','linewidth',1);
            PlotHVLines(0,'v','k','linewidth',1);
%             title([groupNames{i} ' - ' toneNames{j}]);            
            ylabel('df/f');xlabel('Delta peak from tone (frame)');
        end
    end
end

% same, deeper - hist delta peak from tone
if showfig,
    fig=figure;hold on;
    for i=1:nGroups
        for j=1:2 % tones
            subplot(1,2,j);hold on;
            ms = RES{i,j};
            [~,x_deltaPeak] = max(ms(:,15:30),[],2);
            [h,ht] = hist(x_deltaPeak,0:15);
            plot(ht,h/sum(h),'color',colors_groups{i,j},'linewidth',2);
            ylim([0 0.5]);
            PlotHVLines(0,'h','k:','linewidth',1);
            PlotHVLines(0,'v','k','linewidth',1);
            ylabel('Probability');xlabel('Delta peak from tone (frame)');
            
        end
    end
end

%% LEARNING GROUPS - Diff Target & Foil

% Plot
if showfig
    fig=figure;hold on;
    for i=1:nGroups
        subplot(nGroups,2,1+(2*(i-1)));hold on;
        msdiff = RES{i,1}-RES{i,2};
        if i==1
%             [~,order] = sort(max(msdiff(:,15:30),[],2));
            [~,order] = sort(median(msdiff(:,15:30),2));
        end
        PlotColorCurves(msdiff(order,:),[-1 3]);
        PlotHVLines(0,'v','color','w','linewidth',1);
        title([groupNames{i} ' - Diff(T-F)']);
        clim([-0.01 0.06]);
        ylabel('Cells');xlabel('Time (s)');
    end
end

% Plot - mean curves
if showfig
    for i=1:nGroups
        subplot(nGroups,2,2+(2*(i-1)));hold on;
        msdiff = RES{i,1}-RES{i,2};
        mmsdiff = mean(msdiff);
        plot(mmsdiff,'k','linewidth',1);
        ylim([-0.04 0.02]);
        PlotHVLines(pretone*round(acq),'v','color','k','linewidth',1);
        PlotHVLines([-0.02;0;0.02],'h','k:','linewidth',1);
        title([groupNames{i} ' - Diff(T-F)']);
        ylabel('df/f');xlabel('Frames');
    end
end


%% LEARNING GROUPS - ALL CONTEXTS - Target VS Foil

RES = cell(nGroups,4);
for p=1:nPlanes
    signals_thisplane = signals{p};
        
    trace = signals_thisplane(:,:,2);        
    here = ishere{p};
    here(isnan(here) | here~=1) = 0;
    responsive = toneresponsive{p};
    
    for j=1:nGroups         
        groupsok = groups{j}+1;
        cellishere = here(:,groupsok);
        cellhere_group = any(cellishere,2);  
%         okcells = cellhere_group;
        okcells = cellhere_group & responsive;
        
        these_days = ismember(matrix(:,DAY),groupsok);            
        for c = 1:2
            cntxt = matrix(:,CONTEXT)==(c-1);
            for k=1:2
                ok = these_days & cntxt & ks2(:,k);                    
                l = sum(ok);
                matidx = repelem(0:nframes_psth-1,l,1);
                idx = (matidx+repelem(matrix(ok,TONEF(p))- pretone*round(acq),1,nframes_psth))';

                nCells = sum(okcells);
                s = trace(logical(okcells),idx(:))';

                s = reshape(s,nframes_psth,l,nCells);
                ms = median(s,2);
                ms = permute(ms,[1 3 2])';
        
    %             RES{j+(p-1)*4,k} = ms;
                if p==1
                    RES{j,k+(c-1)*2} = ms;
                else
                    RES{j,k+(c-1)*2} = [RES{j,k+(c-1)*2};ms];
                end
            end
        end
    end
end 


% Plot
ctxtNames = {'Probe','Reinf'};
showfig = true;
if showfig
    fig=figure;hold on;
    for i=1:nGroups
        for c=1:2
            for j=1:2 % tones
                subplot(nGroups,4,j+(2*(c-1))+(4*(i-1)));hold on;
                ms = RES{i,j+(c-1)*2};
                if i==1 && j==1 && c==1
                    [~,order] = sort(mean(ms(:,17:26),2));
                end
                PlotColorCurves(ms(order,:),[-1 3]);
                PlotHVLines(0,'v','color','w','linewidth',1);
                title([groupNames{i} ' - ' ctxtNames{c} ' - ' toneNames{j}]);
                clim([-0.04 0.15]);
                ylabel('Cells');xlabel('Time (s)');
            end
        end
    end
end


%% LEARNING GROUPS - ALL CONTEXT - Diff Target & Foil

% Plot
if showfig
    fig=figure;hold on;
    for i=1:nGroups
        for c=1:2
            subplot(nGroups,2,c+(2*(i-1)));hold on;
            msdiff = RES{i,1+(c-1)*2}-RES{i,2+(c-1)*2};
            if i==1 && c==1
                [~,order] = sort(mean(msdiff(:,17:26),2));
            end
            PlotColorCurves(msdiff(order,:),[-1 3]);
            PlotHVLines(0,'v','color','w','linewidth',1);
            title([groupNames{i} ' - ' ctxtNames{c} ' - Diff(T-F)']);
            clim([-0.04 0.15]);
            ylabel('Cells');xlabel('Time (s)');
        end
    end
end

cm = colormap;


% Plot - mean curves - Diff(T-F)
if showfig
    fig=figure;hold on;
    for i=1:nGroups
        for c=1:2
            subplot(nGroups,2,c+(2*(i-1)));hold on;
            msdiff = RES{i,1+(c-1)*2}-RES{i,2+(c-1)*2};
            mmsdiff = mean(msdiff);
            plot(mmsdiff,'k','linewidth',1);
            ylim([-0.04 0.05]);
            PlotHVLines(pretone*round(acq),'v','k','linewidth',1);
            
            title([groupNames{i} ' - ' ctxtNames{c} ' - Diff(T-F)']);
            
            ylabel('df/f');xlabel('Frames');
        end
    end
end

%% Plot - hist mean tone evocked resp diff - Diff(T-F)
colors_probe = {[0 0 0.1];[0 0 0.5];[0 0 0.8];[0 0 1]};
colors_reinforced = {[0.2 0.2 0.2];[0.4 0.4 0.4];[0.6 0.6 0.6];[0.8 0.8 0.8]};
colors_contexts = [colors_probe colors_reinforced];

minmeantoneresp = -0.2;
maxmeantoneresp = 0.2;
minmaxtoneresp = -0.2;
maxmaxtoneresp = 0.2;
for i=1:nGroups
    for c=1:2
        msdiff = RES{i,1+(c-1)*2}-RES{i,2+(c-1)*2};
        meantoneresp = mean(msdiff(:,17:26),2);
        minmeantoneresp = min([minmeantoneresp;min(meantoneresp)]);
        maxmeantoneresp = max([maxmeantoneresp;max(meantoneresp)]);
        maxtoneresp = max(msdiff(:,17:26),[],2);
        minmaxtoneresp = min([minmaxtoneresp;min(maxtoneresp)]);
        maxmaxtoneresp = max([maxmaxtoneresp;max(maxtoneresp)]);
    end
end

binsmax = linspace(minmaxtoneresp,maxmaxtoneresp,1000);
bins = linspace(minmeantoneresp,maxmeantoneresp,1000);
if showfig
    fig=figure;hold on;
    for i=1:nGroups
        for c=1:2
%             subplot(nGroups,2,c+(2*(i-1)));hold on;
            subplot(2,2,c);hold on;
            msdiff = RES{i,1+(c-1)*2}-RES{i,2+(c-1)*2};
            
            meantoneresp = mean(msdiff(:,17:26),2);
%             [~,where] = sort(meantoneresp);
%             [~,where] = max(abs(meantoneresp),[],2);            
%             [h,ht] = hist(msdiff(:,15+where),bins);
            [h,ht] = hist(meantoneresp,bins);
            s = CumSum(h/sum(h));
            toplot = ht>=-0.05 & ht<=0.05;
            plot(ht(toplot),s(toplot),'color',colors_contexts{i,c},'linewidth',2);
            PlotHVLines(0,'v','k:');
            title('Mean in tone evoked window');
            ylim([0 1]);
            
%             subplot(2,2,2+c);hold on;
%             maxtoneresp = max(msdiff(:,17:26),[],2);
%             [h,ht] = hist(maxtoneresp,binsmax);
%             s = CumSum(h/sum(h));
%             toplot = ht>=-0.2 & ht<=0.2;
%             plot(ht(toplot),s(toplot),'color',colors_contexts{i,c},'linewidth',2);
%             PlotHVLines(0,'v','k:');
%             title('Max in tone evoked window');
%             ylim([0 1]);
        end
    end
end

%% Reinf VS probe
% Plot - mean curves - Diff(Reinf-Probe)
if showfig
    fig=figure;hold on;
    for i=1:nGroups
        for t=1:2
            subplot(nGroups,2,t+(2*(i-1)));hold on;
            msdiff = RES{i,1+(t-1)*1}-RES{i,3+(t-1)*1};
            mmsdiff = mean(msdiff);
            plot(mmsdiff,colors{t},'linewidth',1);
            ylim([-0.04 0.02]);
            PlotHVLines(pretone*round(acq),'v','color','k','linewidth',1);
            PlotHVLines([-0.02;0],'h','k:','linewidth',1);
            title([groupNames{i} ' - ' toneNames{t} ' - Diff(Probe-Reinf)']);
            ylabel('df/f');xlabel('Frames');
        end
    end
end

if showfig
    fig=figure;hold on;
    for i=1:nGroups
        for t=1:2
            subplot(nGroups,2,t+(2*(i-1)));hold on;
            msdiff = RES{i,1+(t-1)*1}-RES{i,3+(t-1)*1};
            PlotColorCurves(msdiff);
            clim([-0.04 0.15]);
%             ylim([-0.04 0.02]);
%             PlotHVLines(pretone*round(acq),'v','color','k','linewidth',1);
%             PlotHVLines([-0.02;0],'h','k:','linewidth',1);
            title([groupNames{i} ' - ' toneNames{t} ' - Diff(Probe-Reinf)']);
            ylabel('df/f');xlabel('Frames');
        end
    end
end

%% Plot FOV + toneresp index
% LEARNING GROUPS - ALL CONTEXTS - Target VS Foil

% Compute average FoV for each day
avg_day = cell(nDays,2);
cd(suite2ppath);
for j=1:nPlanes
    cd([suite2ppath 'plane' num2str(j-1)])
%     avg_session = load(['MeanImgPerSessions_Plane' num2str(j) '.mat']);
    avg_session = load([mouse '_MeanImgPerSessions_Plane' num2str(j) '.mat']);
    avg_session = avg_session.tosave;
    for i=1:nDays
        this_day = dayss(i);
        if length(num2str(this_day))==2, toCompare = {[mouse '_0' num2str(this_day)]};
        else, toCompare = {[mouse '_00' num2str(this_day)]}; end
        toCmpr = {};
        for k=1:nFiles, toCmpr = [toCmpr;toCompare]; end       
        ok = cellfun(@strfind,names,toCmpr,'UniformOutput',false); 
        ok2 = nan(nFiles,1);
        for k=1:nFiles, if isempty(ok{k}), continue, end; ok2(k) = ok{k}; end
        n_rec = nansum(ok2);
        list_rec = find(ok2==1);

        m = [];
        for k=1:n_rec
            m(:,:,k) = avg_session{list_rec(k)};
        end
        avg_day{i,j} = mean(m,3);
    end
end
%%
fig=figure;
for p=1:nPlanes
    cd([suite2ppath 'plane' num2str(p-1)]);
    roiName = [mouse '_roi' int2str(p-1) '.zip'];
    rois = ReadImageJROI(roiName); %read imagej rois
    
    signals_thisplane = signals{p};        
    trace = signals_thisplane(:,:,2); % df/f        
    here = ishere{p};
    responsive = toneresponsive{p};
    
    for j=1:nGroups   
%     j=2;
        groupsok = groups{j}+1;
        cellishere = here(:,groupsok);
        cellhere_group = any(cellishere,2);  
        okcells = cellhere_group & responsive;
        
        these_days = ismember(matrix(:,DAY),groupsok);            
        for c = 1:2
            cntxt = matrix(:,CONTEXT)==(c-1);
            ok = these_days & cntxt;                    
            l = sum(ok);
            matidx = repelem(0:nframes_psth-1,l,1);
            idx = (matidx+repelem(matrix(ok,TONEF(p))- pretone*round(acq),1,nframes_psth))';

            nCells = sum(okcells);
            s = trace(logical(okcells),idx(:))';
            s = reshape(s,nframes_psth,l,nCells);
            ms = median(s,2);
            ms = permute(ms,[1 3 2])';

            msdiff = RES{j,1+(c-1)*2}-RES{j,2+(c-1)*2};
            if p==1
              totake = 1:nCells;
            else
               totake = (size(msdiff,1)-nCells+1):size(msdiff,1);
            end
            msdiffok = msdiff(totake',:); 

            % compute mean in tone evoked resp window
            meantoneresp = mean(msdiffok(:,17:26),2);

            % avg fov for group
            m = [];
            n_days_this_group = length(groups{j});
            for nd=1:n_days_this_group
                m(:,:,nd) = avg_day{groupsok(nd),p};
            end
            avggroup = mean(m,3);       

            %                 cm = colormap(spring);
            subplot(nGroups,2,c+(2*(j-1)));hold on;
%             subplot(2,2,p);hold on;
            imgg = imadjust(adapthisteq(uint16(avggroup)));
            hold on;imagesc(imgg);colormap gray;

            %                 minn = -0.04;
            %                 maxx = 0.15;
            %                 nminmax = length(minn:0.01:maxx);

            [~,wh] = sort(meantoneresp);

            hold off;

%                 R = linspace(cm(1,1),cm(end,1),nCells);
%                 G = linspace(cm(1,2),cm(end,2),nCells);
%                 B = linspace(cm(1,3),cm(end,3),nCells);


            % CLASSIC colormap
            x = linspace(1,nCells,size(cm,1));
            R = interp1(x,cm(:,1),1:nCells);
            G = interp1(x,cm(:,2),1:nCells);
            B = interp1(x,cm(:,3),1:nCells);
            ncm = [R' G' B'];
            
%                 R = linspace(cm(1,1),cm(end,1),nCells);
%                 G = linspace(cm(1,2),cm(end,2),nCells);
%                 B = linspace(cm(1,3),cm(end,3),nCells);
%             ncm = 
            
            selectedrois = rois(okcells)';
            for cel=1:nCells
                yroi = selectedrois{cel}.mnCoordinates(:,1); %freehand rois have the outlines in x-y coordinates
                xroi = selectedrois{cel}.mnCoordinates(:,2);
                patch(yroi,xroi,ncm(wh(cel),:),'EdgeColor','none');
            end                
        end
    end
end
 
    
%% TONE PETH Ziyi Code
%get plane 1 dff signal, plot PETH for neurons of that plane
dff1 = signals{1}; % plane 1
dff1 = dff1(:,:,2); % select dff
toneStart1 = start_tone(:,1);
toneEvoked1 = [];
for i = 1:30
    toneEvoked1(:,:,i) = dff1(:,toneStart1 + i -10); % select 30 frames, 10 frames before tone, 20 after tone
end
tonePETH1 = squeeze(mean(toneEvoked1,2));
imagesc(tonePETH1)






