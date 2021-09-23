function ToneEvokedActivityZZ(mouse)
%%

global info;

if strcmp(mouse,'cd017')
    suite2ppath = 'H:\celine\cd017\suite2p\';
    h5path = 'W:\LabData4\celine\cd017\h5\'; % h5path = 'H:\celine\cd017\h5\';
    sbxpath = 'W:\LabData4\celine\cd017\';
    behavpath = 'W:\LabData4\celine\cd017\behavior\';
end

if strcmp(mouse,'cd036')
    suite2ppath = 'H:\celine\cd036\suite2p\';
    h5path = 'H:\celine\cd036\'; 
    sbxpath = 'V:\LabData5\cd036\';
    behavpath = 'V:\LabData5\cd036\behavior\';
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

%sbxread('cd017_000_000',1,1);
sbxread([mouse '_000_000'],1,1);

nFrames = nan(nFiles,1);nFrames_add = nan(nFiles,1);
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
nFrames_oneplane = [[0 0];nFrames_oneplane];

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
        
        toneFramePlane = [toneFrame toneFrame];
        toneFramePlane(logical(mod(toneFrame,2)),:) = [round(toneFrame(logical(mod(toneFrame,2)))/nPlanes) round(toneFrame(logical(mod(toneFrame,2)))/nPlanes)-1];
        toneFramePlane(~mod(toneFrame,2),:) = [toneFrame(~mod(toneFrame,2))/nPlanes toneFrame(~mod(toneFrame,2))/nPlanes];


        if (i==1 && k==4)
%             toneFrame = toneFrame+nFrames_add(list_rec(k));            
            toneFramePlane = toneFramePlane+repelem(nFrames_oneplane(list_rec(k),:),nTrials,1);
        else
            if ~(k==1 && i==1)
                toneFramePlane = toneFramePlane+repelem(nFrames_oneplane(list_rec(k),:),nTrials,1);
            end
        end            
        matrix = [matrix;i*ones(nTrials,1) list_behav(k)*ones(nTrials,1) trialnb contxt resp toneFramePlane];
        
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
cd([suite2ppath 'plane0'])
ishere1 = load('ishere_plane0.mat');
ishere1 = ishere1.ishere;
cd([suite2ppath 'plane1'])
ishere2 = load('ishere_plane1.mat');
ishere2 = ishere2.ishere;
ishere = ishere1;
ishere{2} = ishere2{2};

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
        %signals{i}(:,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i),2) = ...
        %    (submat-median(submat,2))./median(submat,2);%./repmat(max(submat,[],2)-min(submat,[],2),1,size(submat,2));
        tempDff = (submat-quantile(submat,0.5,2))./quantile(submat,0.5,2);
        signals{i}(:,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i),2) = ...
            (submat-quantile(submat,0.5,2))./quantile(submat,0.5,2);%./repmat(max(submat,[],2)-min(submat,[],2),1,size(submat,2));
        
        tempStd = std(tempDff,0,2);
        signals{i}(:,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i),3) = ...
            (tempDff - mean(tempDff,2)) ./ tempStd;%./repmat(max(submat,[],2)-min(submat,[],2),1,size(submat,2));
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

TONEF = [TONEF1 TONEF2];

acq = 30.98;
acq = acq/nPlanes;
pretone = 1; % for PSTH tone, in s
posttone = 3; % for PSTH tone, in s
nframes_psth = pretone*round(acq) + posttone*round(acq);

%% MEDIAN TRACE PER PLANE
% for p=1:nPlanes
%     signals_thisplane = signals{p};
%     trace = signals_thisplane(:,:,2);
%     nCells = size(signals_thisplane,1);     
%     l = size(matrix,1);
%     matidx = repelem(0:nframes_psth-1,l,1);
%     idx = (matidx+repelem(matrix(:,TONEF(p))- pretone*round(acq),1,nframes_psth))';
%     s = trace(:,idx(:))';
%     s = reshape(s,nframes_psth,l*nCells);
%     ms = median(s,2);
%     figure;hold on;
%     plot(ms,'k','linewidth',2);
%     ylim([-0.2 1]);
%     PlotHVLines(pretone*round(acq),'v','color','k','linewidth',1);    
% end


%% MEDIAN TRACE PER PLANE PER DAY
colors = {'g','r'};
for p=1:nPlanes
    signals_thisplane = signals{p};
    trace = signals_thisplane(:,:,2);
%     nCells = size(signals_thisplane,1);   
    figure;hold on;
    for j=1:nDays_behav            
        this_day = matrix(:,DAY)==dayss_behav(j)+1;
        here = ishere{p}(:,dayss_behav(j)+1);
        here(isnan(here)) = 0;
        nCells = sum(here);
        for k=1:2
            ok = this_day & ks2(:,k);                    
            l = sum(ok);
            matidx = repelem(0:nframes_psth-1,l,1);
            idx = (matidx+repelem(matrix(ok,TONEF(p))- pretone*round(acq),1,nframes_psth))';
            s = trace(logical(here),idx(:))';
            s = reshape(s,nframes_psth,l*nCells);
            ms = median(s,2);
            subplot(nDays_behav,2,k+(2*(j-1)));hold on;
            plot(ms,'k','linewidth',2);
%             ylim([-0.2 1]);
            PlotHVLines(pretone*round(acq),'v','color',colors{k},'linewidth',1);                                        
        end
    end
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
%% Find responsive cells
showfig = false;
toneresponsive = cell(nPlanes,1);
for p=1:nPlanes
    signals_thisplane = signals{p};
    nCells = size(signals_thisplane,1); 
    toneresponsive{p} = nan(nCells,1);
    for i=1:nCells
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
            title(['Cell ' num2str(i) ', ' num2str(length(bef)) ', p=' num2str(pttest)]); 
            PlotIntervals([[5 15];[17 26]]);
            pause();
            close(fig);
        end
    end
end

%% LEARNING GROUPS - Target VS Foil
early = 1:3;
acquisition = 4:5;
expression = 6:11;
expert = 12:15;
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
    responsive = toneresponsive{p};
    
    for j=1:nGroups         
        groupsok = groups{j}+1;
        cellishere = here(:,groupsok);
        cellhere_group = any(cellishere,2);  
        okcells = cellhere_group & responsive;
%         okcells = cellhere_group & responsive;
        
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

%% LEARNING GROUPS - Diff Target & Foil

% Plot
if showfig
    fig=figure;hold on;
    for i=1:nGroups
        subplot(nGroups,1,1+(1*(i-1)));hold on;
        msdiff = RES{i,1}-RES{i,2};
        if i==1
            [~,order] = sort(median(msdiff(:,30:60),2));
        end
        PlotColorCurves(msdiff(order,:),[-1 3]);
        PlotHVLines(0,'v','color','w','linewidth',1);
        title([groupNames{i} ' - Diff(T-F)']);
        clim([-0.01 0.06]);
        ylabel('Cells');xlabel('Time (s)');
    end
end

%% LEARNING GROUPS - ALL CONTEXTS - Target VS Foil

RES = cell(nGroups,4);
for p=1:nPlanes
    signals_thisplane = signals{p};
        
    trace = signals_thisplane(:,:,2);        
    here = ishere{p};
    responsive = toneresponsive{p};
    
    for j=1:nGroups         
        groupsok = groups{j}+1;
        cellishere = here(:,groupsok);
        cellhere_group = any(cellishere,2);  
        okcells = cellhere_group & responsive;
%         okcells = cellhere_group & responsive;
        
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
                    [~,order] = sort(median(ms(:,30:60),2));
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
                [~,order] = sort(median(msdiff(:,30:60),2));
            end
            PlotColorCurves(msdiff(order,:),[-1 3]);
            PlotHVLines(0,'v','color','w','linewidth',1);
            title([groupNames{i} ' - ' ctxtNames{c} ' - Diff(T-F)']);
            clim([-0.04 0.15]);
            ylabel('Cells');xlabel('Time (s)');
        end
    end
end



%% LEARNING GROUPS
early = 1:3;
acquisition = 4:5;
expression = 6:11;
expert = 12:15;
groups = {early,acquisition,expression,expert};
nGroups = length(groups);
colors = {'g','r'};
signif = cell(nPlanes,1);
showfig = false;
for p=1:nPlanes
    signals_thisplane = signals{p};
    nCells = size(signals_thisplane,1);      
    signif{p} = nan(nGroups,2*2,nCells);
    for i=1:nCells
        if ~logical(toneresponsive{p}(i)), continue, end
        if showfig
            fig=figure;hold on;
        end
        trace = signals_thisplane(i,:,2);        
        here = ishere{p}(i,:);
        for j=1:nGroups         
            groupsok = groups{j}+1;
            cellishere = find(here==1);
            cellhere_group = cellishere(ismember(cellishere,groupsok));
            these_days = ismember(matrix(:,DAY),cellhere_group);            
            if sum(these_days)==0, continue, end
            for c = 1:2
                cntxt = matrix(:,CONTEXT)==(c-1);
                for k=1:2
                    ok = these_days & cntxt & ks2(:,k);                    
                    l = sum(ok);
                    matidx = repelem(0:nframes_psth-1,l,1);
                    idx = (matidx+repelem(matrix(ok,TONEF(p))- pretone*round(acq),1,nframes_psth))';
                    s = trace(:,idx(:));
                    s = reshape(s,nframes_psth,l);
                    ms = median(s,2);
                    
                    
                    bef=median(s(11:30,:));
                    aft=median(s(36:55,:));
                    [h0,pttest] = ttest(bef(:),aft(:),'alpha',0.01);         
                    signif{p}(j,k+(2*(c-1)),i) = h0;
                    
                    if showfig
                        subplot(nGroups,4,k+(2*(c-1))+(4*(j-1)));hold on;
                        if h0==0, title(['Cell ' num2str(i) ', ' num2str(length(bef)) ', N']); 
                        else, title(['Cell ' num2str(i) ', ' num2str(length(bef)) ', Y']); end                     
                        PlotIntervals([[11 30];[36 55]]);
                        plot(s,'color',[0.8 0.8 0.8]);hold on;
                        plot(ms,'k','linewidth',2);
                        ylim([-0.05 0.1]);
                        PlotHVLines(pretone*round(acq),'v','color',colors{k},'linewidth',1); 
                    end
                end
            end
        end
        if showfig
            pause();
            close(fig);
        end
    end
end 

%% Plot only tone responsive cells

for p=1:nPlanes
    ok = logical(toneresponsive{p});
    nCells = sum(ok);
    trace = signals{p}(ok,:,2);
    here = ishere{p}(ok,:);
    
    for j=1:nGroups         
        groupsok = groups{j}+1;
        cellishere = here(:,groupsok);
        cellhere_group = any(cellishere,2);       
        
        these_days = ismember(matrix(:,DAY),groupsok);   
        for k=1:2
            ok = these_days & ks2(:,k);                    
            l = sum(ok);
            matidx = repelem(0:nframes_psth-1,l,1);
            idx = (matidx+repelem(matrix(ok,TONEF(p))- pretone*round(acq),1,nframes_psth))';
            s = trace(cellhere_group,idx(:));
            s = reshape(s,nframes_psth,l*sum(cellhere_group));
            ms = median(s,2);


            bef=median(s(11:30,:));
            aft=median(s(36:55,:));
            [h0,pttest] = ttest(bef(:),aft(:),'alpha',0.01);         
            signif{p}(j,k+(2*(c-1)),i) = h0;      
            

            if showfig
                subplot(nGroups,4,k+(2*(c-1))+(4*(j-1)));hold on;
                if h0==0, title(['Cell ' num2str(i) ', ' num2str(length(bef)) ', N']); 
                else, title(['Cell ' num2str(i) ', ' num2str(length(bef)) ', Y']); end                     
                PlotIntervals([[11 30];[36 55]]);
                plot(s,'color',[0.8 0.8 0.8]);hold on;
                plot(ms,'k','linewidth',2);
                ylim([-0.05 0.1]);
                PlotHVLines(pretone*round(acq),'v','color',colors{k},'linewidth',1); 
            end
        end
    end
end





%% get 200-400 trials

periodBefAftTone = -15:60;
%periodBefAftTone = 3:12;
toneActTrials = cell(1,2);
selectTrials = 180:375;
for i = 1:nPlanes
    
     for j = 1:length(periodBefAftTone)
         
        toneActTrials{i}(:,:,j) = signals{i}(:,matrix(selectTrials,5+i)+periodBefAftTone(j),3);
        
     end
     
    
    
end
tempAct = [toneActTrials{1};toneActTrials{2}];
actForTCA = tempAct(ishere2plane(:,3)==1,:,:,:);
save('actDay3.mat','actForTCA')


%% get 200-400 trials

periodBefAftTone = -15:60;
%periodBefAftTone = 3:12;
toneActTrials = cell(1,2);
selectTrials = 180:375;
for i = 1:nPlanes
    
    for j = 1:length(periodBefAftTone)

    toneActTrials{i}(:,:,j) = signals{i}(:,matrix(selectTrials,5+i)+periodBefAftTone(j),3);

    end
     

    
end
tempAct = [toneActTrials{1};toneActTrials{2}];
actForTCA = tempAct(ishere2plane(:,3)==1,:,:,:);
save('actDay3.mat','actForTCA')

%% get 200-400 trials

periodBefAftTone = -15:60;
%periodBefAftTone = 3:12;
toneActTrials = cell(1,2);
selectTrials = 180:375;

reinftempcount
reinf
for i = 1:nPlanes
    
    for j = 1:length(periodBefAftTone)
        toneActTrials{i}(:,:,j) = signals{i}(:,matrix(selectTrials,5+i)+periodBefAftTone(j),3);
    end
     

    
end
tempAct = [toneActTrials{1};toneActTrials{2}];
actForTCA = tempAct(ishere2plane(:,3)==1,:,:,:);
save('actDay3.mat','actForTCA')
%%
signals_zscore_days{1} = signals{1}(:,:,2);
signals_zscore_days{1} = signals_zscore_days{1}  - repmat(mean(signals_zscore_days{1},2),[1 size(signals_zscore_days{1},2)]);
signals_zscore_days{1} = signals_zscore_days{1}  ./ repmat(std(signals_zscore_days{1},0,2),[1 size(signals_zscore_days{1},2)]);
signals_zscore_days{2} = signals{2}(:,:,2);
signals_zscore_days{2} = signals_zscore_days{2}  - repmat(mean(signals_zscore_days{2},2),[1 size(signals_zscore_days{2},2)]);
signals_zscore_days{2} = signals_zscore_days{2}  ./ repmat(std(signals_zscore_days{2},0,2),[1 size(signals_zscore_days{2},2)]);

%% make the plots of cells activity acroos days

cmap = redbluecmap;
newCmap = imresize(cmap, [128, 3]);  % original color map contain just 11 colors, this increase it to 64
newCmap = min(max(newCmap, 0), 1);

periodBefAftTone = -15:60;
%periodBefAftTone = 3:12;
toneAct = cell(1,2);
for i = 1:nPlanes
    
     for j = 1:length(periodBefAftTone)
         
        %toneAct{i}(:,:,j) = signals{i}(:,matrix(:,5+i)+periodBefAftTone(j),2);
        toneAct{i}(:,:,j) = signals_zscore_days{i}(:,matrix(:,5+i)+periodBefAftTone(j));
        
     end
     
    
    
    
end

imagesc(squeeze(mean(toneAct{1},2))); caxis([-0.2 0.2]);
days = 2:17;
toneActDay = cell(length(days),2);
targetActDay = cell(length(days),2);
foilActDay = cell(length(days),2);

avgToneActDay = cell(1,2);
avgTargActDay = cell(1,2);
avgFoilActDay = cell(1,2);

avgToneActDayTC = cell(1,2);
avgTargActDayTC = cell(1,2);
avgFoilActDayTC = cell(1,2);

avgReinfActDay = cell(1,2);
avgProbeActDay = cell(1,2);

avgReinfTargActDay = cell(1,2);
avgProbeTargActDay = cell(1,2);
avgReinfFoilActDay = cell(1,2);
avgProbeFoilActDay = cell(1,2);
avgReinfTargActDayTC = cell(1,2);
avgProbeTargActDayTC = cell(1,2);
avgReinfFoilActDayTC = cell(1,2);
avgProbeFoilActDayTC = cell(1,2);

reinfActDay = cell(length(days),2);
probeActDay = cell(length(days),2);

reinfTargetActDay = cell(length(days),2);
probeTargetActDay = cell(length(days),2);
reinfFoilActDay = cell(length(days),2);
probeFoilActDay = cell(length(days),2);
            
%periodAftTone = 12:21;
%periodAftTone = 13:17;
periodAftTone = 19:22;
%periodAftTone = 1;
for i = 1:nPlanes
    for j = 1:length(days)
        for k = 1:length(periodBefAftTone)
            goodTrialFlag = matrix(:,1) == days(j);
            if j == 2
                goodTrialFlag(180:235) = false;
            end
            
            
            toneThisDay = matrix(goodTrialFlag, 5+i);
            targetThisDay = matrix((goodTrialFlag & (matrix(:,5) == 1 | matrix(:,5) == 2)), 5+i);
            foilThisDay = matrix((goodTrialFlag & (matrix(:,5) == 3 | matrix(:,5) == 4)), 5+i);
            
            reinfThisDay = matrix((goodTrialFlag & matrix(:,4)==1), 5+i);
            probeThisDay = matrix((goodTrialFlag & matrix(:,4)==0), 5+i);
            
            reinfTargThisDay = matrix((goodTrialFlag & (matrix(:,5) == 1 | matrix(:,5) == 2)) & matrix(:,4)==1, 5+i);
            probeTargThisDay = matrix((goodTrialFlag & (matrix(:,5) == 1 | matrix(:,5) == 2)) & matrix(:,4)==0, 5+i);
            
            reinfFoilThisDay = matrix((goodTrialFlag & (matrix(:,5) == 3 | matrix(:,5) == 4)) & matrix(:,4)==1, 5+i);
            probeFoilThisDay = matrix((goodTrialFlag & (matrix(:,5) == 3 | matrix(:,5) == 4)) & matrix(:,4)==0, 5+i);
            
%             toneActDay{j,i}(:,:,k) = signals{i}(:,toneThisDay+periodBefAftTone(k),2);
%             targetActDay{j,i}(:,:,k) = signals{i}(:,targetThisDay+periodBefAftTone(k),2);
%             foilActDay{j,i}(:,:,k) = signals{i}(:,foilThisDay+periodBefAftTone(k),2);
            toneActDay{j,i}(:,:,k) = signals_zscore_days{i}(:,toneThisDay+periodBefAftTone(k));
            targetActDay{j,i}(:,:,k) = signals_zscore_days{i}(:,targetThisDay+periodBefAftTone(k));
            foilActDay{j,i}(:,:,k) = signals_zscore_days{i}(:,foilThisDay+periodBefAftTone(k));
            
            %reinfActDay{j,i}(:,:,k) = signals{i}(:,reinfThisDay+periodBefAftTone(k),2);
            %probeActDay{j,i}(:,:,k) = signals{i}(:,probeThisDay+periodBefAftTone(k),2);
            reinfActDay{j,i}(:,:,k) = signals_zscore_days{i}(:,reinfThisDay+periodBefAftTone(k));
            probeActDay{j,i}(:,:,k) = signals_zscore_days{i}(:,probeThisDay+periodBefAftTone(k));
            
            %reinfTargetActDay{j,i}(:,:,k) = signals{i}(:,reinfTargThisDay+periodBefAftTone(k),3);
            %probeTargetActDay{j,i}(:,:,k) = signals{i}(:,probeTargThisDay+periodBefAftTone(k),3);
            %reinfFoilActDay{j,i}(:,:,k) = signals{i}(:,reinfFoilThisDay+periodBefAftTone(k),3);
            %probeFoilActDay{j,i}(:,:,k) = signals{i}(:,probeFoilThisDay+periodBefAftTone(k),3);
            reinfTargetActDay{j,i}(:,:,k) = signals_zscore_days{i}(:,reinfTargThisDay+periodBefAftTone(k));
            probeTargetActDay{j,i}(:,:,k) = signals_zscore_days{i}(:,probeTargThisDay+periodBefAftTone(k));
            reinfFoilActDay{j,i}(:,:,k) = signals_zscore_days{i}(:,reinfFoilThisDay+periodBefAftTone(k));
            probeFoilActDay{j,i}(:,:,k) = signals_zscore_days{i}(:,probeFoilThisDay+periodBefAftTone(k));
            
           
        end
        
        avgToneActDay{i}(:,j) = squeeze(mean(mean(toneActDay{j,i}(:,:,periodAftTone),2),3));
        
        avgToneActDayTC{i}(:,:,j) = squeeze(mean(toneActDay{j,i}(:,:,:),2));
        
        avgReinfTargActDayTC{i}(:,:,j) = squeeze(mean(reinfTargetActDay{j,i}(:,:,:),2));
        avgProbeTargActDayTC{i}(:,:,j) = squeeze(mean(probeTargetActDay{j,i}(:,:,:),2));
        avgReinfFoilActDayTC{i}(:,:,j) = squeeze(mean(reinfFoilActDay{j,i}(:,:,:),2));
        avgProbeFoilActDayTC{i}(:,:,j) = squeeze(mean(probeFoilActDay{j,i}(:,:,:),2));
        
        avgTargActDay{i}(:,j) = squeeze(mean(mean(targetActDay{j,i}(:,:,periodAftTone),2),3));
        avgTargActDayTC{i}(:,:,j) = squeeze(mean(targetActDay{j,i}(:,:,:),2));
        
        avgFoilActDay{i}(:,j) = squeeze(mean(mean(foilActDay{j,i}(:,:,periodAftTone),2),3));
        avgFoilActDayTC{i}(:,:,j) = squeeze(mean(foilActDay{j,i}(:,:,:),2));
        
        avgReinfActDay{i}(:,j) = squeeze(mean(mean(reinfActDay{j,i}(:,:,periodAftTone),2),3));
        avgProbeActDay{i}(:,j) = squeeze(mean(mean(probeActDay{j,i}(:,:,periodAftTone),2),3));
        
        avgReinfTargActDay{i}(:,j) = squeeze(mean(mean(reinfTargetActDay{j,i}(:,:,periodAftTone),2),3));
        avgProbeTargActDay{i}(:,j) = squeeze(mean(mean(probeTargetActDay{j,i}(:,:,periodAftTone),2),3));
        avgReinfFoilActDay{i}(:,j) = squeeze(mean(mean(reinfFoilActDay{j,i}(:,:,periodAftTone),2),3));
        avgProbeFoilActDay{i}(:,j) = squeeze(mean(mean(probeFoilActDay{j,i}(:,:,periodAftTone),2),3));
        
        avgTargActDayTC{i}(:,:,j) = squeeze(mean(targetActDay{j,i}(:,:,:),2));
        avgTargActDayTC{i}(:,:,j) = squeeze(mean(targetActDay{j,i}(:,:,:),2));
        
    end
end
toneRespFlag = [toneresponsive{1};toneresponsive{2}];
avgToneActDay = [avgToneActDay{1}; avgToneActDay{2}];
avgTargActDay = [avgTargActDay{1}; avgTargActDay{2}];
avgFoilActDay = [avgFoilActDay{1}; avgFoilActDay{2}];
avgToneActDayTC = [avgToneActDayTC{1}; avgToneActDayTC{2}];
avgTargActDayTC = [avgTargActDayTC{1}; avgTargActDayTC{2}];
avgFoilActDayTC = [avgFoilActDayTC{1}; avgFoilActDayTC{2}];
avgReinfTargDayTC = [avgReinfTargActDayTC{1}; avgReinfTargActDayTC{2}];
avgProbeTargDayTC = [avgProbeTargActDayTC{1}; avgProbeTargActDayTC{2}];
avgReinfFoilDayTC = [avgReinfFoilActDayTC{1}; avgReinfFoilActDayTC{2}];
avgProbeFoilDayTC = [avgProbeFoilActDayTC{1}; avgProbeFoilActDayTC{2}];

ishereFlag = [ishere{1}; ishere{2}];

%avgToneActDay(~ishereFlag(:,2:17)) = nan;
%avgTargActDay(~ishereFlag(:,2:17)) = nan;
%avgFoilActDay(~ishereFlag(:,2:17)) = nan;
%% save things for TCA
neuronPlane1 = size(signals{1},1);
neuronPlane2 = size(signals{2},1);
tempAct = zeros(neuronPlane1+neuronPlane2, 76, 16,4);
for i = 1:16
tempAct(1:neuronPlane1,:,i,1) = squeeze(mean(reinfTargetActDay{i,1},2));
tempAct(1:neuronPlane1,:,i,2) = squeeze(mean(reinfFoilActDay{i,1},2));
tempAct(1:neuronPlane1,:,i,3) = squeeze(mean(probeTargetActDay{i,1},2));
tempAct(1:neuronPlane1,:,i,4) = squeeze(mean(probeFoilActDay{i,1},2));
tempAct(neuronPlane1+1:end,:,i,1) = squeeze(mean(reinfTargetActDay{i,2},2));
tempAct(neuronPlane1+1:end,:,i,2) = squeeze(mean(reinfFoilActDay{i,2},2));
tempAct(neuronPlane1+1:end,:,i,3) = squeeze(mean(probeTargetActDay{i,2},2));
tempAct(neuronPlane1+1:end,:,i,4) = squeeze(mean(probeFoilActDay{i,2},2));
end
%selectdays = [2:14 16 17];
selectdays = 3:17;
ishere2plane = [ishere{1};ishere{2}];
cellAllDay = sum(ishere2plane(:,selectdays)==1,2) == length(selectdays);
actForTCA = tempAct(cellAllDay,:,2:16,:);
%actForTCA = tempAct(cellAllDay,:,[1:13 15 16],:);
%actForTCA = actForTCA(:,:,2:end,:);
%save('actForTCA_zscore017_allday.mat','actForTCA')


%% smooth for visualization

%avgToneActDay = interpSmoothStuff(avgToneActDay);
%avgTargActDay = interpSmoothStuff(avgTargActDay);
%avgFoilActDay = interpSmoothStuff(avgFoilActDay);

%% sort it

diffAct =  abs(avgTargActDay - avgFoilActDay);
sortIndex = sortStuff(avgToneActDay);

figure;
subplot(1,2,1)
imagesc(diffAct);
caxis([0 0.1])
subplot(1,2,2)
imagesc(diffAct(sortIndex));
caxis([0 0.1])

figure;
subplot(1,4,1)
imagesc(removenan(avgToneActDay(sortIndex,:)))
caxis([-0.05 0.10])
colorbar
subplot(1,4,2)
imagesc(removenan(avgTargActDay(sortIndex,:)))
caxis([-0.05 0.10])
subplot(1,4,3)
imagesc(removenan(avgFoilActDay(sortIndex,:)))
caxis([-0.05 0.10])
subplot(1,4,4)
imagesc(removenan(diffAct(sortIndex,:)))
caxis([0 0.1])

figure;
subplot(1,4,1)
imagesc(removenan(avgToneActDay(sortIndex(logical(toneRespFlag)),:)))
caxis([-0.05 0.10])
colorbar
subplot(1,4,2)
imagesc(removenan(avgTargActDay(sortIndex(logical(toneRespFlag)),:)))
caxis([-0.05 0.10])
subplot(1,4,3)
imagesc(removenan(avgFoilActDay(sortIndex(logical(toneRespFlag)),:)))
caxis([-0.05 0.10])
subplot(1,4,4)
imagesc(removenan(diffAct(sortIndex(logical(toneRespFlag)),:)))
caxis([0 0.1])


figure;
subplot(1,4,1)
imagesc(avgToneActDay(sortIndex,:))
caxis([-0.05 0.10])
colorbar
subplot(1,4,2)
imagesc(avgTargActDay(sortIndex,:))
caxis([-0.05 0.10])
subplot(1,4,3)
imagesc(avgFoilActDay(sortIndex,:))
caxis([-0.05 0.10])
subplot(1,4,4)
imagesc(diffAct(sortIndex,:))
caxis([0 0.1])

figure;
subplot(1,4,1)
imagesc(avgToneActDay(sortIndex(logical(toneRespFlag)),:))
caxis([-0.05 0.10])
colorbar
subplot(1,4,2)
imagesc(avgTargActDay(sortIndex(logical(toneRespFlag)),:))
caxis([-0.05 0.10])
subplot(1,4,3)
imagesc(avgFoilActDay(sortIndex(logical(toneRespFlag)),:))
caxis([-0.05 0.10])
subplot(1,4,4)
imagesc(diffAct(sortIndex(logical(toneRespFlag)),:))
caxis([0 0.1])

%%cd017 old
plotNeuronAcrossDays(595,toneActDay,sortIndex)
plotNeuronAcrossDays(602,toneActDay,sortIndex)
plotNeuronAcrossDays(603,toneActDay,sortIndex)
plotNeuronAcrossDays(604,toneActDay,sortIndex)
plotNeuronAcrossDays(632,toneActDay,sortIndex)
%%
plotNeuronAcrossDays(618,toneActDay,1:651)

%%

figure;
subplot(1,4,1)
imagesc(reshape(avgToneActDayTC(sortIndex,:,:),651,[]))
caxis([-0.05 0.10])
colorbar
subplot(1,4,2)
imagesc(reshape(avgTargActDayTC(sortIndex,:,:),651,[]))
caxis([-0.05 0.10])
subplot(1,4,3)
imagesc(reshape(avgFoilActDayTC(sortIndex,:,:),651,[]))
caxis([-0.05 0.10])

figure;
a = avgToneActDayTC(sortIndex,:,:);
for i = 1:16
    subplot(4,4,i)
    plot(a(625,:,i))
    xticks([0 10 20 30]);
    xticklabels({'-0.67','0','0.67','1.33'})
    ylabel('df/f');
    xlabel('time(s)')
    ylim([-0.1 0.8])
    xlim([0 30])
    title(['day' num2str(i)])
end
%% real plots start here

avgReinfActDay = [avgReinfActDay{1}; avgReinfActDay{2}];
avgProbeActDay = [avgProbeActDay{1}; avgProbeActDay{2}];
        
avgReinfTargActDay = [avgReinfTargActDay{1}; avgReinfTargActDay{2}];
avgProbeTargActDay = [avgProbeTargActDay{1}; avgProbeTargActDay{2}];
avgReinfFoilActDay = [avgReinfFoilActDay{1}; avgReinfFoilActDay{2}];
avgProbeFoilActDay = [avgProbeFoilActDay{1}; avgProbeFoilActDay{2}];

%% plot across days for target and foil
avgReinfActDay(~ishereFlag(:,2:17)) = nan;
avgProbeActDay(~ishereFlag(:,2:17)) = nan;

avgReinfTargActDay(~ishereFlag(:,2:17)) = nan;
avgProbeTargActDay(~ishereFlag(:,2:17)) = nan;
avgReinfFoilActDay(~ishereFlag(:,2:17)) = nan;
avgProbeFoilActDay(~ishereFlag(:,2:17)) = nan;


a = sum(~ishereFlag(:,3:17),2);
notNanFlag = a==0;
%%
colorlimm = [-2 2];

figure;
subplot(1,4,1)
temp1 = avgReinfActDay((notNanFlag & logical(toneRespFlag)),:);
sortIndex = sortStuff(temp1);
imagesc(interpSmoothStuff(temp1(sortIndex,:)))
caxis(colorlimm)
%colorbar
%ylabel('neurons')
%xlabel('days')
%title('all tone')
set(gcf,'Position',[100 100 150 1000])
colormap(newCmap)
set(gca,'YTickLabel',[]);
set(gca,'YTick',[]);
set(gca,'XTickLabel',[]);
set(gca,'XTick',[]);

subplot(1,4,2)
temp2 = avgReinfTargActDay((notNanFlag & logical(toneRespFlag)),:);
%temp = temp(notNanFlag,:);
imagesc(interpSmoothStuff(temp2(sortIndex,:)))
caxis(colorlimm)

ylabel('neurons')
xlabel('days')
title('target')

subplot(1,4,3)

temp3 = avgReinfFoilActDay((notNanFlag & logical(toneRespFlag)),:);
%temp = temp(notNanFlag,:);
imagesc(interpSmoothStuff(temp3(sortIndex,:)))
caxis(colorlimm)

ylabel('neurons')
xlabel('days')
title('foil')

subplot(1,4,4)
%temp = temp(notNanFlag,:);
imagesc(interpSmoothStuff((temp2(sortIndex,:) - temp3(sortIndex,:))))
caxis(colorlimm)
colorbar

ylabel('neurons')
xlabel('days')
title('target-foil')
colormap(newCmap)
%suptitle('reinforced')
%%
figure;
subplot(1,4,1)
temp1 = avgProbeActDay((notNanFlag & logical(toneRespFlag)),:);
sortIndex = sortStuff(temp1);
%sortIndex = sortStuff(temp1);
imagesc(interpSmoothStuff(temp1(sortIndex,:)))
caxis(colorlimm)

colorbar
ylabel('neurons')
xlabel('days')
title('all tone')

subplot(1,4,2)

temp2 = avgProbeTargActDay((notNanFlag & logical(toneRespFlag)),:);
%temp = temp(notNanFlag,:);
imagesc(interpSmoothStuff(temp2(sortIndex,:)))
caxis(colorlimm)

ylabel('neurons')
xlabel('days')
title('target')

subplot(1,4,3)
temp3 = avgProbeFoilActDay((notNanFlag & logical(toneRespFlag)),:);
%temp = temp(notNanFlag,:);
imagesc(interpSmoothStuff(temp3(sortIndex,:)))
caxis(colorlimm)
ylabel('neurons')
xlabel('days')
title('foil')

subplot(1,4,4)

%temp = temp(notNanFlag,:);
imagesc(interpSmoothStuff((temp2(sortIndex,:) - temp3(sortIndex,:))))
caxis(colorlimm)

colorbar
ylabel('neurons')
xlabel('days')
title('target-foil')
colormap(newCmap)
%suptitle('probe')



%% not separating reinf vs. probe


avgToneActDayTC(repmat(reshape(~ishereFlag(:,2:17),651,1,16),[1 31 1])) = nan;
avgTargActDayTC(repmat(reshape(~ishereFlag(:,2:17),651,1,16),[1 31 1])) = nan;
avgFoilActDayTC(repmat(reshape(~ishereFlag(:,2:17),651,1,16),[1 31 1])) = nan;
a = sum(~ishereFlag(:,2:17),2);
notNanFlag = a<3;

figure;
subplot(1,4,1)
temp1 = squeeze(mean(avgToneActDayTC((notNanFlag & logical(toneRespFlag)),periodAftTone,:),2));
sortIndex = sortStuff(temp1);
imagesc(interpSmoothStuff(temp1(sortIndex,:)))
caxis([-0.2 0.2])
colorbar
ylabel('neurons')
xlabel('days')
title('all tone')

subplot(1,4,2)
temp2 = squeeze(mean(avgTargActDayTC((notNanFlag & logical(toneRespFlag)),periodAftTone,:),2));
%temp = temp(notNanFlag,:);
imagesc(interpSmoothStuff(temp2(sortIndex,:)))
caxis([-0.2 0.2])
ylabel('neurons')
xlabel('days')
title('target')

subplot(1,4,3)
temp3 = squeeze(mean(avgFoilActDayTC((notNanFlag & logical(toneRespFlag)),periodAftTone,:),2));
%temp = temp(notNanFlag,:);
imagesc(interpSmoothStuff(temp3(sortIndex,:)))
caxis([-0.2 0.2])
ylabel('neurons')
xlabel('days')
title('foil')

subplot(1,4,4)
%temp = temp(notNanFlag,:);
imagesc(interpSmoothStuff((temp2(sortIndex,:) - temp3(sortIndex,:))))
caxis([-0.20 0.20])
colorbar
ylabel('neurons')
xlabel('days')
title('target-foil')
colormap(newCmap)
%%

figure;
subplot(1,4,1)
imagesc(squeeze(mean(avgToneActDayTC(sortIndex,periodAftTone,:),2)))
caxis([-0.05 0.10])
colorbar

subplot(1,4,2)
imagesc(squeeze(mean(avgTargActDayTC(sortIndex,periodAftTone,:),2)))
caxis([-0.05 0.10])
subplot(1,4,3)
imagesc(squeeze(mean(avgFoilActDayTC(sortIndex,periodAftTone,:),2)))
caxis([-0.05 0.10])

%% 
%plotNeuronAcrossDays(319, toneActDay,1:651)
%plotNeuronAcrossDays(320, toneActDay,1:651)
%plotNeuronAcrossDays(321, toneActDay,1:651)
%plotNeuronAcrossDays(631, toneActDay,1:651)
%plotNeuronAcrossDays(616, toneActDay,1:651)
%plotNeuronAcrossDays(629, toneActDay,1:651)
%plotNeuronAcrossDaysRP(78,avgReinfTargDayTC, avgReinfFoilDayTC)
%plotNeuronAcrossDaysRP(101,avgReinfTargDayTC, avgReinfFoilDayTC)
%plotNeuronAcrossDaysRP(101,avgProbeTargDayTC, avgProbeFoilDayTC)
%plotNeuronAcrossDaysRP(251,avgReinfTargDayTC, avgReinfFoilDayTC)
%plotNeuronAcrossDaysRP(251,avgProbeTargDayTC, avgProbeFoilDayTC)

%plotNeuronAcrossDaysRP(376,avgReinfTargDayTC, avgReinfFoilDayTC)
%plotNeuronAcrossDaysRP(376,avgProbeTargDayTC, avgProbeFoilDayTC) % good cell for cd017

%plotNeuronAcrossDaysRP(131,avgReinfTargDayTC, avgReinfFoilDayTC,[-1 1])
plotNeuronAcrossDaysRP(618,avgReinfTargDayTC, avgReinfFoilDayTC,[-1 1])
end
%% functions

function sortIndex = sortStuff(mat)
    
[~, index] = max(mat,[],2);
[~, sortIndex] = sort(index);

end

function mat = interpSmoothStuff(mat)

for i = 1:size(mat,1)
    x = mat(i,:);
    nanx = isnan(x);
    t    = 1:numel(x);
    try
        x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx));
    catch
        disp('neuron not there')
    end
    mat(i,:) = smoothdata(x,2,'movmean',3);
end
%a = sum(mat,2);
%mat(isnan(a),:) = [];

end

function mat = removenan(mat)
a = sum(mat,2);
mat(isnan(a),:) = [];
end

function plotNeuronAcrossDays(neuronNum,toneActDay, sortIndex)
    figure;
    for i = 1:16
        temp = [toneActDay{i,1};toneActDay{i,2}]; 
        act = squeeze(mean(temp(sortIndex==neuronNum, :,:),2));
        subplot(1,16,i)
        
        plot(act)
        ylim([-0.1 0.3])
        xlim([0 31])
    end

end


function plotNeuronAcrossDaysRP(neuronNum,toneActDay1, toneActDay2,ylimm)
    %figure;
    act11 = [];
    act22 = [];
    for i = 1:16
        
        act1 = squeeze(toneActDay1(neuronNum, :,i));
        act11 = [act11 act1];
        act2 = squeeze(toneActDay2(neuronNum, :,i));
        act22 = [act22 act2];
        %subplot(1,16,i)
        
%         plot(act1, 'LineWidth',2)
%         hold on; plot(act2, 'LineWidth',2)
%         hold on; plot([10 10], [-0.1 0.4], 'Color', [0.8 0.8 0.8]);
%         ylim([-0.1 0.4])
%         xlim([0 31])
        
        

%         set(gca,'YTickLabel',[]);
%         set(gca,'YTick',[]);
%         set(gca,'XTickLabel',[]);
%         set(gca,'XTick',[]);
    end
    figure;
    plot(act11, 'LineWidth',2); hold on; plot(act22, 'LineWidth',2)
    for i = 1:16
        plot([1+(i-1)*31 1+(i-1)*31], [-0.1 0.4], 'Color', [0 0 0], 'LineWidth',1.5);
        plot([10+(i-1)*31 10+(i-1)*31], [-0.1 0.4], 'Color', [0.8 0.8 0.8]);
    end
    ylim(ylimm)
    xlim([1 31*16])
    set(gca,'YTickLabel',[]);
    set(gca,'YTick',[]);
    set(gca,'XTickLabel',[]);
    set(gca,'XTick',[]);
end

