function results = GetPreprocessedInfo(mouse)
% results = {nFrames_oneplane,signals,matrix,ishere,ks,ks2,ctx,acq,TONEF,REWF,acq,dayss_behav,startfrom};

global info;

if strcmp(mouse,'cd017')
    suite2ppath = 'H:\celine\cd017\suite2p\';
    h5path = 'U:\LabData4\celine\cd017\h5\'; % h5path = 'H:\celine\cd017\h5\';
    sbxpath = 'U:\LabData4\celine\cd017\';
    behavpath = 'U:\LabData4\celine\cd017\behavior\';
end

if strcmp(mouse,'cd036')
    suite2ppath = 'H:\celine\cd036\suite2p\';
    h5path = 'T:\LabData5\cd036\'; 
    sbxpath = 'T:\LabData5\cd036\';
    behavpath = 'T:\LabData5\cd036\behavior\';
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
        
        if (i==1 && k==4)
%             toneFrame = toneFrame+nFrames_add(list_rec(k));            
            toneFramePlane = toneFramePlane+repelem(nFrames_oneplane(list_rec(k),:),nTrials,1);
            rewardFramePlane = rewardFramePlane+repelem(nFrames_oneplane(list_rec(k),:),nTrials,1);
        else
            if ~(k==1 && i==1)
                toneFramePlane = toneFramePlane+repelem(nFrames_oneplane(list_rec(k),:),nTrials,1);
                rewardFramePlane = rewardFramePlane+repelem(nFrames_oneplane(list_rec(k),:),nTrials,1);
            end
        end            
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

TONEF1 = 6; % tone frame plane 1
TONEF2 = 7; % tone frame plane 2
REW1 = 8;
REW2 = 9;


%% Load cell activity + tracking matrix
cd(suite2ppath);
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

results = {nFrames_oneplane,signals,matrix,ishere,ks,ks2,ctx,acq,TONEF,REWF,acq,dayss_behav,startfrom};