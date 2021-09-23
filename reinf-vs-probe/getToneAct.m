filename = {'cd036_002_001','cd036_002_002','cd036_002_003',...
    'cd036_003_001','cd036_003_002','cd036_003_004',...
    'cd036_004_001','cd036_004_002','cd036_004_003'};
behavname = {'cd036_2v1', 'cd036_2v2', 'cd036_2v3',...
'cd036_3v1', 'cd036_3v2', 'cd036_3v3',...
'cd036_4v1', 'cd036_4v2', 'cd036_4v3'};
fileindex = [10 11 12  15 16 17 19 20 21];
days = [3 3 3 4 4 4 5 5 5];

targDff = {[],[]};
foilDff = {[],[]};
for i = 1:length(behavname)
    targReinf = [];
    foilReinf = [];
    targFrame = [];
    foilFrame = [];
    toneWindow = -15:45;
    behavMat = load(['U:\LabData5\cd036\behavior\' behavname{i} '.txt']);
    load([ filename{i} '_dff.mat']);
    
    targFrame = behavMat(behavMat(:,13)==1 & behavMat(:,3)==13454,12);
    foilFrame = behavMat(behavMat(:,13)==1 & behavMat(:,3)==22627,12);  
    for j = 1:length(targFrame)
        targReinf(:,:,j) = sessionDff(:,targFrame(j)/2+toneWindow);
        foilReinf(:,:,j) = sessionDff(:,foilFrame(j)/2+toneWindow);
    end
    targDff{1} = cat(3,targDff{1}, targReinf);
    foilDff{1} = cat(3,foilDff{1}, foilReinf);
    
    targReinf = [];
    foilReinf = [];
    targFrame = behavMat(behavMat(:,13)==0 & behavMat(:,3)==13454,12);
    foilFrame = behavMat(behavMat(:,13)==0 & behavMat(:,3)==22627,12);
    if ~isempty(targFrame)
        for j = 1:length(targFrame)
            targReinf(:,:,j) = sessionDff(:,targFrame(j)/2+toneWindow);
            foilReinf(:,:,j) = sessionDff(:,foilFrame(j)/2+toneWindow);
        end
        targDff{2} = cat(3,targDff{2}, targReinf);
        foilDff{2} = cat(3,foilDff{2}, foilReinf);
    end
end

%%
baselineFrame = 1:15;
for i=1:2
    targDff{i} = targDff{i} - repmat(mean(targDff{i}(:,baselineFrame,:),2),1,size(targDff{i},2),1);
    foilDff{i} = foilDff{i} - repmat(mean(foilDff{i}(:,baselineFrame,:),2),1,size(targDff{i},2),1);
end

%%

filename = {'cd017_004_001','cd017_004_002','cd017_004_003','cd017_004_004',...
    'cd017_005_002','cd017_005_003','cd017_005_004','cd017_005_005',...
    'cd017_006_001','cd017_006_002','cd017_006_003','cd017_006_004'};
behavname = {'cd017_4v1', 'cd017_4v2', 'cd017_4v3','cd017_4v4',...
'cd017_5v1', 'cd017_5v2', 'cd017_5v3','cd017_5v4',...
'cd017_6v1', 'cd017_6v2', 'cd017_6v3','cd017_6v4'};
fileindex = [17 18 19 20  22 23 24 25  27 28 29 30];
days = [4 4 4 4 5 5 5 5 6 6 6 6];

targDff = {[],[]};
foilDff = {[],[]};
for i = 1:length(behavname)
    targReinf = [];
    foilReinf = [];
    targFrame = [];
    foilFrame = [];
    toneWindow = -15:45;
    behavMat = load(['V:\LabData4\celine\cd017\behavior\' behavname{i} '.txt']);
    load([ filename{i} '_dff.mat']);
    
    targFrame = behavMat(behavMat(:,13)==1 & behavMat(:,3)==9514,12);
    foilFrame = behavMat(behavMat(:,13)==1 & behavMat(:,3)==16000,12);  
    for j = 1:length(targFrame)
        targReinf(:,:,j) = sessionDff(:,targFrame(j)/2+toneWindow);
        foilReinf(:,:,j) = sessionDff(:,foilFrame(j)/2+toneWindow);
    end
    targDff{1} = cat(3,targDff{1}, targReinf);
    foilDff{1} = cat(3,foilDff{1}, foilReinf);
    
    targReinf = [];
    foilReinf = [];
    targFrame = behavMat(behavMat(:,13)==0 & behavMat(:,3)==9514,12);
    foilFrame = behavMat(behavMat(:,13)==0 & behavMat(:,3)==16000,12);
    if ~isempty(targFrame)
        for j = 1:length(targFrame)
            targReinf(:,:,j) = sessionDff(:,targFrame(j)/2+toneWindow);
            foilReinf(:,:,j) = sessionDff(:,foilFrame(j)/2+toneWindow);
        end
        targDff{2} = cat(3,targDff{2}, targReinf);
        foilDff{2} = cat(3,foilDff{2}, foilReinf);
    end
end
save('cd017_stimDff.mat','targDff','foilDff');
%%
baselineFrame = 1:15;
for i=1:2
    targDff{i} = targDff{i} - repmat(mean(targDff{i}(:,baselineFrame,:),2),1,size(targDff{i},2),1);
    foilDff{i} = foilDff{i} - repmat(mean(foilDff{i}(:,baselineFrame,:),2),1,size(targDff{i},2),1);
end
save('cd017_stimDffCorrected.mat','targDff','foilDff');

