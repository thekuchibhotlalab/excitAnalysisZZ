%% this section is for cd036
filename = {'cd036_002_001','cd036_002_002','cd036_002_003',...
    'cd036_003_001','cd036_003_002','cd036_003_004',...
    'cd036_004_001','cd036_004_002','cd036_004_003'};
behavname = {'cd036_2v1', 'cd036_2v2', 'cd036_2v3',...
'cd036_3v1', 'cd036_3v2', 'cd036_3v3',...
'cd036_4v1', 'cd036_4v2', 'cd036_4v3'};
fileindex = [10 11 12  15 16 17 19 20 21];
days = [3 4 5];

for i = 1:length(fileindex)
    sessionStart = nFrames_oneplane(fileindex(i),:)+1;
    sessionEnd = nFrames_oneplane(fileindex(i)+1,:);
    
    session_plane1 = signals{1}(:,sessionStart(1):sessionEnd(1),2);
    session_plane2 = signals{2}(:,sessionStart(2):sessionEnd(2),2);
    try
        sessionAct = [session_plane1; session_plane2];
    catch
        sessionAct = [session_plane1; [session_plane2 nan(size(session_plane2,1),1)] ];
    end
    sessionIshere = [ishere{1}; ishere{2}];
    sessionIshere = sum(sessionIshere(:,days),2)==length(days);
    sessionDff = sessionAct(sessionIshere,:);
    save([ filename{i} '_dff.mat'],'sessionDff')

    behavMat = load(['U:\LabData5\cd036\behavior\' behavname{i} '.txt']);
    targFrame{1} = behavMat(behavMat(:,13)==1 & behavMat(:,3)==13454,12);
    targFrame{2} = behavMat(behavMat(:,13)==0 & behavMat(:,3)==13454,12);
    foilFrame{2} = behavMat(behavMat(:,13)==1 & behavMat(:,3)==22627,12);
    foilFrame{2} = behavMat(behavMat(:,13)==0 & behavMat(:,3)==22627,12);
    foilFrame{1} = behavMat(behavMat(:,13)==1 & behavMat(:,3)==22627,12);
    save([filename{i} '_frame.mat'],'targFrame','foilFrame')
    %save([filename{i} '_frame.mat'],'targFrame','foilFrame')

end

%% this section is for cd017

filename = {'cd017_004_001','cd017_004_002','cd017_004_003','cd017_004_004',...
    'cd017_005_002','cd017_005_003','cd017_005_004','cd017_005_005',...
    'cd017_006_001','cd017_006_002','cd017_006_003','cd017_006_004'};
behavname = {'cd017_4v1', 'cd017_4v2', 'cd017_4v3','cd017_4v4',...
'cd017_5v1', 'cd017_5v2', 'cd017_5v3','cd017_5v4',...
'cd017_6v1', 'cd017_6v2', 'cd017_6v3','cd017_6v4'};
fileindex = [17 18 19 20  22 23 24 25  27 28 29 30];
days = [4 5 6];

for i = 1:length(fileindex)
    sessionStart = nFrames_oneplane(fileindex(i),:)+1;
    sessionEnd = nFrames_oneplane(fileindex(i)+1,:);
    
    session_plane1 = signals{1}(:,sessionStart(1):sessionEnd(1),2);
    session_plane2 = signals{2}(:,sessionStart(2):sessionEnd(2),2);
    try
        sessionAct = [session_plane1; session_plane2];
    catch
        sessionAct = [session_plane1; [session_plane2 nan(size(session_plane2,1),1)] ];
    end
    sessionIshere = [ishere{1}; ishere{2}];
    sessionIshere = sum(sessionIshere(:,days),2)==length(days);
    sessionDff = sessionAct(sessionIshere,:);
    save([ filename{i} '_dff.mat'],'sessionDff')

    behavMat = load(['V:\LabData4\celine\cd017\behavior\' behavname{i} '.txt']);
    targFrame{1} = behavMat(behavMat(:,13)==1 & behavMat(:,3)==13454,12);
    targFrame{2} = behavMat(behavMat(:,13)==0 & behavMat(:,3)==13454,12);
    foilFrame{2} = behavMat(behavMat(:,13)==1 & behavMat(:,3)==22627,12);
    foilFrame{2} = behavMat(behavMat(:,13)==0 & behavMat(:,3)==22627,12);
    foilFrame{1} = behavMat(behavMat(:,13)==1 & behavMat(:,3)==22627,12);
    save([filename{i} '_frame.mat'],'targFrame','foilFrame')
    %save([filename{i} '_frame.mat'],'targFrame','foilFrame')

end
