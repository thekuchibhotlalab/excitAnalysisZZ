function getBaselinePCA(mouse,varExpThre,zscoreFlag)
%mouse = 'cd041';
datapath = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' mouse];
savingpath = [datapath '\corrected\']; mkdir(savingpath);
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
tic;
%%
for j = 1:size(baselineS,3) % the first and last session 
    [basis{j}, varExp{j}, proj{j}, covMat{j}] = fn_pca(baselineS(:,:,j),'zscore',zscoreFlag);
end

for i = 1:length(varExp)
    tempComp(i) = find(cumsum(varExp{i})>varExpThre/100,1);
end
nComp = round(mean(tempComp));
varExpX = zeros(size(baselineS,3));
for i = 1:size(baselineS,3)
    b = basis{i}(:,1:nComp);
    for j = 1:size(baselineS,3)
        z = baselineS(:,:,j) - repmat(mean(baselineS(:,:,j),2),1,size(baselineS(:,:,j),2));
        projz = b'*z;
        varz = sum(var(projz,[],2));
        covmat = cov(baselineS(:,:,j)');
        [d,v] = eig(covmat);
        v = diag(v);
        varExpX(i,j) = varz / sum(v);
    end
end

if ~zscoreFlag
    save([savingpath '\' mouse '_baseline_data_pca' int2str(varExpThre) '.mat'],'selectNeuron','baselineS',...
        'basis', 'varExp', 'proj', 'covMat','nComp','varExpX')
else
    save([savingpath '\' mouse '_baseline_data_pca' int2str(varExpThre) '_zscore.mat'],'selectNeuron','baselineS',...
        'basis', 'varExp', 'proj', 'covMat','nComp','varExpX')
end
t = toc;
disp([mouse ' is finished. Time = ' num2str(t,'%.2f') ' sec'])
end