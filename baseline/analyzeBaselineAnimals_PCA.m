clear;
experimental = {'cd017','cd036','cd037','cd042','cd044'};
for i = 1:length(experimental)
    expDatapath{i} = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' experimental{i} '\corrected'];
    tempData = load([expDatapath{i} '\' experimental{i} '_baseline_data_pca90.mat']);
    expData{i} = tempData;
end

ctrl = {'cd019','cd041'};
for i = 1:length(ctrl)
    ctrlDatapath{i} = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' ctrl{i} '\corrected'];
    tempData = load([ctrlDatapath{i} '\' ctrl{i} '_baseline_data_pca90.mat']);
    ctrlData{i} = tempData;
end
%%
varExpX_exp = []; rangeDay = 8;
for i = 1:length(experimental)
    varExpX = expData{i}.varExpX;
    varExpX_expAni = [];
    if strcmp(experimental{i},'cd017')
        varExpX(6,:) = nan; varExpX(:,6) = nan; 
        varExpX(16,:) = nan; varExpX(:,16) = nan; 
    elseif strcmp(experimental{i},'cd037')
        varExpX(12,:) = nan; varExpX(:,12) = nan; 
    elseif strcmp(experimental{i},'cd044')
        varExpX(17,:) = nan; varExpX(:,17) = nan; 
    end
    nDays = size(expData{i}.varExpX,1);
    for j = 3:nDays-2
        %tempVarExp = zeros(rangeDay*2+1,1);
        tempIdx = (j-rangeDay):1:(j+rangeDay); newTempIdx = tempIdx; 
        newTempIdx(tempIdx<1) = 1; 
        newTempIdx(tempIdx>nDays) = nDays; 
        tempVarExp = varExpX(j,newTempIdx) / varExpX(j,j);
        tempVarExp(tempIdx<1) = nan; 
        tempVarExp(tempIdx>nDays) = nan; 
        varExpX_expAni = cat(1, varExpX_expAni, tempVarExp);
    end
    varExpX_exp = cat(1,varExpX_exp,nanmean(varExpX_expAni));
end

varExpX_ctrl = []; rangeDay = 8;
for i = 1:length(ctrl)
    nDays = size(ctrlData{i}.varExpX,1);
    varExpX_ctrlAni = [];
    for j = 3:nDays-2
        %tempVarExp = zeros(rangeDay*2+1,1);
        tempIdx = (j-rangeDay):1:(j+rangeDay); newTempIdx = tempIdx; 
        newTempIdx(tempIdx<1) = 1; newTempIdx(tempIdx>nDays) = nDays; 
        tempVarExp = ctrlData{i}.varExpX(j,newTempIdx) / ctrlData{i}.varExpX(j,j);
        tempVarExp(tempIdx<1) = nan; tempVarExp(tempIdx>nDays) = nan; 
        varExpX_ctrlAni = cat(1, varExpX_ctrlAni, tempVarExp);
    end
    varExpX_ctrl = cat(1,varExpX_ctrl,nanmean(varExpX_ctrlAni));
end

  
dayAxis = -rangeDay:rangeDay;
yMean = nanmean(varExpX_exp,1); ySEM = [];
for i = 1:size(varExpX_exp,2)
    ySEM(i) = nanstd(varExpX_exp(:,i),0,1) ./ sqrt(sum(~isnan(varExpX_exp(:,i)),1));
end
figure; hold on;
f = fn_plotFillErrorbar(dayAxis,yMean, ySEM,matlabColors(1),'LineStyle','none');
f.FaceAlpha = 0.3;
p1 = plot(dayAxis,yMean,'Color',matlabColors(1),'LineWidth',2);


yMean = nanmean(varExpX_ctrl,1); ySEM = [];
for i = 1:size(varExpX_ctrl,2)
    ySEM(i) = nanstd(varExpX_ctrl(:,i),0,1) ./ sqrt(sum(~isnan(varExpX_ctrl(:,i)),1));
end
f = fn_plotFillErrorbar(dayAxis,yMean, ySEM,matlabColors(2),'LineStyle','none');
f.FaceAlpha = 0.3;
p2 = plot(dayAxis,yMean,'Color',matlabColors(2),'LineWidth',2);
xlabel('Diff Days'); ylabel('Similarity'); 
xlim([-8 8]);legend([p1 p2],'behav','ctrl')

figure; hold on;
for i = 1:size(varExpX_exp,1)
    plot(dayAxis,varExpX_exp(i,:),'Color',matlabColors(1)*0.2 + [1 1 1]*0.8); 
end
plot(dayAxis, nanmean(varExpX_exp,1),'Color',matlabColors(1),'LineWidth',2)

for i = 1:size(varExpX_ctrl,1)
    plot(dayAxis,varExpX_ctrl(i,:),'Color',matlabColors(2)*0.2 + [1 1 1]*0.8); 
end
plot(dayAxis, nanmean(varExpX_ctrl,1),'Color',matlabColors(2),'LineWidth',2)
xlabel('Diff Days'); ylabel('Similarity'); xlim([-8 8]);

%% plot 1.2 - aligned to day 1
dayRange =10; varExpX1_exp = zeros(length(experimental),dayRange);
for i = 1:length(experimental)
    varExpX = expData{i}.varExpX(1,1:dayRange);
    varExpX = varExpX ./ varExpX(1);
    varExpX1_exp(i,:) = varExpX;
end
dayRange =10;varExpX1_ctrl = zeros(length(ctrl),dayRange);
for i = 1:length(ctrl)
    varExpX = ctrlData{i}.varExpX(1,1:dayRange);
    varExpX = varExpX ./ varExpX(1);
    varExpX1_ctrl(i,:) = varExpX; 
end
dayAxis = 1:dayRange; yMean = nanmean(varExpX1_exp,1);
ySEM = nanstd(varExpX1_exp,0,1) ./ sqrt(size(varExpX1_exp,1));
figure; hold on;
f = fn_plotFillErrorbar(dayAxis,yMean, ySEM,matlabColors(1),'LineStyle','none');
f.FaceAlpha = 0.3;
plot(dayAxis,yMean,'Color',matlabColors(1))

dayAxis = 1:dayRange; yMean = nanmean(varExpX1_ctrl,1);
ySEM = nanstd(varExpX1_ctrl,0,1) ./ sqrt(size(varExpX1_ctrl,1));
f = fn_plotFillErrorbar(dayAxis,yMean, ySEM,matlabColors(2),'LineStyle','none');
f.FaceAlpha = 0.3;
plot(dayAxis,yMean,'Color',matlabColors(2))