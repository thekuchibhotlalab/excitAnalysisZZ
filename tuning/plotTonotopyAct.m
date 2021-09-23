clear;

allMouse = 'cd044';
savePath = 'D:\labData\excitatory\tuning\masterData\';
load([savePath '\' allMouse '\' 'prePostTuning.mat']);

%% plot1.6: analyis on line of tonotopy changes
colormapIndex = round(linspace(1,64,17));

f1 = plotMapIdx(img, rois, neuronEachPlane,jet,colormapIndex,...
   {preTuning.tuningPeak, postTuning.tuningPeak}, ...
    {responsiveFlagPre, responsiveFlagPost});
h = subplot(2,2,3);
disp('waiting to draw line')
lineRoi = drawline(h);
pause();disp('waiting key press')
linePos= lineRoi.Position;
close(f1);
plotMapIdxLine(img, rois, neuronEachPlane,jet,colormapIndex,...
   {preTuning.tuningPeak, postTuning.tuningPeak}, ...
    {responsiveFlagPre, responsiveFlagPost}, linePos);
saveas(gcf,[savePath mouse '\tonotopy\f2_mapLine.png']);

%% plot1.65: plot the distribution of best peaks on the line
cmap = jet;
cmap = cmap(colormapIndex,:);
projPre = projectCoordPlaneTone(preTuning.tuningPeak, responsiveFlagPre, roiCentroid,linePos,peakSignifPre, neuronPlane);
projPost = projectCoordPlaneTone(postTuning.tuningPeak, responsiveFlagPre, roiCentroid,linePos,peakSignifPost, neuronPlane);
[toneStat, allCounts] = plotTonotopyProj({projPre,projPost}, {peakSignifPre,peakSignifPost},...
    allToneIdx,[targIdx, foilIdx],cmap,false,3);
%%
projPre = projectCoordPlaneToneAct(responsiveTonePre, responsiveFlagPre, roiCentroid,linePos,1:17, neuronPlane);
projPost = projectCoordPlaneToneAct(responsiveTonePost, responsiveFlagPre, roiCentroid,linePos,1:17, neuronPlane);

[toneStat, allCounts] = plotTonotopyProj({projPre,projPost}, {{1:17,1:17},{1:17,1:17}},...
    1:17,[targIdx, foilIdx],cmap,false,3);
saveas(gcf,[savePath mouse '\tonotopy\f3_actProj.png']);
%% plot1.7: quantify the shift of peak tones by the spread of the distribution
[f1, targWidth, foilWidth] = plotTonotopyProjChange(toneStat,{peakSignifPre,peakSignifPost},tones,targIdx,foilIdx);














function [f1] = plotMapIdx(img, rois, neuronEachPlane,C,colormapIndex,plotIndex, plotFlag)
    f1 = figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.04, 0.6, 0.9]);
    if ~iscell(plotIndex)
        plotIndex = {plotIndex};
    end
    if ~iscell(plotFlag)
        plotFlag = {plotFlag};
    end
    nPlanes = length(neuronEachPlane);
    for i = 1:nPlanes
        neuronPlane = [0 cumsum(neuronEachPlane)];
        for k = 1:length(plotIndex)
            subplot(length(plotIndex),2,(k-1)*2+i);             
            imagesc(img{i});colormap gray;hold on;
            ylim([0 size(img{i},1)]);xlim([0 size(img{1},2)]);
            for j = 1:neuronEachPlane(i)
                cellIndex = neuronPlane(i) + j;
                x = rois{i}{j}(:,1); %freehand rois have the outlines in x-y coordinates
                y = rois{i}{j}(:,2); %matlab matrices are inverted so x values are the 2nd column of mn coordinates, and y is the 1st columna
                if plotFlag{k}(cellIndex)
                    %plot(x,y,'.','color',C(colormapIndex(peakIndex(cellIndex)),:),'MarkerSize',1);
                    patch(x,y,C(colormapIndex(plotIndex{k}(cellIndex)),:),'EdgeColor','none');
                else
                    patch(x,y,[0.8 0.8 0.8],'EdgeColor','none');
                    %plot(x,y,'.','color',[0.8 0.8 0.8],'MarkerSize',1);
                end
            end
            title(['Plane' int2str(i)])
            axis off;
        end
    end
end





function [f1] = plotMapIdxLine(img, rois, neuronEachPlane,C,colormapIndex,plotIndex, plotFlag, linePos)
    f1 = figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.04, 0.6, 0.9]);
    if ~iscell(plotIndex)
        plotIndex = {plotIndex};
    end
    if ~iscell(plotFlag)
        plotFlag = {plotFlag};
    end
    nPlanes = length(neuronEachPlane);
    for i = 1:nPlanes
        neuronPlane = [0 cumsum(neuronEachPlane)];
        for k = 1:length(plotIndex)
            h = subplot(length(plotIndex),2,(k-1)*2+i);             
            imagesc(ones(size(img{i})));colormap gray;hold on;
            ylim([0 size(img{i},1)]);xlim([0 size(img{1},2)]);
            for j = 1:neuronEachPlane(i)
                cellIndex = neuronPlane(i) + j;
                x = rois{i}{j}(:,1); %freehand rois have the outlines in x-y coordinates
                y = rois{i}{j}(:,2); %matlab matrices are inverted so x values are the 2nd column of mn coordinates, and y is the 1st columna
                if plotFlag{k}(cellIndex)
                    %plot(x,y,'.','color',C(colormapIndex(peakIndex(cellIndex)),:),'MarkerSize',1);
                    patch(x,y,C(colormapIndex(plotIndex{k}(cellIndex)),:),'EdgeColor','none');
                else
                    patch(x,y,[0.8 0.8 0.8],'EdgeColor','none');
                    %plot(x,y,'.','color',[0.8 0.8 0.8],'MarkerSize',1);
                end
            end
            plot(linePos(:,1),linePos(:,2),'LineWidth',2,'Color',[0 0 0]); 
            set(h, 'Ydir', 'reverse')
            title(['Plane' int2str(i)])
            axis off;
        end
    end
end

function proj = projectCoord(roiCentroid,linePos)
    nNeuron = size(roiCentroid,1);
    lineSlope = (linePos(2,2)-linePos(1,2)) / (linePos(2,1)-linePos(1,1));
    tempVec = [1; lineSlope]; tempVec = tempVec / norm(tempVec);
    proj = roiCentroid * tempVec;
end

function proj = projectCoordPlaneTone(peak, responsiveFlag, roiCentroid,linePos,allToneIdx, neuronPlane)
    nPlanes = length(roiCentroid);
    for i = 1:nPlanes
        nTones = length(allToneIdx{i});
        tempIdx = neuronPlane(i)+1:neuronPlane(i+1);
        tempPeak = peak(tempIdx); 
        tempCent = roiCentroid{i};
        for j = 1:nTones
            cellIdx = (tempPeak == allToneIdx{i}(j));
            proj{i}{j} = projectCoord( tempCent(cellIdx & responsiveFlag(tempIdx),:), linePos);
        end
    end
end


function proj = projectCoordPlaneToneAct(respIdx, responsiveFlag, roiCentroid,linePos,allToneIdx, neuronPlane)
    nPlanes = length(roiCentroid);
    for i = 1:nPlanes
        nTones = length(allToneIdx);
        tempIdx = neuronPlane(i)+1:neuronPlane(i+1);
        tempResp = respIdx(:,tempIdx); 
        tempCent = roiCentroid{i};
        for j = 1:nTones
            cellIdx = (tempResp (allToneIdx(j),:)==1) ;
            proj{i}{j} = projectCoord( tempCent(cellIdx & responsiveFlag(tempIdx),:), linePos);
        end
    end
end

function [toneStat,allCounts] = plotTonotopyProj(proj,peakSignif,allToneIdx,behavToneIdx,cmap,normalized,smoothWindow )
    if nargin == 6
        smoothWindow = 1;
    end
    nData = length(proj);
    nPlanes = length(proj{1});
    toneStat = [];
    figure;
    temp = [];
    for i = 1:nData     
        for j = 1:nPlanes
            nTones = length(proj{i}{j});
            for k = 1:nTones
                temp = [temp;proj{i}{j}{k}];
            end
        end
    end
    edgeMax = max(temp); edgeMax = (ceil(edgeMax/100)+1)*100;
    edgeMin = min(temp); edgeMin = (floor(edgeMin/100)-1)*100;
    edges = edgeMin:50:edgeMax;
    for i = 1:nData     
        for j = 1:nPlanes
            subplot(nData,nPlanes,j+(i-1)*nPlanes); hold on;
            temp = proj{i}{j};
            nTones = length(temp);
            toneStat{i}{j} = zeros(nTones,3);
            for k = 1:nTones
                toneTemp = temp{k};
                toneID = find(peakSignif{i}{j}(k)==allToneIdx);
                if any(toneID == behavToneIdx)
                    linewid = 3;
                else
                    linewid = 1;
                end
                [counts,~] = histcounts(toneTemp,edges);
                allCounts{i}{j}(:,k) = counts;
                if normalized
                    counts = counts / sum(counts);
                end
                edgeMed = edges(1:end-1)/2 + edges(2:end)/2;
                plot(edgeMed, smoothdata(counts,'gaussian',smoothWindow),'Color',cmap(toneID,:),'LineWidth',linewid);
                toneStat{i}{j}(k,1) = median(toneTemp);
                toneStat{i}{j}(k,2) = prctile(toneTemp,75) - prctile(toneTemp,25);
                toneStat{i}{j}(k,3) = length(toneTemp);
            end
            title(['plane' int2str(j)]);
            xlabel('position (pixels)'); ylabel('frequency')
        end
        
    end
end

function [f1, targWidth, foilWidth] = plotTonotopyProjChange(toneStat,peakSignif,toneFreq,targIdx,foilIdx)
%myFun - Description
%
% Syntax: f1 = plotTonotopyProjChange(input)
%
% Long description
    nData = length(toneStat);
    nPlanes = length(toneStat{1});
    nTones = 17;
    f1 = figure;
    set(f1, 'Units', 'Normalized', 'OuterPosition', [0.025, 0.2, 0.95, 0.6]);
    for j = 1:nPlanes
        subplot(1,nPlanes,j)
        preTemp = toneStat{1}{j}; postTemp = toneStat{2}{j};
        preWidth = zeros(nTones,1); postWidth = zeros(nTones,1);
        preWidth(peakSignif{1}{j}) = preTemp(:,2);
        postWidth(peakSignif{2}{j}) = postTemp(:,2);
        plot(1:nTones, preWidth,'LineWidth',2); hold on;
        plot(1:nTones,postWidth,'LineWidth',2);
        ylimm = ylim;
        plot([targIdx targIdx],ylimm,'--g')
        plot([foilIdx foilIdx],ylimm,'--r')
        legend('preLearning','postLearning');
        title(['Plane' int2str(j)]);
        ticksLoc = 1:4:17;
        xlim([1 17])
        xticks(ticksLoc)
        xticklabels(strsplit(num2str(toneFreq(ticksLoc))));
        xlabel('frequency'); ylabel('width (pixels)')

        targWidth{j,1} = preWidth(targIdx);targWidth{j,2} = postWidth(targIdx);
        foilWidth{j,1} = preWidth(foilIdx);foilWidth{j,2} = postWidth(foilIdx);
    end

end