function [f] = func_tonotopyPlot(img, rois,C,colormapIndex,plotIndex, varargin)

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('title', [])
p.addParameter('axes', [])
p.addParameter('planeFlag', [])
p.addParameter('plotFlag', [])
p.addParameter('colorMod', [])
p.parse(varargin{:});

nPlot = length(plotIndex);


if isempty(p.Results.axes)
    f = figure; 
    set(f, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.04, 0.6, 0.9]);
    for i = 1:nPlot
        h{i} = subplot_tight(1,nPlot,i);
    end
else 
    h = p.Results.axes;
    f = ancestor(h{1},'figure');
end


for k = 1:nPlot
    set(f,'CurrentAxes',h{k})
    xmax = 0; ymax = 0; xmin = 1000; ymin = 1000;
    imagesc(h{k},img);colormap gray;hold on;
    
    if ~isempty(p.Results.planeFlag)
        tempPlotIndex = plotIndex{k}(p.Results.planeFlag);
        if ~isempty(p.Results.plotFlag); tempPlotFlag = p.Results.plotFlag{k}(p.Results.planeFlag); end
    else 
        tempPlotIndex = plotIndex{k};
        if ~isempty(p.Results.plotFlag); tempPlotFlag = p.Results.plotFlag{k}; end
    end
    nNeuron = length(tempPlotIndex);
    
    if isempty(p.Results.colorMod)
        colorMod = ones(1,nNeuron);
    else
        maxColorMod = prctile(abs(cell2mat(p.Results.colorMod)),95);
        colorMod = abs(p.Results.colorMod{k}) ./ maxColorMod;
        colorMod(colorMod>1) = 1;
    end

    for j = 1:nNeuron
        x = rois{j}(:,1); %freehand rois have the outlines in x-y coordinates
        y = rois{j}(:,2); %matlab matrices are inverted so x values are the 2nd column of mn coordinates, and y is the 1st columna

        if isempty(p.Results.plotFlag) || tempPlotFlag(j) 
            thisColor = C(colormapIndex(tempPlotIndex(j)),:);
            thisColor = thisColor + (1-thisColor) * (1-colorMod(j));
            %plot(x,y,'.','color',C(colormapIndex(peakIndex(cellIndex)),:),'MarkerSize',1);
            patch(h{k},x,y, thisColor, 'EdgeColor','none');
        else
            patch(h{k},x,y,[0.8 0.8 0.8],'EdgeColor','none');
            %plot(x,y,'.','color',[0.8 0.8 0.8],'MarkerSize',1);
        end       
        xmax = max([xmax;x]); ymax = max([ymax;y]); xmin = min([xmin;x]); ymin = min([ymin;y]); 
    end
    %set(h{k},'visible','off')
    if ~isempty (p.Results.title); title(h{k},p.Results.title{k}); end
    axis(h{k},[xmin xmax ymin ymax]);
    axis off;
    
end

end