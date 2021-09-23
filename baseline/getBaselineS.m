function getBaselineS(mouse)
    %mouse = 'cd041';
    datapath = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' mouse];
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
    tic;
    %figure; scatter(mean(s{6},2),mean(c{6},2));
    configPath = 'C:\Users\zzhu34\Documents\tempdata\excitData\config\mouse';
    filename = [configPath '\' mouse '_config.csv'];
    configTable = readtable(filename);
    baselineStemp = s(strcmp(configTable.BehavType,'Baseline'));
    baselineS = cell2mat(baselineStemp);
    baselineS = reshape(baselineS,size(baselineStemp{1},1),size(baselineStemp{1},2),length(baselineStemp));
    baselineDay = 2:(size(ishere,2)-2);
    ishereS = ishere(2:end-2);
    save([datapath '\corrected_plane0\' mouse '_baselineS.mat'],'baselineS','ishereS');
    t = toc;
    disp([mouse ' is finished. Time = ' num2str(t,'%.2f') ' sec'])
end