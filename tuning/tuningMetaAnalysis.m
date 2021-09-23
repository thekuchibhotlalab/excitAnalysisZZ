dataPath = 'D:\labData\excitatory\tuning\masterData\';
%mouseGNG = { 'cd017', 'cd036', 'cd037', 'cd042','cd044'};
mouseGNG = { 'cd037'};
tStay =[];tSwitch =[];tAway =[];tLose = [];
fStay =[];fSwitch =[];fAway =[];fLose = [];
midF =[];midT =[];midStay =[];midLose = [];
tSideStay =[];tSideT =[];tSideF =[];tSideLose = [];
fSideStay =[];fSideT =[];fSideF =[];fSideLose = [];
noT = []; noF = [];

tUp =[];tDown =[];tStayD =[];tDownSw =[];
fUp =[];fDown =[];fStayD =[];fUpSw =[];
nUp =[];nDown =[];nStay =[];
nUpT =[];nUpF =[];nDownT =[];nDownF =[];

tUp_abs =[];tDown_abs =[];tStayD_abs =[];
fUp_abs =[];fDown_abs =[];fStayD_abs =[];
nUp_abs =[];nDown_abs =[];nStay_abs =[];
nUpT_abs =[];nUpF_abs =[];nDownT_abs =[];nDownF_abs =[];

nRespCell = 0;
nDecodeCell = 0;

for i = 1:length(mouseGNG)
    load([dataPath mouseGNG{i} '_prePostTuning.mat'],'peak','decoder');

    tStay = [tStay peak.tPre.stay];
    tSwitch = [tSwitch peak.tPre.switch];
    tAway = [tAway peak.tPre.away];
    tLose = [tLose peak.tPre.lose];

    fStay = [fStay peak.fPre.stay];
    fSwitch = [fSwitch peak.fPre.switch];
    fAway = [fAway peak.fPre.away];
    fLose = [fLose peak.fPre.lose];

    midF = [midF peak.middlePre.F];
    midT = [midT peak.middlePre.T];
    midStay =  [midStay peak.middlePre.stay];
    midLose = [midLose peak.middlePre.lose];

    tSideStay = [tSideStay peak.tSidePre.stay];
    tSideT = [tSideT peak.tSidePre.closeT];
    tSideF = [tSideF peak.tSidePre.closeF];
    tSideLose = [tSideLose peak.tSidePre.lose];

    fSideStay =[fSideStay peak.fSidePre.stay];
    fSideT = [fSideT peak.fSidePre.closeT];
    fSideF = [fSideF peak.fSidePre.closeF];
    fSideLose = [fSideLose peak.fSidePre.lose];
    
    noT = [noT peak.noPre.closeT];
    noF = [noF peak.noPre.closeF];

    tUp = [tUp decoder.tChange.up];
    tDown = [tDown decoder.tChange.down];
    tStayD = [tStayD decoder.tChange.stay];
    tDownSw = [tDownSw decoder.tChange.downSwitch];

    fUp = [fUp decoder.fChange.up];
    fDown = [fDown decoder.fChange.down];
    fStayD = [fStayD decoder.fChange.stay];
    fUpSw = [fUpSw decoder.fChange.upSwitch];

    nUp = [nUp decoder.nChange.up];
    nDown = [nDown decoder.nChange.down];
    nStay = [nStay decoder.nChange.stay];

    nUpT = [nUpT decoder.nChange.upT];
    nUpF = [nUpF decoder.nChange.upF];
    nDownT = [nDownT decoder.nChange.downT];
    nDownF = [nDownF decoder.nChange.downF];

    tUp_abs = [tUp_abs decoder.tChangeAbs.up];
    tDown_abs = [tDown_abs decoder.tChangeAbs.down];
    tStayD_abs = [tStayD_abs decoder.tChangeAbs.stay];

    fUp_abs = [fUp_abs decoder.fChangeAbs.up];
    fDown_abs = [fDown_abs decoder.fChangeAbs.down];
    fStayD_abs = [fStayD_abs decoder.fChangeAbs.stay];

    nUp_abs = [nUp_abs decoder.nChangeAbs.up];
    nDown_abs = [nDown_abs decoder.nChangeAbs.down];
    nStay_abs = [nStay_abs decoder.nChangeAbs.stay];

    nUpT_abs = [nUpT_abs decoder.nChangeAbs.upT];
    nUpF_abs = [nUpF_abs decoder.nChangeAbs.upF];
    
    nRespCell = nRespCell + peak.tuningUp + ...
        peak.tuningDown + peak.tuningStay;
    nDecodeCell = nDecodeCell + decoder.nDecode;
end


figure; 
subplot(2,3,1)
tPre = [sum(tSwitch) sum(tAway) sum(tStay) sum(tLose)];
tPre(tPre==0) = 0.1;
pie(tPre / sum(tPre),[1 1 1 1])
legend({'Switch2F', 'Away', 'Stay','Lose'},'Location','Best'); title ('Peak at Target')
subplot(2,3,2)
fPre = [sum(fSwitch) sum(fAway) sum(fStay) sum(fLose)];
fPre(fPre==0) = 0.1;
pie(fPre/sum(fPre),[1 1 1 1])
legend({'Switch2T', 'Away', 'Stay','Lose'},'Location','Best'); title ('Peak at Foil')
subplot(2,3,3)
midPre = [sum(midT) sum(midF) sum(midStay) sum(midLose)];
midPre(midPre==0) = 0.1;
pie(midPre / sum(midPre),[1 1 1 1])
legend({'move2T', 'move2F', 'Stay','Lose'},'Location','Best'); title ('Peak in Middle')
subplot(2,3,4)
tSidePre = [sum(tSideT) sum(tSideF) sum(tSideStay) sum(tSideLose)];
tSidePre(tSidePre==0)= 0.1;
pie(tSidePre / sum(tSidePre),[1 1 1 1])
legend({'Close2T', 'Close2F', 'Stay','Lose'},'Location','Best'); title ('Peak on Tside')
subplot(2,3,5)
fSidePre = [sum(fSideT) sum(fSideF) sum(fSideStay) sum(fSideLose)];
fSidePre(fSidePre==0) = 0.1;
pie(fSidePre/sum(fSidePre),[1 1 1 1])
legend({'Close2T', 'Close2F', 'Stay','Lose'},'Location','Best'); title ('Peak on Fside')
subplot(2,3,6)
noPre = [sum(noT) sum(noF)];
noPre(noPre==0) = 0.1;
pie(noPre/sum(noPre),[1 1 ])
legend({'Close2T', 'Close2F'},'Location','Best'); title ('Gain Selectivity')


figure; 
subplot(1,3,1)
tDecode = [sum(tUp) sum(tDown-tDownSw) sum(tDownSw) sum(tStayD)];
tDecode(tDecode==0) = 0.1;
pie(tDecode / sum(tDecode),[1 1 1 1])
legend({'Incease', 'Decrease&No-Switch','Decrease&Switch', 'Stay'},'Location','Best'); title ('Target Preferring')
subplot(1,3,2)
fDecode = [sum(fDown) sum(fUp-fUpSw) sum(fUpSw) sum(fStayD)];
fDecode(fDecode==0) = 0.1;
pie(fDecode/sum(fDecode),[1 1 1 1])
legend({'Incease', 'Decrease&No-Switch','Decrease&Switch', 'Stay'}); title ('Foil Preferring')
subplot(1,3,3)
nDecode = [ sum(nUpT) sum(nUpF) sum(nDownT) sum(nDownF)  sum(nStay) ];
nDecode (nDecode == 0) = 0.1;
pie(nDecode / sum(nDecode),[1 1 1 1 1])
legend({'Move-To-T&End-State-T', 'Move-To-T&End-State-F','Move-To-F&End-State-T',...
    'Move-To-F&End-State-F','Stay'},'Location','Best'); title ('No Preference')


figure; 
subplot(1,3,1)
tDecode = [sum(tUp_abs) sum(tDown_abs) sum(tStayD_abs)];
tDecode(tDecode==0) = 0.1;
pie(tDecode / sum(tDecode),[1 1 1])
legend({'Incease', 'Decrease', 'Stay'},'Location','Best'); title ('Target Preferring')
subplot(1,3,2)
fDecode = [sum(fDown_abs) sum(fUp_abs) sum(fStayD_abs)];
fDecode(fDecode==0) = 0.1;
pie(fDecode/sum(fDecode),[1 1 1])
legend({'Incease', 'Decrease', 'Stay'}); title ('Foil Preferring')
subplot(1,3,3)
nDecode = [ sum(nUpT_abs) sum(nUpF_abs)  sum(nStay_abs) ];
nDecode (nDecode == 0) = 0.1;
pie(nDecode / sum(nDecode),[1 1 1])
legend({'Gain selectivity to T', 'Gain selectivity to F','Stay'},'Location','Best'); title ('No Preference')

    