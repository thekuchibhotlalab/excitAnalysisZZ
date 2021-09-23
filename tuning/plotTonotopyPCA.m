clear;

allMouse = 'cd037';
savePath = 'D:\labData\excitatory\tuning\masterData\';
load([savePath '\' allMouse '\' 'prePostTuning.mat']);

%%
pre = squeeze(mean(peakActPre,2));
post = squeeze(mean(peakActPost,2));
C = jet;
colormapIndex = round(linspace(1,64,17));
allColorPre = zeros(nNeuron,3);
for j = 1:nNeuron
    if responsiveFlagPre(j)
        thisColor = C(colormapIndex(preTuning.tuningPeak(j)),:);
    else; thisColor = [0.8 0.8 0.8];
    end
    allColorPre(j,:) = thisColor;
end

allColorPost = zeros(nNeuron,3);
for j = 1:nNeuron
    if responsiveFlagPost(j)
        thisColor = C(colormapIndex(postTuning.tuningPeak(j)),:);
    else; thisColor = [0.8 0.8 0.8];
    end
    allColorPost(j,:) = thisColor;
end
[basis, varExp, proj, covMat] = func_pca([pre;post]);

figure; plot(cumsum(varExp))


figure; subplot(1,3,1); scatter(proj(1,:),proj(2,:),20, allColorPre,'filled'); 
xlabel('pc1'); ylabel('pc2')
subplot(1,3,2); scatter(proj(1,:),proj(3,:),20, allColorPre,'filled')
xlabel('pc1'); ylabel('pc3')
subplot(1,3,3); scatter(proj(2,:),proj(3,:),20, allColorPre,'filled')
xlabel('pc2'); ylabel('pc3')
suptitle('pre learning')

figure; subplot(1,3,1); scatter(proj(1,:),proj(2,:),20, allColorPost,'filled');
xlabel('pc1'); ylabel('pc2')
subplot(1,3,2); scatter(proj(1,:),proj(3,:),20, allColorPost,'filled')
xlabel('pc1'); ylabel('pc3')
subplot(1,3,3); scatter(proj(2,:),proj(3,:),20, allColorPost,'filled')
xlabel('pc2'); ylabel('pc3')
suptitle('post learning')

figure;subplot(1,3,1);plot(basis(1:17,1)); hold on; plot(basis(18:end,1));
xlabel('tones'); ylabel('weight'); title('pc1'); xlim([1 17])
subplot(1,3,2);plot(basis(1:17,2)); hold on; plot(basis(18:end,2));
xlabel('tones'); ylabel('weight'); title('pc2'); xlim([1 17])
subplot(1,3,3);plot(basis(1:17,3)); hold on; plot(basis(18:end,3));
xlabel('tones'); ylabel('weight'); title('pc3'); xlim([1 17])