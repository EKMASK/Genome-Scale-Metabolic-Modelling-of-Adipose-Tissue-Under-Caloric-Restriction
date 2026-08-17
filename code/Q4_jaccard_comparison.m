clear; clc; close all;

%% ---- Load models and metadata ----
setup_models;

%% ---- 1. Reaction and metabolite counts per individual ----
figure('Name','Model Sizes','Position',[100 100 900 400]);

% Reactions bar chart
subplot(1,2,1); hold on;
for i = 1:nParticipants
    if idxVLCD(i)
        bar(i, nRxns(i), 'FaceColor', [0.85 0.33 0.1]);
    else
        bar(i, nRxns(i), 'FaceColor', [0.0 0.45 0.74]);
    end
end
set(gca, 'XTick', 1:nParticipants, 'XTickLabel', participantIDs, 'XTickLabelRotation', 45);
ylabel('Number of reactions');
title('Reactions per baseline model');
% Dummy plots for correct legend
h1 = bar(NaN, NaN, 'FaceColor', [0.85 0.33 0.1]);
h2 = bar(NaN, NaN, 'FaceColor', [0.0 0.45 0.74]);
legend([h1 h2], {'VLCD','LCD'}, 'Location','best');
hold off;

% Metabolites bar chart
subplot(1,2,2); hold on;
for i = 1:nParticipants
    if idxVLCD(i)
        bar(i, nMets(i), 'FaceColor', [0.85 0.33 0.1]);
    else
        bar(i, nMets(i), 'FaceColor', [0.0 0.45 0.74]);
    end
end
set(gca, 'XTick', 1:nParticipants, 'XTickLabel', participantIDs, 'XTickLabelRotation', 45);
ylabel('Number of metabolites');
title('Metabolites per baseline model');
hold off;

sgtitle(sprintf('Baseline model sizes (%s)', reconstructionType));
saveas(gcf, sprintf('Q4_model_sizes_%s.png', reconstructionType));

% Print table for report
fprintf('\n%-6s %-5s %-4s %8s %8s\n', 'ID', 'Diet', 'Sex', 'Rxns', 'Mets');
fprintf('%s\n', repmat('-', 1, 38));
for i = 1:nParticipants
    fprintf('%-6s %-5s %-4s %8d %8d\n', ...
        participantIDs{i}, dietGroup{i}, sex{i}, nRxns(i), nMets(i));
end

%% ---- 2. Compute pairwise Jaccard similarity ----
% Jaccard(A,B) = |A ? B| / |A ? B| using reaction ID sets.

allRxnSets = cell(1, nParticipants);
for i = 1:nParticipants
    allRxnSets{i} = models{i}.rxns;
end

allRxnsUnion = unique(vertcat(allRxnSets{:}));
fprintf('\nTotal unique reactions across all 20 models: %d\n', length(allRxnsUnion));

jaccardMatrix = zeros(nParticipants);
for i = 1:nParticipants
    for j = i:nParticipants
        nInt   = length(intersect(allRxnSets{i}, allRxnSets{j}));
        nUnion = length(union(allRxnSets{i}, allRxnSets{j}));
        jaccardMatrix(i,j) = nInt / nUnion;
        jaccardMatrix(j,i) = jaccardMatrix(i,j);
    end
end

jaccardDist = 1 - jaccardMatrix;

%% ---- 3. Jaccard heatmap (original order) ----
figure('Name','Jaccard Heatmap','Position',[100 100 800 700]);
imagesc(jaccardMatrix);
colormap(parula); colorbar;
caxis([min(jaccardMatrix(:)) 1]);
set(gca, 'XTick', 1:nParticipants, 'XTickLabel', participantIDs, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:nParticipants, 'YTickLabel', participantIDs);
title(sprintf('Jaccard Similarity — Baseline (%s)', reconstructionType));
saveas(gcf, sprintf('Q4_jaccard_heatmap_%s.png', reconstructionType));

%% ---- 4. Hierarchical clustering + dendrogram ----
distVec = squareform(jaccardDist);
linkTree = linkage(distVec, 'average');

figure('Name','Dendrogram','Position',[100 100 800 400]);
[~, ~, outperm] = dendrogram(linkTree, 0, 'Labels', participantIDs, ...
    'ColorThreshold', 'default');
title(sprintf('Hierarchical Clustering (Jaccard distance, average linkage) — %s', reconstructionType));
ylabel('Jaccard distance');
set(gca, 'XTickLabelRotation', 45);
saveas(gcf, sprintf('Q4_dendrogram_%s.png', reconstructionType));

%% ---- 5. Clustered heatmap with diet annotation ----
orderedJaccard = jaccardMatrix(outperm, outperm);
orderedIDs     = participantIDs(outperm);
orderedDiet    = dietGroup(outperm);

figure('Name','Clustered Jaccard','Position',[100 100 800 700]);
imagesc(orderedJaccard);
colormap(parula);
cb = colorbar; cb.Label.String = 'Jaccard Similarity';
set(gca, 'XTick', 1:nParticipants, 'XTickLabel', orderedIDs, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:nParticipants, 'YTickLabel', orderedIDs);
title(sprintf('Clustered Jaccard Similarity — %s', reconstructionType));

% Diet annotation bar along top
hold on;
for i = 1:nParticipants
    if strcmp(orderedDiet{i}, 'VLCD')
        rectangle('Position', [i-0.5, 0.0, 1, 0.5], 'FaceColor', [0.85 0.33 0.1], 'EdgeColor','none');
    else
        rectangle('Position', [i-0.5, 0.0, 1, 0.5], 'FaceColor', [0.0 0.45 0.74], 'EdgeColor','none');
    end
end
hold off;
saveas(gcf, sprintf('Q4_clustered_jaccard_%s.png', reconstructionType));

%% ---- 6. Core reactions per diet group ----
rxnsVLCD_core = allRxnSets{find(idxVLCD, 1)};
for i = find(idxVLCD)
    rxnsVLCD_core = intersect(rxnsVLCD_core, allRxnSets{i});
end

rxnsLCD_core = allRxnSets{find(idxLCD, 1)};
for i = find(idxLCD)
    rxnsLCD_core = intersect(rxnsLCD_core, allRxnSets{i});
end

rxnsOnlyVLCD = setdiff(rxnsVLCD_core, rxnsLCD_core);
rxnsOnlyLCD  = setdiff(rxnsLCD_core, rxnsVLCD_core);
rxnsShared   = intersect(rxnsVLCD_core, rxnsLCD_core);

fprintf('\n--- Core reactions per diet group ---\n');
fprintf('VLCD core: %d | LCD core: %d | Shared: %d\n', ...
    length(rxnsVLCD_core), length(rxnsLCD_core), length(rxnsShared));
fprintf('VLCD-only core: %d | LCD-only core: %d\n', ...
    length(rxnsOnlyVLCD), length(rxnsOnlyLCD));

%% ---- 7. MDS coloured by diet, sex, BMI ----
[Y, eigvals] = cmdscale(jaccardDist);
explVar = eigvals / sum(abs(eigvals)) * 100;

figure('Name','MDS Jaccard','Position',[100 100 1200 400]);

% By diet
subplot(1,3,1); hold on;
hV = []; hL = [];
for i = 1:nParticipants
    if idxVLCD(i)
        hV = scatter(Y(i,1), Y(i,2), 80, [0.85 0.33 0.1], 'filled');
    else
        hL = scatter(Y(i,1), Y(i,2), 80, [0.0 0.45 0.74], 'filled');
    end
    text(Y(i,1)+0.003, Y(i,2)+0.003, participantIDs{i}, 'FontSize', 7);
end
xlabel(sprintf('MDS1 (%.1f%%)', explVar(1)));
ylabel(sprintf('MDS2 (%.1f%%)', explVar(2)));
title('By diet group');
legend([hV hL], {'VLCD','LCD'}, 'Location','best');
hold off;

% By sex
subplot(1,3,2); hold on;
hM = []; hF = [];
for i = 1:nParticipants
    if strcmp(sex{i}, 'M')
        hM = scatter(Y(i,1), Y(i,2), 80, [0.47 0.67 0.19], 'filled');
    else
        hF = scatter(Y(i,1), Y(i,2), 80, [0.64 0.08 0.18], 'filled');
    end
    text(Y(i,1)+0.003, Y(i,2)+0.003, participantIDs{i}, 'FontSize', 7);
end
xlabel(sprintf('MDS1 (%.1f%%)', explVar(1)));
ylabel(sprintf('MDS2 (%.1f%%)', explVar(2)));
title('By sex');
legend([hM hF], {'Male','Female'}, 'Location','best');
hold off;

% By BMI
subplot(1,3,3);
scatter(Y(:,1), Y(:,2), 80, BMI_baseline, 'filled');
colormap(gca, hot);
cb2 = colorbar; cb2.Label.String = 'BMI (baseline)';
for i = 1:nParticipants
    text(Y(i,1)+0.003, Y(i,2)+0.003, participantIDs{i}, 'FontSize', 7);
end
xlabel(sprintf('MDS1 (%.1f%%)', explVar(1)));
ylabel(sprintf('MDS2 (%.1f%%)', explVar(2)));
title('By BMI');

sgtitle(sprintf('MDS on Jaccard distance — %s', reconstructionType));
saveas(gcf, sprintf('Q4_MDS_jaccard_%s.png', reconstructionType));

%% ---- 8. Within-group vs between-group Jaccard ----
withinVLCD = []; withinLCD = []; betweenGroups = [];
vlcdIdx = find(idxVLCD); lcdIdx = find(idxLCD);

for i = 1:length(vlcdIdx)
    for j = i+1:length(vlcdIdx)
        withinVLCD(end+1) = jaccardMatrix(vlcdIdx(i), vlcdIdx(j)); %#ok
    end
end
for i = 1:length(lcdIdx)
    for j = i+1:length(lcdIdx)
        withinLCD(end+1) = jaccardMatrix(lcdIdx(i), lcdIdx(j)); %#ok
    end
end
for i = 1:length(vlcdIdx)
    for j = 1:length(lcdIdx)
        betweenGroups(end+1) = jaccardMatrix(vlcdIdx(i), lcdIdx(j)); %#ok
    end
end

fprintf('Mean Jaccard:    %.4f +/- %.4f\n', ...
    mean(jaccardMatrix(triu(true(nParticipants),1))), ...
    std(jaccardMatrix(triu(true(nParticipants),1))));
fprintf('Within VLCD:     %.4f +/- %.4f\n', mean(withinVLCD), std(withinVLCD));
fprintf('Within LCD:      %.4f +/- %.4f\n', mean(withinLCD), std(withinLCD));
fprintf('Between groups:  %.4f +/- %.4f\n', mean(betweenGroups), std(betweenGroups));
