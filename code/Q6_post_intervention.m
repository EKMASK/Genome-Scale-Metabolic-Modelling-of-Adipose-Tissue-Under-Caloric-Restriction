

clear; clc; close all;

changeCobraSolver('gurobi', 'all');

%% ---- Load BASELINE models ----
setup_models;  % loads baseline into 'models'
baselineModels = models;
baselineRxns   = nRxns;
baselineMets   = nMets;

%% ---- Load POST-INTERVENTION models ----
postDir = 'iMAT_Foguet_post_models';
fprintf('\nLoading post-intervention models from: %s\n', postDir);

postModels = cell(1, nParticipants);
postRxns   = zeros(1, nParticipants);
postMets   = zeros(1, nParticipants);

% Filenames have typos matching original data columns
% P01: post_intervetion, P25: post_intervnetion, rest: post_intervention
for i = 1:nParticipants
    % Try the standard naming first
    fname = fullfile(postDir, sprintf('model_%s_post_intervention_%s.mat', ...
        participantIDs{i}, dietGroup{i}));
    
    % If not found, try known typo variants
    if ~isfile(fname)
        fname = fullfile(postDir, sprintf('model_%s_post_intervetion_%s.mat', ...
            participantIDs{i}, dietGroup{i}));
    end
    if ~isfile(fname)
        fname = fullfile(postDir, sprintf('model_%s_post_intervnetion_%s.mat', ...
            participantIDs{i}, dietGroup{i}));
    end
    
    if ~isfile(fname)
        warning('Post model not found for %s', participantIDs{i});
        continue;
    end
    
    tmp = load(fname);
    if isfield(tmp, 'tissueModel')
        postModels{i} = tmp.tissueModel;
    elseif isfield(tmp, 'model')
        postModels{i} = tmp.model;
    else
        fn = fieldnames(tmp);
        postModels{i} = tmp.(fn{1});
    end
    
    postRxns(i) = length(postModels{i}.rxns);
    postMets(i) = length(postModels{i}.mets);
    fprintf('  %s (%s): %d rxns, %d mets\n', participantIDs{i}, dietGroup{i}, postRxns(i), postMets(i));
end

nPostLoaded = sum(~cellfun(@isempty, postModels));
fprintf('Loaded %d/%d post-intervention models.\n', nPostLoaded, nParticipants);

%% ---- 1. Model size comparison: baseline vs post ----
figure('Name','Model Size Comparison','Position',[100 100 1000 400]);

subplot(1,2,1); hold on;
x = 1:nParticipants;
bar(x - 0.15, baselineRxns, 0.3, 'FaceColor', [0.3 0.6 0.9]);
bar(x + 0.15, postRxns, 0.3, 'FaceColor', [0.9 0.4 0.3]);
set(gca, 'XTick', x, 'XTickLabel', participantIDs, 'XTickLabelRotation', 45);
ylabel('Number of reactions');
title('Reactions: baseline vs post');
legend({'Baseline','Post-intervention'}, 'Location','best');
hold off;

subplot(1,2,2); hold on;
bar(x - 0.15, baselineMets, 0.3, 'FaceColor', [0.3 0.6 0.9]);
bar(x + 0.15, postMets, 0.3, 'FaceColor', [0.9 0.4 0.3]);
set(gca, 'XTick', x, 'XTickLabel', participantIDs, 'XTickLabelRotation', 45);
ylabel('Number of metabolites');
title('Metabolites: baseline vs post');
legend({'Baseline','Post-intervention'}, 'Location','best');
hold off;

sgtitle('Model sizes: baseline vs post-intervention');
saveas(gcf, 'Q6_model_size_comparison.png');

% Print table
fprintf('\n%-6s %-5s %8s %8s %8s %8s\n', 'ID', 'Diet', 'Rxn_BL', 'Rxn_Post', 'Met_BL', 'Met_Post');
fprintf('%s\n', repmat('-', 1, 52));
for i = 1:nParticipants
    fprintf('%-6s %-5s %8d %8d %8d %8d\n', participantIDs{i}, dietGroup{i}, ...
        baselineRxns(i), postRxns(i), baselineMets(i), postMets(i));
end

%% ---- 2. Jaccard: baseline vs post per participant ----
jaccardPerParticipant = zeros(1, nParticipants);
for i = 1:nParticipants
    if isempty(postModels{i}), continue; end
    nInt   = length(intersect(baselineModels{i}.rxns, postModels{i}.rxns));
    nUnion = length(union(baselineModels{i}.rxns, postModels{i}.rxns));
    jaccardPerParticipant(i) = nInt / nUnion;
end

figure('Name','Jaccard BL vs Post','Position',[100 100 700 350]);
hold on;
for i = 1:nParticipants
    if idxVLCD(i)
        bar(i, jaccardPerParticipant(i), 'FaceColor', [0.85 0.33 0.1]);
    else
        bar(i, jaccardPerParticipant(i), 'FaceColor', [0.0 0.45 0.74]);
    end
end
set(gca, 'XTick', 1:nParticipants, 'XTickLabel', participantIDs, 'XTickLabelRotation', 45);
ylabel('Jaccard similarity');
title('Jaccard similarity: baseline vs post-intervention (per participant)');
ylim([0 1]);
hV = bar(NaN, NaN, 'FaceColor', [0.85 0.33 0.1]);
hL = bar(NaN, NaN, 'FaceColor', [0.0 0.45 0.74]);
legend([hV hL], {'VLCD','LCD'}, 'Location','best');
hold off;
saveas(gcf, 'Q6_jaccard_BL_vs_post.png');

fprintf('\nJaccard (baseline vs post) per participant:\n');
for i = 1:nParticipants
    fprintf('  %s (%s): %.4f\n', participantIDs{i}, dietGroup{i}, jaccardPerParticipant(i));
end
fprintf('VLCD mean: %.4f +/- %.4f\n', mean(jaccardPerParticipant(idxVLCD)), std(jaccardPerParticipant(idxVLCD)));
fprintf('LCD mean:  %.4f +/- %.4f\n', mean(jaccardPerParticipant(idxLCD)), std(jaccardPerParticipant(idxLCD)));

%% ---- 3. Reactions unique to post-intervention per diet group ----
% For each participant: which reactions are NEW in post (not in baseline)?
% And which disappeared (in baseline, not in post)?

newInPostVLCD = {};  lostInPostVLCD = {};
newInPostLCD  = {};  lostInPostLCD  = {};

for i = 1:nParticipants
    if isempty(postModels{i}), continue; end
    newRxns  = setdiff(postModels{i}.rxns, baselineModels{i}.rxns);
    lostRxns = setdiff(baselineModels{i}.rxns, postModels{i}.rxns);
    
    if idxVLCD(i)
        newInPostVLCD  = [newInPostVLCD; newRxns];
        lostInPostVLCD = [lostInPostVLCD; lostRxns];
    else
        newInPostLCD  = [newInPostLCD; newRxns];
        lostInPostLCD = [lostInPostLCD; lostRxns];
    end
end

% Count frequency of each new/lost reaction
[newVLCD_unique, ~, ic] = unique(newInPostVLCD);
newVLCD_counts = accumarray(ic, 1);
[lostVLCD_unique, ~, ic] = unique(lostInPostVLCD);
lostVLCD_counts = accumarray(ic, 1);

[newLCD_unique, ~, ic] = unique(newInPostLCD);
newLCD_counts = accumarray(ic, 1);
[lostLCD_unique, ~, ic] = unique(lostInPostLCD);
lostLCD_counts = accumarray(ic, 1);

fprintf('\n--- Reactions gained/lost post-intervention ---\n');
fprintf('VLCD: %d unique new reactions, %d unique lost reactions\n', ...
    length(newVLCD_unique), length(lostVLCD_unique));
fprintf('LCD:  %d unique new reactions, %d unique lost reactions\n', ...
    length(newLCD_unique), length(lostLCD_unique));

% Reactions gained in BOTH diet groups
newInBoth = intersect(newVLCD_unique, newLCD_unique);
newOnlyVLCD = setdiff(newVLCD_unique, newLCD_unique);
newOnlyLCD  = setdiff(newLCD_unique, newVLCD_unique);

fprintf('\nNew reactions shared by both diets: %d\n', length(newInBoth));
fprintf('New reactions VLCD-only: %d\n', length(newOnlyVLCD));
fprintf('New reactions LCD-only:  %d\n', length(newOnlyLCD));

% Print top new reactions (most frequent across participants)
fprintf('\nTop new reactions in VLCD post-intervention:\n');
[~, sortIdx] = sort(newVLCD_counts, 'descend');
for k = 1:min(15, length(sortIdx))
    fprintf('  %-40s (in %d/10 participants)\n', newVLCD_unique{sortIdx(k)}, newVLCD_counts(sortIdx(k)));
end

fprintf('\nTop new reactions in LCD post-intervention:\n');
[~, sortIdx] = sort(newLCD_counts, 'descend');
for k = 1:min(15, length(sortIdx))
    fprintf('  %-40s (in %d/10 participants)\n', newLCD_unique{sortIdx(k)}, newLCD_counts(sortIdx(k)));
end

%% ---- 4. Jaccard heatmap: all 40 models together ----
% Combine baseline + post models and compute full Jaccard matrix
allModels40 = [baselineModels, postModels];
labels40    = cell(1, 40);
group40     = cell(1, 40);  % 'BL_VLCD', 'BL_LCD', 'Post_VLCD', 'Post_LCD'

for i = 1:nParticipants
    labels40{i}    = [participantIDs{i} '_BL'];
    labels40{i+20} = [participantIDs{i} '_Post'];
    if idxVLCD(i)
        group40{i}    = 'BL_VLCD';
        group40{i+20} = 'Post_VLCD';
    else
        group40{i}    = 'BL_LCD';
        group40{i+20} = 'Post_LCD';
    end
end

jaccard40 = zeros(40);
for i = 1:40
    for j = i:40
        if isempty(allModels40{i}) || isempty(allModels40{j})
            jaccard40(i,j) = NaN;
        else
            nInt   = length(intersect(allModels40{i}.rxns, allModels40{j}.rxns));
            nUnion = length(union(allModels40{i}.rxns, allModels40{j}.rxns));
            jaccard40(i,j) = nInt / nUnion;
        end
        jaccard40(j,i) = jaccard40(i,j);
    end
end

% Replace NaN with 0 for clustering
jaccard40_clean = jaccard40;
jaccard40_clean(isnan(jaccard40_clean)) = 0;

dist40 = 1 - jaccard40_clean;
dist40 = (dist40 + dist40') / 2;
dist40(1:41:end) = 0;
distVec40 = squareform(dist40);
linkTree40 = linkage(distVec40, 'average');

figure('Name','Jaccard 40 models','Position',[100 100 1000 500]);
[~, ~, outperm40] = dendrogram(linkTree40, 0, 'Labels', labels40, 'ColorThreshold', 'default');
title('Hierarchical Clustering â€” All 40 models (baseline + post)');
ylabel('Jaccard distance');
set(gca, 'XTickLabelRotation', 90, 'FontSize', 6);
saveas(gcf, 'Q6_dendrogram_40models.png');

% Clustered heatmap
ord40 = jaccard40_clean(outperm40, outperm40);
ordLabels40 = labels40(outperm40);
ordGroup40  = group40(outperm40);

figure('Name','Clustered Jaccard 40','Position',[100 100 1000 900]);
imagesc(ord40); colormap(parula);
cb = colorbar; cb.Label.String = 'Jaccard Similarity';
set(gca, 'XTick', 1:40, 'XTickLabel', ordLabels40, 'XTickLabelRotation', 90, 'FontSize', 6);
set(gca, 'YTick', 1:40, 'YTickLabel', ordLabels40, 'FontSize', 6);
title('Clustered Jaccard â€” Baseline + Post-intervention');

% Annotation bar: colour by group
hold on;
groupColors = containers.Map({'BL_VLCD','BL_LCD','Post_VLCD','Post_LCD'}, ...
    {[0.85 0.33 0.1], [0.0 0.45 0.74], [1.0 0.6 0.3], [0.4 0.7 1.0]});
for i = 1:40
    clr = groupColors(ordGroup40{i});
    rectangle('Position', [i-0.5, 0.0, 1, 0.8], 'FaceColor', clr, 'EdgeColor','none');
end
hold off;
saveas(gcf, 'Q6_clustered_jaccard_40.png');

%% ---- 5. FBA on post-intervention models ----
fprintf('\nRunning FBA on post-intervention models...\n');

% Set objective
for i = 1:nParticipants
    if isempty(postModels{i}), continue; end
    model = postModels{i};
    
    patterns = {'DM_atp', 'ATPM', 'ATPhyd', 'ATPase'};
    for k = 1:length(patterns)
        matches = find(contains(model.rxns, patterns{k}, 'IgnoreCase', true));
        if ~isempty(matches)
            model = changeObjective(model, model.rxns{matches(1)});
            break;
        end
    end
    postModels{i} = model;
end

% Build union of ALL reactions (baseline + post)
allRxnsAll = {};
for i = 1:nParticipants
    allRxnsAll = union(allRxnsAll, baselineModels{i}.rxns);
    if ~isempty(postModels{i})
        allRxnsAll = union(allRxnsAll, postModels{i}.rxns);
    end
end
nAllRxns = length(allRxnsAll);

% FBA on baseline (re-run to get aligned flux vectors)
fprintf('\nFBA on baseline models...\n');
fluxBL = zeros(nAllRxns, nParticipants);
objBL  = zeros(1, nParticipants);
feasBL = true(1, nParticipants);

for i = 1:nParticipants
    model = baselineModels{i};
    % Set objective
    patterns = {'DM_atp', 'ATPM', 'ATPhyd', 'ATPase'};
    for k = 1:length(patterns)
        matches = find(contains(model.rxns, patterns{k}, 'IgnoreCase', true));
        if ~isempty(matches)
            model = changeObjective(model, model.rxns{matches(1)});
            break;
        end
    end
    
    sol = optimizeCbModel(model, 'max');
    if sol.stat == 1
        objBL(i) = sol.f;
        [~, ia, ib] = intersect(allRxnsAll, model.rxns);
        fluxBL(ia, i) = sol.v(ib);
        fprintf('  %s BL: obj=%.4f\n', participantIDs{i}, sol.f);
    else
        feasBL(i) = false;
        fprintf('  %s BL: INFEASIBLE\n', participantIDs{i});
    end
end

% FBA on post
fprintf('\nFBA on post-intervention models...\n');
fluxPost = zeros(nAllRxns, nParticipants);
objPost  = zeros(1, nParticipants);
feasPost = true(1, nParticipants);

for i = 1:nParticipants
    if isempty(postModels{i})
        feasPost(i) = false;
        continue;
    end
    
    sol = optimizeCbModel(postModels{i}, 'max');
    if sol.stat == 1
        objPost(i) = sol.f;
        [~, ia, ib] = intersect(allRxnsAll, postModels{i}.rxns);
        fluxPost(ia, i) = sol.v(ib);
        fprintf('  %s Post: obj=%.4f\n', participantIDs{i}, sol.f);
    else
        feasPost(i) = false;
        fprintf('  %s Post: INFEASIBLE\n', participantIDs{i});
    end
end

%% ---- 6. Compare objective values ----
feasBoth = feasBL & feasPost;
fprintf('\nBoth feasible: %d/%d\n', sum(feasBoth), nParticipants);

figure('Name','Obj Values BL vs Post','Position',[100 100 800 400]);
hold on;
bIdx = find(feasBoth);
for i = 1:length(bIdx)
    k = bIdx(i);
    if idxVLCD(k)
        bar(i - 0.15, objBL(k), 0.3, 'FaceColor', [0.85 0.33 0.1]);
        bar(i + 0.15, objPost(k), 0.3, 'FaceColor', [1.0 0.6 0.3]);
    else
        bar(i - 0.15, objBL(k), 0.3, 'FaceColor', [0.0 0.45 0.74]);
        bar(i + 0.15, objPost(k), 0.3, 'FaceColor', [0.4 0.7 1.0]);
    end
end
set(gca, 'XTick', 1:length(bIdx), 'XTickLabel', participantIDs(bIdx), 'XTickLabelRotation', 45);
ylabel('ATPM flux');
title('FBA Objective: Baseline vs Post-intervention');
legend({'BL (VLCD)','Post (VLCD)','BL (LCD)','Post (LCD)'}, 'Location','best');
hold off;
saveas(gcf, 'Q6_objective_BL_vs_post.png');

%% ---- 7. Cosine similarity: baseline vs post per participant ----
cosPerParticipant = zeros(1, nParticipants);
for i = 1:nParticipants
    if ~feasBoth(i), continue; end
    a = fluxBL(:,i); b = fluxPost(:,i);
    d = norm(a) * norm(b);
    if d > 0
        cosPerParticipant(i) = dot(a,b) / d;
    end
end

figure('Name','Cosine BL vs Post','Position',[100 100 700 350]);
hold on;
for i = 1:nParticipants
    if ~feasBoth(i), continue; end
    if idxVLCD(i)
        bar(i, cosPerParticipant(i), 'FaceColor', [0.85 0.33 0.1]);
    else
        bar(i, cosPerParticipant(i), 'FaceColor', [0.0 0.45 0.74]);
    end
end
set(gca, 'XTick', 1:nParticipants, 'XTickLabel', participantIDs, 'XTickLabelRotation', 45);
ylabel('Cosine similarity');
title('Flux similarity: baseline vs post (per participant)');
ylim([0 1]);
hold off;
saveas(gcf, 'Q6_cosine_BL_vs_post.png');

fprintf('\nCosine (BL vs Post) per participant:\n');
for i = 1:nParticipants
    if feasBoth(i)
        fprintf('  %s (%s): %.4f\n', participantIDs{i}, dietGroup{i}, cosPerParticipant(i));
    end
end
fprintf('VLCD mean: %.4f +/- %.4f\n', mean(cosPerParticipant(idxVLCD & feasBoth)), std(cosPerParticipant(idxVLCD & feasBoth)));
fprintf('LCD mean:  %.4f +/- %.4f\n', mean(cosPerParticipant(idxLCD & feasBoth)), std(cosPerParticipant(idxLCD & feasBoth)));

%% ---- 8. Significantly different fluxes per diet group ----
% For each diet group, paired comparison: baseline vs post
% Use Wilcoxon signed-rank test (non-parametric, paired, small n)

fprintf('\n--- Significantly different fluxes (Wilcoxon signed-rank, p<0.05) ---\n');

% VLCD
vlcdIdx_feas = find(idxVLCD & feasBoth);
fluxDiffVLCD = fluxPost(:, vlcdIdx_feas) - fluxBL(:, vlcdIdx_feas);
pVLCD = ones(nAllRxns, 1);
for r = 1:nAllRxns
    vals = fluxDiffVLCD(r, :);
    if any(vals ~= 0)
        pVLCD(r) = signrank(fluxBL(r, vlcdIdx_feas), fluxPost(r, vlcdIdx_feas));
    end
end
sigVLCD = find(pVLCD < 0.05);
fprintf('VLCD: %d reactions with p < 0.05\n', length(sigVLCD));

% LCD
lcdIdx_feas = find(idxLCD & feasBoth);
fluxDiffLCD = fluxPost(:, lcdIdx_feas) - fluxBL(:, lcdIdx_feas);
pLCD = ones(nAllRxns, 1);
for r = 1:nAllRxns
    vals = fluxDiffLCD(r, :);
    if any(vals ~= 0)
        pLCD(r) = signrank(fluxBL(r, lcdIdx_feas), fluxPost(r, lcdIdx_feas));
    end
end
sigLCD = find(pLCD < 0.05);
fprintf('LCD:  %d reactions with p < 0.05\n', length(sigLCD));

% Overlap
sigBoth = intersect(sigVLCD, sigLCD);
sigOnlyVLCD = setdiff(sigVLCD, sigLCD);
sigOnlyLCD  = setdiff(sigLCD, sigVLCD);

fprintf('\nSignificant in both diets: %d\n', length(sigBoth));
fprintf('Significant VLCD-only: %d\n', length(sigOnlyVLCD));
fprintf('Significant LCD-only:  %d\n', length(sigOnlyLCD));

% Print top significantly different reactions with mean flux change
fprintf('\nTop significant reactions (VLCD) â€” sorted by |mean flux change|:\n');
meanDiffVLCD = mean(fluxDiffVLCD, 2);
[~, sortIdx] = sort(abs(meanDiffVLCD(sigVLCD)), 'descend');
for k = 1:min(20, length(sortIdx))
    r = sigVLCD(sortIdx(k));
    fprintf('  %-45s  mean_diff=%+.4f  p=%.4f\n', allRxnsAll{r}, meanDiffVLCD(r), pVLCD(r));
end

fprintf('\nTop significant reactions (LCD) â€” sorted by |mean flux change|:\n');
meanDiffLCD = mean(fluxDiffLCD, 2);
[~, sortIdx] = sort(abs(meanDiffLCD(sigLCD)), 'descend');
for k = 1:min(20, length(sortIdx))
    r = sigLCD(sortIdx(k));
    fprintf('  %-45s  mean_diff=%+.4f  p=%.4f\n', allRxnsAll{r}, meanDiffLCD(r), pLCD(r));
end

% Shared significant reactions
if ~isempty(sigBoth)
    fprintf('\nReactions significant in BOTH diets:\n');
    for k = 1:min(20, length(sigBoth))
        r = sigBoth(k);
        fprintf('  %-45s  VLCD=%+.4f (p=%.4f)  LCD=%+.4f (p=%.4f)\n', ...
            allRxnsAll{r}, meanDiffVLCD(r), pVLCD(r), meanDiffLCD(r), pLCD(r));
    end
end

%%  9. Summary 
fprintf('Post models loaded: %d/20\n', nPostLoaded);
fprintf('Both BL+Post feasible: %d/20\n', sum(feasBoth));
fprintf('\nJaccard (BL vs Post):\n');
fprintf('  VLCD: %.4f +/- %.4f\n', mean(jaccardPerParticipant(idxVLCD)), std(jaccardPerParticipant(idxVLCD)));
fprintf('  LCD:  %.4f +/- %.4f\n', mean(jaccardPerParticipant(idxLCD)), std(jaccardPerParticipant(idxLCD)));
fprintf('\nCosine (BL vs Post flux):\n');
fprintf('  VLCD: %.4f +/- %.4f\n', mean(cosPerParticipant(idxVLCD & feasBoth)), std(cosPerParticipant(idxVLCD & feasBoth)));
fprintf('  LCD:  %.4f +/- %.4f\n', mean(cosPerParticipant(idxLCD & feasBoth)), std(cosPerParticipant(idxLCD & feasBoth)));
fprintf('\nSignificant flux changes (p<0.05):\n');
fprintf('  VLCD: %d reactions\n', length(sigVLCD));
fprintf('  LCD:  %d reactions\n', length(sigLCD));
fprintf('  Both: %d | VLCD-only: %d | LCD-only: %d\n', ...
    length(sigBoth), length(sigOnlyVLCD), length(sigOnlyLCD));
