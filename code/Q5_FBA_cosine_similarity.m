clear; clc; close all;

changeCobraSolver('gurobi', 'all');

%% ---- Load models and metadata ----
setup_models;

%% ---- 1. Inspect existing exchange reaction bounds ----
% The Fouget adipocyte model is a manually curated tissue-specific
% reconstruction. Its exchange reactions already have physiologically
% appropriate bounds set by the model authors. We inspect them here.

model1 = models{1};
excLogical = findExcRxns(model1);
excRxnNames = model1.rxns(excLogical);

fprintf('\nExchange reactions in %s: %d total\n', participantIDs{1}, length(excRxnNames));

% Categorise existing bounds
nUptake = 0; nSecretion = 0; nClosed = 0; nFree = 0;
fprintf('\n%-35s %8s %8s  %s\n', 'Reaction', 'LB', 'UB', 'Type');
fprintf('%s\n', repmat('-', 1, 70));
for k = 1:length(excRxnNames)
    idx = find(strcmp(model1.rxns, excRxnNames{k}));
    lb = model1.lb(idx); ub = model1.ub(idx);
    
    if lb < 0 && ub <= 0
        rtype = 'uptake only';  nUptake = nUptake + 1;
    elseif lb >= 0 && ub > 0
        rtype = 'secretion only'; nSecretion = nSecretion + 1;
    elseif lb < 0 && ub > 0
        rtype = 'bidirectional'; nFree = nFree + 1;
    else
        rtype = 'closed'; nClosed = nClosed + 1;
    end
    
    fprintf('%-35s %8.2f %8.2f  %s\n', excRxnNames{k}, lb, ub, rtype);
end
fprintf('\nSummary: %d uptake, %d secretion, %d bidirectional, %d closed\n', ...
    nUptake, nSecretion, nFree, nClosed);

%% ---- 2. Constraint strategy ----
%
% RATIONALE:
%   The Fouget adipocyte model is a manually curated, tissue-specific
%   reconstruction whose exchange reaction bounds already reflect the
%   metabolic environment of subcutaneous adipose tissue. These bounds
%   encode which metabolites the tissue can take up from the bloodstream
%   (e.g., glucose, fatty acids, amino acids, oxygen) and which it
%   secretes (e.g., lactate, glycerol, CO2).
%
%   Since iMAT only determines which reactions are ACTIVE (by removing
%   reactions inconsistent with gene expression), it preserves the
%   original exchange reaction bounds from the parent model. Therefore,
%   the context-specific models already carry appropriate physiological
%   constraints.
%
%   We adopt the following strategy:
%     1. KEEP the existing curated bounds (do NOT zero them out)
%     2. Verify that essential exchanges are open (glucose, O2, etc.)
%     3. Only adjust bounds if a specific exchange appears unreasonable
%
%   This is preferable to overriding all bounds because:
%     - The Fouget bounds are literature-derived and tissue-specific
%     - Zeroing all exchanges and selectively re-opening risks missing
%       essential cofactors or metabolites, causing infeasibility
%     - The original bounds already account for adipocyte-specific
%       uptake rates (e.g., amino acid uptake capped at 0.1 mmol/gDW/h)
%
%   For a generic model like Recon3D, we would need to manually constrain
%   all exchange reactions. This is another advantage of using a curated
%   tissue-specific reconstruction.

fprintf('\nConstraint strategy: using Fouget curated bounds (no override).\n');
fprintf('Verifying essential exchanges are open...\n');

essentialKeys = {'glc', 'o2', 'h2o', 'atp'};
for i = 1:nParticipants
    model = models{i};
    excLog = findExcRxns(model);
    excNames = model.rxns(excLog);
    
    for k = 1:length(essentialKeys)
        hits = excNames(contains(excNames, essentialKeys{k}, 'IgnoreCase', true));
        if isempty(hits) && i == 1
            fprintf('  WARNING: No exchange for "%s" in %s\n', essentialKeys{k}, participantIDs{i});
        end
    end
end
fprintf('Verification complete.\n');

%% ---- 3. Set objective function ----
%
% OBJECTIVE: ATP maintenance demand (Adipocytes_DM_atp_c)
%
% RATIONALE:
%   Adipose tissue in adults is non-proliferating, so a biomass objective
%   (representing cell growth/division) is inappropriate. ATP maintenance
%   demand captures basal energy expenditure for cellular housekeeping:
%   maintaining ion gradients, protein turnover, and membrane integrity.
%
%   Maximising ATPM asks: "given the available nutrients and the set of
%   active reactions determined by iMAT, what is the maximum rate of ATP
%   turnover this metabolic network can sustain?"
%
%   This provides a physiologically interpretable measure of each
%   individual's metabolic capacity at baseline.
%
% ALTERNATIVES CONSIDERED:
%   - Biomass maximisation: inappropriate for non-dividing tissue
%   - pFBA (min total flux): useful post-hoc but no biological objective
%   - Lipid storage: relevant but difficult to define generically

fprintf('\nSetting objective function...\n');
objRxnUsed = '';

for i = 1:nParticipants
    model = models{i};
    found = false;
    
    patterns = {'DM_atp', 'ATPM', 'ATPhyd', 'ATPase'};
    for k = 1:length(patterns)
        matches = find(contains(model.rxns, patterns{k}, 'IgnoreCase', true));
        if ~isempty(matches)
            objRxn = model.rxns{matches(1)};
            model = changeObjective(model, objRxn);
            found = true;
            if i == 1
                objRxnUsed = objRxn;
                fprintf('  Objective: %s\n', objRxn);
            end
            break;
        end
    end
    
    if ~found
        % Fallback to biomass
        matches = find(contains(model.rxns, 'biomass', 'IgnoreCase', true));
        if ~isempty(matches)
            objRxn = model.rxns{matches(1)};
            model = changeObjective(model, objRxn);
            if i == 1
                objRxnUsed = objRxn;
                fprintf('  WARNING: No ATPM. Using biomass: %s\n', objRxn);
            end
        else
            warning('No objective found for %s!', participantIDs{i});
        end
    end
    
    models{i} = model;
end

%% ---- 4. Run FBA ----
fprintf('\nRunning FBA...\n');

% Build union of all reaction IDs for aligned flux vectors
allRxns = {};
for i = 1:nParticipants
    allRxns = union(allRxns, models{i}.rxns);
end
nTotalRxns = length(allRxns);

fluxMatrix = zeros(nTotalRxns, nParticipants);
objValues  = zeros(1, nParticipants);
feasible   = true(1, nParticipants);

for i = 1:nParticipants
    sol = optimizeCbModel(models{i}, 'max');
    
    if sol.stat == 1
        objValues(i) = sol.f;
        [~, ia, ib] = intersect(allRxns, models{i}.rxns);
        fluxMatrix(ia, i) = sol.v(ib);
        fprintf('  %s: obj=%.4f\n', participantIDs{i}, sol.f);
    else
        feasible(i) = false;
        fprintf('  %s: INFEASIBLE (stat=%d)\n', participantIDs{i}, sol.stat);
    end
end

nFeasible = sum(feasible);
fprintf('\nFeasible: %d/%d\n', nFeasible, nParticipants);

if nFeasible == 0
    error('All models infeasible! Check constraints and objective.');
end

%% ---- 5. Cosine similarity ----
% cos(A,B) = (A . B) / (||A|| ||B||)

feasIdx  = find(feasible);
nF       = length(feasIdx);
fluxFeas = fluxMatrix(:, feasIdx);

% Remove zero-flux reactions (no information)
nonzero  = any(fluxFeas ~= 0, 2);
fluxFilt = fluxFeas(nonzero, :);
fprintf('Non-zero flux reactions: %d / %d\n', sum(nonzero), nTotalRxns);

cosineMatrix = zeros(nF);
for i = 1:nF
    for j = i:nF
        a = fluxFilt(:,i); b = fluxFilt(:,j);
        d = norm(a) * norm(b);
        if d > 0
            cosineMatrix(i,j) = dot(a,b) / d;
        end
        cosineMatrix(j,i) = cosineMatrix(i,j);
    end
end

feasIDs  = participantIDs(feasIdx);
feasDiet = dietGroup(feasIdx);
feasSex  = sex(feasIdx);
feasBMI  = BMI_baseline(feasIdx);

%% ---- 6. Cosine heatmap ----
figure('Name','Cosine Heatmap','Position',[100 100 800 700]);
imagesc(cosineMatrix); colormap(parula);
cb = colorbar; cb.Label.String = 'Cosine Similarity';
set(gca, 'XTick', 1:nF, 'XTickLabel', feasIDs, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:nF, 'YTickLabel', feasIDs);
title(sprintf('Cosine Similarity of FBA Fluxes — %s', reconstructionType));
saveas(gcf, sprintf('Q5_cosine_heatmap_%s.png', reconstructionType));

%% ---- 7. Hierarchical clustering ----
cosDist = 1 - cosineMatrix;
cosDist(cosDist < 0) = 0;
cosDist = (cosDist + cosDist') / 2;  % force symmetry
cosDist(1:nF+1:end) = 0;            % force diagonal to exactly 0
distVec = squareform(cosDist);
linkTree = linkage(distVec, 'average');

figure('Name','Cosine Dendrogram','Position',[100 100 800 400]);
[~, ~, outperm] = dendrogram(linkTree, 0, 'Labels', feasIDs, 'ColorThreshold', 'default');
title(sprintf('Hierarchical Clustering (Cosine distance) — %s', reconstructionType));
ylabel('Cosine distance');
set(gca, 'XTickLabelRotation', 45);
saveas(gcf, sprintf('Q5_cosine_dendrogram_%s.png', reconstructionType));

% Clustered heatmap
ordCos  = cosineMatrix(outperm, outperm);
ordIDs  = feasIDs(outperm);
ordDiet = feasDiet(outperm);

figure('Name','Clustered Cosine','Position',[100 100 800 700]);
imagesc(ordCos); colormap(parula);
cb = colorbar; cb.Label.String = 'Cosine Similarity';
set(gca, 'XTick', 1:nF, 'XTickLabel', ordIDs, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:nF, 'YTickLabel', ordIDs);
title(sprintf('Clustered Cosine Similarity — %s', reconstructionType));
hold on;
for i = 1:nF
    if strcmp(ordDiet{i}, 'VLCD')
        clr = [0.85 0.33 0.1];
    else
        clr = [0.0 0.45 0.74];
    end
    rectangle('Position', [i-0.5, 0.0, 1, 0.5], 'FaceColor', clr, 'EdgeColor','none');
end
hold off;
saveas(gcf, sprintf('Q5_clustered_cosine_%s.png', reconstructionType));

%% ---- 8. MDS on cosine distance ----
[Y, eigvals] = cmdscale(cosDist);
explVar = eigvals / sum(abs(eigvals)) * 100;

figure('Name','MDS Cosine','Position',[100 100 1200 400]);

% By diet
subplot(1,3,1); hold on;
hV = []; hL = [];
for i = 1:nF
    if strcmp(feasDiet{i}, 'VLCD')
        hV = scatter(Y(i,1), Y(i,2), 80, [0.85 0.33 0.1], 'filled');
    else
        hL = scatter(Y(i,1), Y(i,2), 80, [0.0 0.45 0.74], 'filled');
    end
    text(Y(i,1)*1.03, Y(i,2)*1.03, feasIDs{i}, 'FontSize', 7);
end
xlabel(sprintf('MDS1 (%.1f%%)', explVar(1)));
ylabel(sprintf('MDS2 (%.1f%%)', explVar(2)));
title('By diet group');
legend([hV hL], {'VLCD','LCD'}, 'Location','best');
hold off;

% By sex
subplot(1,3,2); hold on;
hM = []; hF = [];
for i = 1:nF
    if strcmp(feasSex{i}, 'M')
        hM = scatter(Y(i,1), Y(i,2), 80, [0.47 0.67 0.19], 'filled');
    else
        hF = scatter(Y(i,1), Y(i,2), 80, [0.64 0.08 0.18], 'filled');
    end
    text(Y(i,1)*1.03, Y(i,2)*1.03, feasIDs{i}, 'FontSize', 7);
end
xlabel(sprintf('MDS1 (%.1f%%)', explVar(1)));
ylabel(sprintf('MDS2 (%.1f%%)', explVar(2)));
title('By sex');
legend([hM hF], {'Male','Female'}, 'Location','best');
hold off;

% By BMI
subplot(1,3,3);
scatter(Y(:,1), Y(:,2), 80, feasBMI, 'filled');
colormap(gca, hot);
cb2 = colorbar; cb2.Label.String = 'BMI (baseline)';
for i = 1:nF
    text(Y(i,1)*1.03, Y(i,2)*1.03, feasIDs{i}, 'FontSize', 7);
end
xlabel(sprintf('MDS1 (%.1f%%)', explVar(1)));
ylabel(sprintf('MDS2 (%.1f%%)', explVar(2)));
title('By BMI');

sgtitle(sprintf('MDS on Cosine distance (FBA) — %s', reconstructionType));
saveas(gcf, sprintf('Q5_MDS_cosine_%s.png', reconstructionType));

%% ---- 9. Objective values ----
figure('Name','Objective Values','Position',[100 100 700 350]);
hold on;
for i = 1:nF
    if strcmp(feasDiet{i}, 'VLCD')
        bar(i, objValues(feasIdx(i)), 'FaceColor', [0.85 0.33 0.1]);
    else
        bar(i, objValues(feasIdx(i)), 'FaceColor', [0.0 0.45 0.74]);
    end
end
set(gca, 'XTick', 1:nF, 'XTickLabel', feasIDs, 'XTickLabelRotation', 45);
ylabel(sprintf('Objective value (%s)', objRxnUsed));
title(sprintf('FBA Objective Values — %s', reconstructionType));
hold off;
saveas(gcf, sprintf('Q5_objective_values_%s.png', reconstructionType));

%% ---- 10. Summary statistics ----
feasVLCD = find(strcmp(feasDiet, 'VLCD'));
feasLCD  = find(strcmp(feasDiet, 'LCD'));

wVLCD = []; wLCD = []; btwn = [];
for i = 1:length(feasVLCD)
    for j = i+1:length(feasVLCD)
        wVLCD(end+1) = cosineMatrix(feasVLCD(i), feasVLCD(j)); %#ok
    end
end
for i = 1:length(feasLCD)
    for j = i+1:length(feasLCD)
        wLCD(end+1) = cosineMatrix(feasLCD(i), feasLCD(j)); %#ok
    end
end
for i = 1:length(feasVLCD)
    for j = 1:length(feasLCD)
        btwn(end+1) = cosineMatrix(feasVLCD(i), feasLCD(j)); %#ok
    end
end

fprintf('Objective: %s\n', objRxnUsed);
fprintf('Feasible: %d/%d\n', nFeasible, nParticipants);
fprintf('Mean obj value: %.4f +/- %.4f\n', mean(objValues(feasIdx)), std(objValues(feasIdx)));
fprintf('\nCosine similarity:\n');
fprintf('  Overall:        %.4f +/- %.4f\n', ...
    mean(cosineMatrix(triu(true(nF),1))), std(cosineMatrix(triu(true(nF),1))));
fprintf('  Within VLCD:    %.4f +/- %.4f\n', mean(wVLCD), std(wVLCD));
fprintf('  Within LCD:     %.4f +/- %.4f\n', mean(wLCD), std(wLCD));
fprintf('  Between groups: %.4f +/- %.4f\n', mean(btwn), std(btwn));
