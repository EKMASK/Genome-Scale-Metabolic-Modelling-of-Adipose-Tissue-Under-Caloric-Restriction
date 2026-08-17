clear; clc; close all;
%% Group assignment 1.3
%% COBRA initialization
initCobraToolbox(false);
changeCobraSolver('gurobi', 'all');
%% Question 2: Build 20 context-specific models with recon3D using iMAT
% load generic model
recon3 = readCbModel('../data/Recon3DModel.mat');
% load expression
exp = readtable('../data/Adipose_Expression_Recon3.txt', 'FileType', 'text');
geneIDs_raw = exp{:,1};
sampleNames = exp.Properties.VariableNames(2:end);
X = table2array(exp(:,2:end)); % genes x samples

% keep baseline (20 participants)
baseline = contains(sampleNames, 'baseline');
baselineNames = sampleNames(baseline);
X_base = X(:, baseline);
fprintf("Baseline samples: %d\n", numel(baselineNames));
%% building context-specific models
% output models into a folder
outdir = '../data/iMAT_recon3_baseline_models';
if ~exist(outdir, 'dir')
    mkdir(outdir); 
end

for k = 1:numel(baselineNames)
    name = baselineNames{k};
    fprintf('\nBuilding model for %s (%d/%d) ...\n', name, k, numel(baselineNames));
    % gene expression vector per participant
    gexp = X_base(:,k);
    % handle missing values
    gexp(isinf(gexp)) = NaN;
    if any(isnan(gexp))
        gexp = fillmissing(gexp,'constant',0);
    end
    
    % map gene expression to reaction expression
    exp_data = struct();
    exp_data.gene = cellstr(string(geneIDs_raw));
    exp_data.value = double(gexp(:));
    rxnExp = mapExpressionToReactions(recon3, exp_data);
    % thresholds (percentile)
    low_thres = prctile(rxnExp(~isnan(rxnExp)), 25);
    high_thres = prctile(rxnExp(~isnan(rxnExp)), 75);
    % iMAT settings
    options = struct();
    options.solver = 'iMAT';
    options.expressionRxns = rxnExp;
    options.threshold_lb = low_thres;
    options.threshold_ub = high_thres;
    options.tol = 1e-4;
    options.logfile = 'MILPlog';
    % build tissue model with iMAT
    tissueModel = createTissueSpecificModel(recon3, options);
    fprintf('tissueModel: rxns=%d, mets=%d, genes=%d\n', ...
        numel(tissueModel.rxns), numel(tissueModel.mets), numel(tissueModel.genes));
    % save models
    safeName = regexprep(name, '[^a-zA-Z0-9_]','_');
    save(fullfile(outdir, ['model_' safeName '.mat']), 'tissueModel');
end
fprintf("all context-specific models are generated in %s\n", outdir);
