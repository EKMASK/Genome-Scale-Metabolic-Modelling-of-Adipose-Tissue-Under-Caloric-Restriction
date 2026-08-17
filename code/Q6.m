clear; clc; close all;
%% Group assignment 1.3
%% COBRA initialization
initCobraToolbox(false);
changeCobraSolver('gurobi', 'all');
%% Question 6: create context-specific models for post-intervention data with Fouget using iMAT
% load Foguet Adipocyte model
fa = readCbModel('../data/Foguet_adipocyte.mat');
% load expression
exp = readtable('../data/Adipose_Expression_Fouget.txt', 'FileType', 'text');
geneIDs_raw = exp{:,1};
sampleNames = exp.Properties.VariableNames(2:end);
X = table2array(exp(:,2:end)); % genes x samples

% keep post-intervention
post = contains(sampleNames, 'post');
postNames = sampleNames(post);
X_post = X(:, post);
fprintf("Post-intervention samples: %d\n", numel(postNames));
%% building context-specific models
% output models into a folder
outdir = '../data/iMAT_Foguet_post_models';
if ~exist(outdir, 'dir')
    mkdir(outdir); 
end

for k = 1:numel(postNames)
    name = postNames{k};
    fprintf('\nBuilding model for %s (%d/%d) ... \n', name, k, numel(postNames));
    % gene expression vector per participant
    gexp = X_post(:,k);
    % handle missing values
    gexp(isinf(gexp)) = NaN;
    if any(isnan(gexp))
        gexp = fillmissing(gexp,'constant',0);
    end
    
    % map gene expression to reaction expression
    exp_data = struct();
    exp_data.gene = cellstr(string(geneIDs_raw));
    exp_data.value = double(gexp(:));
    rxnExp = mapExpressionToReactions(fa, exp_data);
    fprintf('rxnExp: n=%d, -1=%d, NaN=%d, >=0=%d, min=%g, max=%g\n', ...
    numel(rxnExp), sum(rxnExp==-1), sum(isnan(rxnExp)), sum(rxnExp>=0), ...
    min(rxnExp(rxnExp>=0)), max(rxnExp(rxnExp>=0)));
    % thresholds (percentile)
    vaild = rxnExp(~isnan(rxnExp) & rxnExp >=0);
    low_thres = prctile(vaild, 25);
    high_thres = prctile(vaild, 75);
    % iMAT settings
    options = struct();
    options.solver = 'iMAT';
    options.expressionRxns = rxnExp;
    options.threshold_lb = low_thres;
    options.threshold_ub = high_thres;
    options.tol = 1e-4;
    options.logfile = 'MILPlog';
    % build tissue model with iMAT
    tissueModel = createTissueSpecificModel(fa, options);
    fprintf('tissueModel: rxns=%d, mets=%d, genes=%d\n', ...
        numel(tissueModel.rxns), numel(tissueModel.mets), numel(tissueModel.genes));
    % save models
    safeName = regexprep(name, '[^a-zA-Z0-9_]','_');
    save(fullfile(outdir, ['model_' safeName '.mat']), 'tissueModel');
end

fprintf("all context-specific models are generated in %s\n", outdir);