clear; clc; close all;
%% Group assignment 1.3
%% COBRA initialization
initCobraToolbox(false);
changeCobraSolver('gurobi', 'all');
%% Question 1
%% Adipose Whole Gene Expression Dataset
% load data (whole gene expression)
WholeGenome = readtable('../data/Adipose_Expression_WholeGenome.txt', 'FileType', 'text');
sampleNames = WholeGenome.Properties.VariableNames(3:end);
X = table2array(WholeGenome(:,3:end));
fprintf('Matrix size: %d genes x %d samples\n;', size(X,1), size(X,2));

% perform PCA
X_t = X'; %transpose to samples x genes
% Diagnose missing or invalid values
fprintf('NaNs in X_t: %d\n', sum(isnan(X_t(:))));
fprintf('Infs in X_t: %d\n', sum(isinf(X_t(:))));

% Impute: simplest (0) or use column mean
X_t(isinf(X_t)) = NaN; %treat inf as missing
X_t = fillmissing(X_t,'constant',0);  

% PCA automatically mean-centers variables
[coeff, score, latent, ~, explained] = pca(X_t);
fprintf('Explained variance PC1 = %.2f%%\n', explained(1));
fprintf('Explained variance PC2 = %.2f%%\n', explained(2));

% extract group information from column names
n = length(sampleNames);
isBaseline = contains(sampleNames, 'baseline');
isPost = contains(sampleNames, 'post');
isVLCD = contains(sampleNames, 'VLCD');
isLCD = contains(sampleNames, 'LCD') & ~isVLCD;

% baseline vs post-intervention
figure;
gscatter(score(:,1), score(:,2), isBaseline, 'br', 'o', 8);
xlabel(sprintf('PC1 (%.2f%%)', explained(1)));
ylabel(sprintf('PC2 (%.2f%%)', explained(2)));
title('Whole Genome PCA: Baseline vs Post');
legend({'Post','Baseline'});
grid on;

% LCD vs VLCD
figure;
groupDiet = categorical(isVLCD,[0 1],{'LCD','VLCD'});
gscatter(score(:,1), score(:,2), groupDiet, 'gm', 'o', 8);
xlabel(sprintf('PC1 (%.2f%%)', explained(1)));
ylabel(sprintf('PC2 (%.2f%%)', explained(2)));
title('Whole Genome PCA: Diet Groups');
grid on;
%% Yoyo Dataset
% load data (yoyo extra features)
yoyo = readtable('../data/Yoyo_additional_measurements.txt');
% use sex feature
sex = yoyo.sex;
sex_full = repelem(sex,2);
figure;
gscatter(score(:,1), score(:,2), sex_full);
title('PCA coloured by Sex');

% use weight loss feature (baseline - post intervention)
weight = yoyo.weight_baseline - yoyo.weight_post_intervetnion;
weight_full = repelem(weight,2);
[r,p] = corr(score(:,1), weight_full);
fprintf('Correlation PC1 vs weight loss: r = %.3f, p = %.3f\n', r, p);

% use BMI feature
BMI_full = repelem(yoyo.BMI_baseline,2);
figure;
scatter(score(:,1), score(:,2), 80, BMI_full, 'filled');
colorbar;
title('PCA coloured by Baseline BMI');
