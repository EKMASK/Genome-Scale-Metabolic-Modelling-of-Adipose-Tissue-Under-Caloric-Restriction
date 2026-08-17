%% setup_models.m
% Shared setup: loads iMAT baseline models + participant metadata.
% Place this script in the same folder as iMAT_Foguet_baseline_models/
% and iMAT_recon3_baseline_models/.

%% ---- Configuration ----
reconstructionType = 'fouget';  % 'recon3d' or 'fouget'

switch reconstructionType
    case 'recon3d'
        modelDir = 'iMAT_recon3_baseline_models';
    case 'fouget'
        modelDir = 'iMAT_Foguet_baseline_models';
    otherwise
        error('Unknown reconstruction type: %s', reconstructionType);
end

%% ---- Participant metadata (from Yoyo_additional_measurements.txt) ----
participantIDs = {'P01','P02','P03','P04','P05','P06','P08','P09','P10','P11', ...
                  'P12','P13','P14','P19','P20','P21','P22','P23','P25','P26'};

dietGroup = {'VLCD','VLCD','VLCD','VLCD','VLCD','LCD','LCD','LCD','LCD','LCD', ...
             'VLCD','LCD','LCD','LCD','LCD','VLCD','VLCD','LCD','VLCD','VLCD'};

sex = {'F','F','M','F','M','F','F','F','F','M', ...
       'F','F','M','F','M','F','M','M','M','M'};

age = [46, 45, 52, 41, 56, 66, 48, 53, 40, 43, ...
       51, 53, 63, 67, 60, 42, 48, 66, 50, 65];

BMI_baseline = [29.9, 30.1, 31.3, 33.9, 30.2, 32.7, 28.7, 33.9, 33.8, 29.1, ...
                31.2, 36.0, 34.8, 35.6, 33.1, 34.8, 30.96, 31.5, 30.6, 32.1];

weight_baseline = [83.4, 80.9, 99.6, 80.3, 96.7, 82.1, 78.3, 97.6, 100, 94.3, ...
                   85, 102.9, 119.6, 86, 95, 87.4, 107.1, 99.9, 94.2, 114.1];

nParticipants = length(participantIDs);
idxVLCD = strcmp(dietGroup, 'VLCD');
idxLCD  = strcmp(dietGroup, 'LCD');

%% ---- Load baseline context-specific models ----
fprintf('Loading %d baseline models from: %s\n', nParticipants, modelDir);

models = cell(1, nParticipants);
nRxns  = zeros(1, nParticipants);
nMets  = zeros(1, nParticipants);

for i = 1:nParticipants
    % Filename: model_P01_baseline_VLCD.mat (includes diet group)
    fname = fullfile(modelDir, sprintf('model_%s_baseline_%s.mat', ...
        participantIDs{i}, dietGroup{i}));
    
    if ~isfile(fname)
        warning('Not found: %s', fname);
        continue;
    end
    
    tmp = load(fname);
    
    % Models saved as 'tissueModel'
    if isfield(tmp, 'tissueModel')
        models{i} = tmp.tissueModel;
    elseif isfield(tmp, 'model')
        models{i} = tmp.model;
    else
        fn = fieldnames(tmp);
        models{i} = tmp.(fn{1});
    end
    
    nRxns(i) = length(models{i}.rxns);
    nMets(i) = length(models{i}.mets);
    fprintf('  %s (%s): %d rxns, %d mets\n', participantIDs{i}, dietGroup{i}, nRxns(i), nMets(i));
end

nLoaded = sum(~cellfun(@isempty, models));
fprintf('\nLoaded %d/%d models.\n', nLoaded, nParticipants);
fprintf('Reactions:   min=%d, max=%d, mean=%.1f\n', min(nRxns), max(nRxns), mean(nRxns));
fprintf('Metabolites: min=%d, max=%d, mean=%.1f\n', min(nMets), max(nMets), mean(nMets));

if nLoaded == 0
    error('No models loaded! Check that %s/ contains the .mat files.', modelDir);
end