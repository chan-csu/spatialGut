function finish = simSpatialGutExample(dataFile, simIDs)
% finish = simSpatialGutExample(dataFile, nSimulations)
%
% Example setup of the parameters and call of spatialGut
%
% INPUTS:
%    dataFile         .mat file containing the following variables:
%              model:      a COBRA community model with *.infoCom and *.indCom. See getMultiSpeciesModelId.m 
%              rxnEXsp:    a cell array of reaction IDs for the organism-community exchange reaction
%                          whose uptake bounds are changed according to rxnEXspLBs in each simulation
%              rxnEXspLBs: #rxnEXsp-by-#simulations matrix. Each column is one set of uptake bounds for rxnEXsp 
%              diet:       cell array of nutrients available to the microbiota. 
%                          1st column: community metabolite IDs
%                          2nd column: available amount
%    simIDs    index vector, simulations to perform using bounds in rxnEXspLBs(:, simIDs)
%
% OUTPUT:
%    finish           true if all simulations have finished

if nargin == 0 || isempty(dataFile)
    dataFile = 'simulationData.mat';
end
if nargin < 2
    simIDs = 1:200;
end

%% data needed from the dataFile
data = load(dataFile);
% 5-member community model with *.infoCom and *.indCom
model5 = data.model;
% organism-community exchange reaction IDs 
rxnEXsp = data.rxnEXsp;
% sets of uptake bounds for rxnEXsp, one column one set of uptake bounds
rxnEXspLBs = data.rxnEXspLBs;
% cell array of nutrients available to the microbiota. 
% 1st column: community metabolite IDs; 2nd column: available amount
chowDiet = data.diet;

if size(rxnEXspLBs, 1) ~= numel(rxnEXsp)
    error('rxnEXspLBs must have the number of rows equal to the number of elements in rxnEXsp')
end

% mucosal microbial biomass density
XmucExp = [6.75e-6,2.7e-6,1.5e-5,1.95e-4,5.25e-4,7.5e-6,9e-6];
% total oxygen availability in each section
o2Total = 0.159411 * 1e-5 * ones(1, 7);
% fraction of oxygen available to the luminal microbiota in each intestinal section
o2LumFrac = [0.20, 0.15, 0.15, 0.010, 0.050, 0.050, 0];
% specific oxygen uptake rate for each member
o2spRate = [0 0 0 2 2];

% avoid magnitude close to tolerance level by scaling up substrate available and biomass density of the community
scalingFactor = 1000;

parallel = 0;

options = struct();
options.X = XmucExp * scalingFactor;
% oxygen available to the mucosal microbiota in each section
options.o2fluxMuc = o2Total(:)' .* (1 - o2LumFrac(:)') * scalingFactor;
% oxygen available to the luminal microbiota in each section
options.o2fluxLum = o2Total(:)' .* o2LumFrac(:)' * scalingFactor;

% initial nutrients (diet)
[yn, id] = ismember(model5.infoCom.Mcom, chowDiet(:, 1));
options.C = zeros(numel(model5.infoCom.Mcom), 1);
options.C(yn) = cell2mat(chowDiet(id(yn), 2)) * scalingFactor;

% oxygen met ID
options.o2Id = 'o2[u]';
% index of oxygen among community metabolites
o2indCom = strcmp(model5.infoCom.Mcom, options.o2Id);
% oxygen balance is not accounted. Oxygen availability is controlled by other parameters
options.C(o2indCom) = 0.00; 

% transit time in each intestinal section
options.T = [2 2 2 3 3 3 3];
% name for the files to be saved
options.saveName = ['simSpatialGutExample' filesep 'default_params'];

% SteadyCom parameters
options.algorithm = 2;
options.GRtol = 1e-4;
options.GRguess = 0.2;
options.verbFlag = 3;

% the default 1e-9 tolerance tends to cause numerical instability in large
% community models
changeCobraSolverParams('LP', 'feasTol', 1e-8);
solverParam = struct();

% IDs for organism-community exchange reactions
rxnEXspId = findRxnIDs(model5, rxnEXsp);

% check previous progress
undone = false(numel(simIDs), 1);
for k = 1:numel(simIDs)
    if ~exist(sprintf('%s_sim%d_sect7save01.mat',options.saveName, simIDs(k)),'file')
        undone(k) = true;
    end
end

% check if all simulations have already been finished
finish = true;
if ~any(undone)
    finish = true;
    fprintf('All %d simulations finished already\n', numel(options.nSim));
    return
elseif isfield(options,'infoOnly') && options.infoOnly
    % return options structure only
    finish = false;
    return
end
simIDs = simIDs(undone);

if ~exist('parallel', 'var') || parallel == 0
    % single thread computation
    for j = 1:numel(simIDs)
        optionsJ = options;
        optionsJ.nSim = simIDs(j);
        fprintf('Sim #%d:\n', optionsJ.nSim);
        % assign bounds
        % (many random double numbers using all 16 digits seem to cause trouble to the solver)
        model5.lb(rxnEXspId) = round(rxnEXspLBs(:, simIDs(j)), 4, 'significant');
        % members that have oxygen uptake reactions
        spWtO2 = model5.indCom.EXsp(o2indCom, :) > 0;
        % assign oxygen specific uptake rate assigned in the option
        model5.lb(model5.indCom.EXsp(o2indCom, spWtO2)) = -abs(o2spRate(spWtO2));
        optionsJ.saveName = [optionsJ.saveName sprintf('_sim%d', simIDs(j))];
        spatialGut(model5, optionsJ, solverParam);
    end
else
    % parallel computation
    if isempty(gcp('nocreate'))
        parpool;
    end
    parfor j = 1:numel(simIDs)
        optionsJ = options;
        optionsJ.nSim = simIDs(j);
        fprintf('Sim #%d:\n', optionsJ.nSim);
        modelJ = model5;
        % assign bounds
        % (many random double numbers using all 16 digits seem to cause trouble to the solver)
        modelJ.lb(rxnEXspId) = round(rxnEXspLBs(:, simIDs(j)), 4, 'significant');
        spWtO2 = modelJ.indCom.EXsp(o2indCom, :) > 0;
        % assign oxygen specific uptake rate assigned in the option
        modelJ.lb(modelJ.indCom.EXsp(o2indCom, spWtO2)) = -abs(o2spRate(spWtO2));
        optionsJ.saveName = [optionsJ.saveName sprintf('_sim%d', simIDs(j))];
        spatialGut(modelJ, optionsJ, solverParam);
    end
end
end