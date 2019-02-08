function [finish, options] = simMouseChowDiet(modelCom, diet, options, paramFile)
o2Id = 'o2[u]';
undone = false(numel(options.nSim), 1);
for k = 1:numel(options.nSim)
    if ~exist(sprintf('%s_sim%d_sect7save01.mat',options.saveName,options.nC(j),options.nR(k)),'file')
        undone(k) = true;
    end
end


%% Random bound
load(paramFile,'LB','nRand');

%%



% a scaling factor to avoid tolerance issues for small values
sf = 1000;
% initial nutrients (diet)
[yn, id] = ismember(modelCom.infoCom.Mcom, diet(:, 1));
options.C = zeros(numel(modelCom.infoCom.Mcom), 1);
options.C(yn) = cell2mat(diet(id(yn), 2)) * sf;
% oxygen balance is not accounted. Oxygen availability is controlled by other parameters
options.C(strcmp(modelCom.infoCom.Mcom, o2Id)) = 0.00; 
options.X = options.X * sf;
options.O2fluxMuc = options.O2fluxMuc * sf;
options.O2fluxLum = options.O2fluxLum * sf;
options.T = [2 2 2 3 3 3 3];  % transit time in each intestinal section
options.dtMuc = 0.5;
options.dtLum = 1/12;
options.O2 = o2Id;
options.algorithm = 2;
options.GRtol = 1e-4;
options.verbFlag = 2;
options.saveFre = 21;
options.GRguess = 0.2;
options.verbFlag = 3;

finish = true;
if nnz(undone) == 0
    finish = true;
    fprintf('All %d simulations finished already\n', numel(options.nSim));
    return
elseif isfield(options,'infoOnly') && options.infoOnly
    finish = false;
    return
end

modelCom = rules2grRules(modelCom,1);
modelCom.indCom.geneSps = 5 * ones(numel(modelCom.genes), 1);
modelCom.infoCom = infoCom2indCom(modelCom, modelCom.indCom, true, modelCom.sps, modelCom.sps);
if ~isfield(options,'parallel') || options.parallel == 0
    for j = 1:numel(options.nC)
        options.nR = options.nR0(undone(:,j));
        for k = 1:numel(options.nR)
            fprintf('nC.nR: %d.%d\n',options.nC(j),options.nR(k));
            modelCom.lb(modelCom.EXsp(modelCom.EXsp~=0)) = round(LB(:,options.nR(k),options.nC(j)),4,'significant');
            modelCom.lb(modelCom.EXsp(o2IdCom,[1 2 4 5])) = -abs(options.O2spRate([1 2 4 5]));
            options.nCR = [options.nC(j) options.nR(k)];
            options.saveDsteady = [options.saveName sprintf('_c%dr%d',options.nC(j),options.nR(k))];
            finish(j,k) = dSteadyOxy7(modelCom, options, solverParam);
        end
    end
else
    for j = 1:numel(options.nC)
        options.nR = options.nR0(undone(:,j));
        parfor k = 1:numel(options.nR)
            fprintf('nC.nR: %d.%d\n',options.nC(j),options.nR(k));
            modelJ = modelCom;
            modelJ.lb(modelJ.EXsp(modelCom.EXsp~=0)) = round(LB(:,options.nR(k),options.nC(j)),4,'significant');
            modelJ.lb(modelJ.EXsp(o2IdCom,[1 2 4 5])) = -abs(options.O2spRate([1 2 4 5]));
          
            optionsJ = options;
            optionsJ.saveDsteady = [options.saveName sprintf('_c%dr%d',options.nC(j),options.nR(k))];
            optionsJ.nCR = [options.nC(j) options.nR(k)];
            finish(j,k) = dSteadyOxy7(modelJ, optionsJ, solverParam);
        end
    end
end
