function [finish,pre] = simMouseChowDiet(modelCom, diet, param, paramFile)
o2Id = 'o2[u]';
scfaId = {'ac[u]', 'but[u]', 'ppa[u]'};

undone = false(numel(param.nSim), 1);
for k = 1:numel(param.nSim)
    if ~exist(sprintf('%s_sim%d_sect7save01.mat',param.saveName,param.nC(j),param.nR(k)),'file')
        undone(k) = true;
    end
end

finish = true;
if nnz(undone) == 0
    finish = true;
    pre = struct();
    fprintf('All %d simulations finished already\n', numel(param.nSim));
    return
end

if isfield(param,'scaleFactor')
    sf = param.scaleFactor;
else
    sf = 1;
end

%% Random bound
load(paramFile,'LB','nRand');

%%

options = struct();
options.GRguess = 0.2;
options.algorithm = 2;
options.GRtol = 1e-5;
options.verbFlag = 3;

saveName = 'spatialGut_mouseChow/run';
saveDir = strsplit(saveName,filesep);
if numel(saveDir) > 1
    saveDir = strjoin(saveDir(1:end-1),filesep);
    if ~exist(saveDir,'dir')
        mkdir(saveDir);
    end
end

% initial nutrients (diet)
[yn, id] = ismember(modelCom.infoCom.Mcom, diet(:, 1));
options.C = zeros(numel(modelCom.infoCom.Mcom), 1);
options.C(yn) = cell2mat(diet(id(yn), 2)); %amount for one meal (mmol)
options.C(o2IdCom) = 0.00; %assume low oxygen level in diet

% optionsJ.C(optionsJ.C>0) = 1000;
options.T = param.T;%[0.5 1];
options.dtMuc = param.dt; %every half hour?
options.dtLum = 1/12; %every 5 min
options.dtOxy = 1/60/60/100; %every second
options.O2 = 'o2[u]';
options.O2s = 1;%[1; 0.1];
%wikipedia: mass diffusivity (cm^2/s)
options.O2d = param.O2d * 10^2 * 60^2;%2.1e-5 * 5.^(0:-1:-3) * 10^2 * 60^2;%[2.1e-5; 2.1e-7] * 10^2 * 60^2;
options.R = param.R;%[5;15];
options.Rmuc = param.Rmuc;%[2;2];
options.RmucCom = param.RmucCom;%[2;2];
%Assume same density on the mucosal layer, so total biomass proportional to area.
%assume 0.005 gdw / unit area
options.X = param.X * sf;%0.005 * 2 * pi * optionsJ.RmucCom .* (optionsJ.R-optionsJ.RmucCom/2);
options.L = param.L;
options.saveDsteady = [param.saveName sprintf('_test%02d',jTest)];%sprintf('%s/ab''testSteadyOxy6_08_08_2016/cUT' num2str(nCarbon(jCarbon)) '_x05_Sug01_test' sprintf('%03d',jTest)];
options.algorithm = 2;
options.GRtol = 1e-4;
options.verbFlag = 2;
options.saveFre = 21;
%40 mmHg ~ 0.05 atm from Albenberg Gastroenterology 2014,
%3% ~ 0.03 from Zeitouni, Mole Cell Pediatrics, 2016
options.pO2 = param.pO2; %[1e-4 * 5.^(0:-1:-3);%[1e-2;1e-3];
%E coli OUR from iAF1260, C glut OUR from papers, divided by the
%partial pressure in typical fermentor, 95% air
options.O2ut = param.O2ut;%[0 0 0 18.5 4] / (0.20946 * 0.95) / 10;
options.O2total = param.O2total;
%at 36.9 oC: 770 * exp(1700 * (1/310.05 - 1/289.15)) [atm /(mol/L)]
options.Kpc = 619/1000;
options.Len = param.Len;
options.O2fluxMuc = param.O2fluxMuc * sf;
options.O2fluxLum = param.O2fluxLum * sf;

proton = find(strcmp(metCom,'h[u]'));
pi = find(strcmp(metCom,'pi[u]'));

pre = struct();
pre.options = options;
pre.idUT = idUT;
pre.UT = UT;
pre.metCom = metCom;
pre.metComMW = MW(modelCom.metFormulas(modelCom.metSps==0));
if isfield(param,'infoOnly') && param.infoOnly
    finish = false;
    return
end
param.nR0 = param.nR;
[modelCom.Mcom, ~] = find(modelCom.S(:, modelCom.EXcom(:, 1)));
modelCom.Msp = zeros(size(modelCom.EXcom, 1), numel(modelCom.sps));
for j = 1:size(modelCom.EXcom, 1)
    for k = 1:numel(modelCom.sps)
        if modelCom.EXsp(j, k) > 0
            mJK = find(modelCom.S(:, modelCom.EXsp(j, k)));
            if ~any(mJK == modelCom.Mcom(j)) || numel(mJK) ~= 2
                error('Incorrect metabolite index')
            end
            mJK = mJK(mJK ~= modelCom.Mcom(j));
            modelCom.Msp(j, k) = mJK;
        end
    end
end
comField = {'spATPM', 'spBm', 'rxnSps', 'metSps', 'EXsp', 'EXcom', 'Mcom', 'Msp'};
for j = 1:numel(comField)
    modelCom.indCom.(comField{j}) = modelCom.(comField{j});
end
modelCom = rules2grRules(modelCom,1);
modelCom.indCom.geneSps = 5 * ones(numel(modelCom.genes), 1);
modelCom.infoCom = infoCom2indCom(modelCom, modelCom.indCom, true, modelCom.sps, modelCom.sps);
if ~isfield(param,'parallel') || param.parallel == 0
    for j = 1:numel(param.nC)
        param.nR = param.nR0(undone(:,j));
        for k = 1:numel(param.nR)
            fprintf('nC.nR: %d.%d\n',param.nC(j),param.nR(k));
            modelCom.lb(modelCom.EXsp(modelCom.EXsp~=0)) = round(LB(:,param.nR(k),param.nC(j)),4,'significant');
            modelCom.lb(modelCom.EXsp(o2IdCom,[1 2 4 5])) = -abs(param.O2spRate([1 2 4 5]));
            if isfield(param,'Hex') && ~param.Hex
                modelCom.lb(modelCom.EXsp(proton,modelCom.EXsp(proton,:)~=0)) = 0;
            end
            if isfield(param,'pi') && ~isempty(param.pi)
                modelCom.lb(modelCom.EXsp(pi,modelCom.EXsp(pi,:)~=0)) = -abs(param.pi);
            end
            options.nCR = [param.nC(j) param.nR(k)];
            options.saveDsteady = [param.saveName sprintf('_c%dr%d',param.nC(j),param.nR(k))];
            finish(j,k) = dSteadyOxy7(modelCom, options, solverParam);
        end
    end
else
    for j = 1:numel(param.nC)
        param.nR = param.nR0(undone(:,j));
        parfor k = 1:numel(param.nR)
            fprintf('nC.nR: %d.%d\n',param.nC(j),param.nR(k));
            modelJ = modelCom;
            modelJ.lb(modelJ.EXsp(modelCom.EXsp~=0)) = round(LB(:,param.nR(k),param.nC(j)),4,'significant');
            modelJ.lb(modelJ.EXsp(o2IdCom,[1 2 4 5])) = -abs(param.O2spRate([1 2 4 5]));
            if isfield(param,'Hex') && ~param.Hex
                modelJ.lb(modelJ.EXsp(proton,modelCom.EXsp(proton,:)~=0)) = 0;
            end
            if isfield(param,'pi') && ~isempty(param.pi)
                modelJ.lb(modelJ.EXsp(pi,modelCom.EXsp(pi,:)~=0)) = -abs(param.pi);
            end
            optionsJ = options;
            optionsJ.saveDsteady = [param.saveName sprintf('_c%dr%d',param.nC(j),param.nR(k))];
            optionsJ.nCR = [param.nC(j) param.nR(k)];
            finish(j,k) = dSteadyOxy7(modelJ, optionsJ, solverParam);
        end
    end
end
