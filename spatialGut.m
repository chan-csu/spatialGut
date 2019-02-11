function spatialGut(model, options, solverParam)
% Simulate the luminal and the mocusally adherent microbiota and their metabolism 
% along the intestines using a community metabolic model. 
% In each section, the luminal microbiota changes dynamically with time and grow 
% while the mocusal microbiota is simulated at a community steady-state
% where the newly grown cells are shed into the lumen and become luminal microbes.
% No output for this function. Results are saved in a set of files. Use
% plotSpatialGutResults.m to retrieve and plot the data.
%
% USAGE:
%    spatialGut(model, options, solverParam)
%
% INPUTS:
%    model:      COBRA community model (template for all communities)
%    options:    structure with the following fields:
%                - parameters that must be supplied
%                  * C - mCom-by-1 vector of the initial amounts of community metabolites 
%                        entering into the first intestinal section array
%                  * X - nSect-by-1 vector of mucosal community biomass for each section of the intestines
%                  * T - nSect-by-1 vector of transit time for each section of the intestines
%                  * o2fluxMuc - nSect-by-1 vector of oxygen community uptake bound for 
%                                the mucosal microbiota in each section of the intestines
%                  * o2fluxLum - nSect-by-1 vector of oxygen community uptake bound for 
%                                the luminal microbiota in each section of the intestines
%               - optional parameters
%                  * dtMuc - time step for each consecutive simulation of the mucosal microbiota (0.5);
%                  * dtLum - time step for each consecutive simulation of the luminal microbiota (1/12)
%                  * o2Id - oxygen community metabolite ID (o2[u])
%                  * nSim - the number of this simulation, for printing information only (1)
%                  * saveDetail - true to save all details including fluxes, solutions and 
%                                 concentration vectors at each step (false)
%                  * saveName - filename for saving. Can include directory ('spatialGutSim/test') 
%                  * saveFre - steps per savefile (21)
%
%    solverParam: COBRA solver parameters

%% Check inputs
tReal = tic;
tol = 1e-10;        
[resultTmp, dtMuc, dtLum, C, T, X, saveName, saveFre, o2Id, o2fluxMuc, o2fluxLum, nSim, saveDetail] = ...
    getSpatialGutParams({'resultTmp','dtMuc','dtLum','C','T','X','saveName','saveFre', ...
    'o2Id','o2fluxMuc','o2fluxLum','nSim','saveDetail'}, options, model);
if nargin < 3
    solverParam = struct();
end
[feasTol, ~] = getCobraSolverParams('LP', {'feasTol'; 'optTol'}, solverParam);

% check the directory for saving
saveDir = fileparts(saveName);
if ~isempty(saveDir) && ~exist(saveDir, 'dir')
    mkdir(saveDir)
end
% check if necessary parameters are provided
field = {'C', 'T', 'X', 'o2fluxMuc', 'o2fluxLum'};
for j = 1:numel(field)
    if ~isfield(options,field{j})
        error('%s must be provided in the option structure.',field{j});
    end
end
% number of community metabolites
nCom = size(model.infoCom.Mcom, 1);
if numel(C) ~= nCom
    error('Size of options.C (%d) not equal to the number of community metabolites (%d).', numel(C), nCom);
elseif numel(X) ~= numel(T)
    error('Size of options.X not equal to size of options.T, which is the number of sections in the intestines to be simulated.');
end
nSect = numel(X);
if ~isfield(model, 'infoCom')
    error('model.infoCom does not exist. Get model.infoCom using getMultiSpeciesModelId.m')
end
if ~isfield(model, 'indCom')
    model.indCom = infoCom2indCom(model);
end
%% Initialization
% ensure dtMuc is an integer multiple of dtLum 
dtLum = dtMuc / ceil(dtMuc / dtLum);
% number of simulations for the luminal microbiota per simulation for the mucosal microbiota
nStepLum = round(dtMuc / dtLum);
% index of O2 among the community metabolites
O2 = find(strcmp(model.infoCom.Mcom, o2Id));
% number of organisms in the community
nSp = numel(model.infoCom.spAbbr);
[m, n] = size(model.S);
% solution structure for mucosal and luminal microbiota
resMuc = repmat(resultTmp,1,saveFre);
resLum = repmat(resultTmp,nSp,nStepLum,saveFre);
% O2 flux bounds and values
o2utMuc = zeros(2,saveFre);
o2utLum = zeros(2,nStepLum,saveFre);
% flux from resMuc may need rescaling after optimization. Record separately.
fluxMuc = zeros(size(model.S,2),saveFre);
% time at which the result was got
time = zeros(saveFre,1);
timeLum = zeros(nStepLum,saveFre); 
% metabolite level
[Ct] = deal(zeros(size(C,1), nStepLum, saveFre));

% biomass (relative abundance)
Xmuc = zeros(nSp, saveFre);
[Xlum,optOrder] = deal(zeros(nSp, nStepLum, saveFre));
% growth rate
GRmuc = zeros(saveFre, 1);
GRlum = zeros(nSp, nStepLum, saveFre);
% whether the species is dying (if it is not growing and cannot satisfy non-growth associated ATP maintenance
dying = false(nSp, nStepLum, saveFre);

% Solve SteadyCom or not in each step (if the uptake condition remains the same, no need)
scSolve = false(saveFre,1);
scFinish = true(saveFre, 1);

% whether the entire simulation has finished
finish = false;

%Variabile for the current step
%nothing in the lumen initially
XlumCur = zeros(nSp,1);

%O2 level capture by the two specific variables o2MucCur and o2LumCur
C(O2) = 10000;
%simulation time
t = 0;
%dt0 = dt;
%no. of total steps
kTotal = 1;
%expected total number of steps
Nstep = ceil(sum(T) / dtMuc);
digit = floor(log10(Nstep))+1;
%previous uptake bound for mucosal community
ubMucPrev = inf(nCom,1);
%options
[optionsMuc, optionsLum] = deal(options);

%% Detect previous saved progress if any
j0 = 0;
[kTotal0, kStep0] = deal(1);
kSave0 = 0;
while true
    if exist(sprintf(['%s_sect%dsave%0' num2str(digit) 'd.mat'], ...
            saveName, j0, kSave0+1), 'file')
        kSave0 = kSave0 + 1;
    elseif exist(sprintf(['%s_sect%dsave%0' num2str(digit) 'd.mat'], ...
            saveName, j0+1, 1), 'file')
        j0 = j0 + 1;
        kSave0 = 0;
    else
        break
    end
end
if kSave0 > 0 || j0 >= 1
    %this load all current variables from the latest simulation
    while j0 >= 1
        try
            load(sprintf(['%s_sect%dsave%0' num2str(digit) 'd.mat'], ...
                saveName, j0, kSave0),'t','XlumCur','C', ...
                'j0','kTotal0','kSave0','nextJ','resMuc','kStep0');
            break
        catch ME
            ME.message
            j0 = j0 - 1;
        end
    end
    if j0 >= 1
        kStepPrev = kStep0;
        kTotal = kTotal0 + 1;
        if nextJ
            % if solving for the next intestinal section
            j0 = j0 + 1;
            kSave0 = 1;
        else
            kSave0 = kSave0 + 1;
        end
        dataPrev = load(sprintf('%s_pre.mat', saveName),'SpRate','SpRateUt','SpRateEx');
        Id = findRxnIDs(model,dataPrev.SpRate);
        model.lb(Id) = dataPrev.SpRateUt;
        model.ub(Id) = dataPrev.SpRateEx;
    end
end
if j0 == 0
    j0 = 1;
    kSave0 = 1;
    nextJ = true;
    SpRate = model.rxns(model.indCom.EXsp(model.indCom.EXsp~=0));
    SpRateUt = model.lb(model.indCom.EXsp(model.indCom.EXsp~=0));
    SpRateEx = model.ub(model.indCom.EXsp(model.indCom.EXsp~=0));
    save(sprintf('%s_pre.mat', saveName), 'options','solverParam',...
        'SpRate', 'SpRateUt', 'SpRateEx', 'nStepLum', 'dtMuc', 'dtLum',...
        'digit', 'Nstep', 'feasTol', 'tol', 'o2fluxMuc', 'o2fluxLum');
    clear SpRate SpRateUt SpRateEx
end

%% Simulate the mucosal microbiota
init = true;
fprintf('Sim #%d. Start from section %d, save %d...\n', nSim, j0,kSave0);
% for each section in the intestines
for j = j0:nSect
    
    %counter for step and number of save files
    kStep = 1;
    if init
        % first run of the current call of the function. Start with the detected latest saved progress
        kSave = kSave0;
        init = false;
    else
        % >= 2nd run of the current call. Start with the save file #1
        kSave = 1;
    end
    
    % store the initial time
    time(kStep) = t;
    % assign the mucosal microbial biomass
    optionsMuc.BMweight = X(j);
    %simulation time for the current section
    maxT = sum(T(1:j)) - dtMuc * 1e-5;
    
    while true
        % store the O2 uptake bound
        o2utMuc(1,kStep) = o2fluxMuc(j);
        % community uptake bounds for the mucosal microbiota
        ubMucCur = min(10000 * ones(nCom, 1), C./dtMuc);
        %
        ubMucCur(O2) = o2fluxMuc(j);
        ubMucCur(ubMucCur < feasTol) = 0; %may cause numerical problem for such small ub
        % difference in uptake bounds
        ubDiff = abs(ubMucCur - ubMucPrev);
        if ~(isequal(ubMucCur ~= 0, ubMucPrev ~= 0)) || sum(ubDiff(ubMucPrev ~= 0) ./ ubMucPrev(ubMucPrev ~= 0)) > 0.01
            model.ub(model.indCom.EXcom(:,1)) = ubMucCur;
            % there are differences. Solve SteadyCom
            scSolve(kStep) = true;
            [~, resMuc(kStep)] = SteadyCom(model, optionsMuc, solverParam);
            if ~strcmp(resMuc(kStep).stat,'optimal')
                scFinish(kStep) = false;
                resMuc(kStep).Ex = zeros(nCom,1);
                resMuc(kStep).Ut = zeros(nCom,1);
                resMuc(kStep).GRmax = 0;
                resMuc(kStep).BM = zeros(nSp,1);
                resMuc(kStep).vBM = zeros(nSp,1);
                resMuc(kStep).flux = zeros(n,1);
            end
            if sum(resMuc(kStep).BM) > 0
                % scale the optimal solution
                XmucCur = resMuc(kStep).BM * optionsMuc.BMweight  / sum(resMuc(kStep).BM);
                fluxMuc(:,kStep) = resMuc(kStep).flux * optionsMuc.BMweight / sum(resMuc(kStep).BM);
            else
                % infeasible growth for the mucosal community, which should not happen in realistic gut microbiota.
                % The model parameters for the simulation are probably too restrictive.
                XmucCur = zeros(nSp,1);
                fluxMuc(:,kStep) = 0;
            end
        else
            % no difference in uptake bounds. Use the previous solution
            if kStep == 1
                if kSave == 1
                    % the first save file
                    fprintf('Sim $%d:\n', nSim);
                    resMuc(kStep) = resMuc(kStepPrev);
                    fluxMuc(:,kStep) = fluxMuc(:,kStepPrev);
                else
                    % new save file
                    resMuc(kStep) = resMuc(saveFre);
                    fluxMuc(:,kStep) = fluxMuc(:,saveFre);
                end
            else
                % not a new save file, just use the previous solution
                resMuc(kStep) = resMuc(kStep - 1);
                fluxMuc(:,kStep) = fluxMuc(:,kStep-1);
            end
            if sum(resMuc(kStep).BM) > 0
                % current biomass composition of the mucosal microbiota
                XmucCur = resMuc(kStep).BM * optionsMuc.BMweight  / sum(resMuc(kStep).BM);
            else
                XmucCur = zeros(nSp,1);
            end
        end
        
        if any(isnan(XmucCur .* resMuc(kStep).GRmax * dtLum))
            error('Sim #%d: Error in predicting growth in the mucosal microbiota!', nSim);
        end
        
        % store the current uptake bounds for comparison at the next iteration
        ubMucPrev = ubMucCur;
        % store the mucosal biomass composition
        Xmuc(:,kStep) = XmucCur;
        % store the actual oxygen uptake by the mucosal microbiota
        o2utMuc(2,kStep) = fluxMuc(model.indCom.EXcom(O2,1),kStep) - fluxMuc(model.indCom.EXcom(O2,2),kStep);
        % store the community growth rate
        GRmuc(kStep) = resMuc(kStep).GRmax;
        % concentration change vector. Actuate in each dtLum step
        C_changeByMuc = fluxMuc(model.indCom.EXcom(:,2),kStep) - fluxMuc(model.indCom.EXcom(:,1),kStep);
        C_changeByMuc(O2) = 0;
        % Concentration vector after exchanges by the mucosal microbiota.
        % Use linear change because the biomass of mucosal community is constant
        Cmuc = C + C_changeByMuc * dtMuc;
        Cmuc(Cmuc < tol) = 0;
        
        %% Simulate the luminal microbiota
        kStepLum = 1;
        % store the current luminal microbial biomass composition
        Xlum(:,kStepLum,kStep) = XlumCur;
        % store the concentration vector in each step
        Ct(:,kStepLum,kStep) = C;
        % store the simulation time
        tLum = t;
        timeLum(kStepLum, kStep) = tLum;
        
        while true
            % oxygen uptake bounds
            o2utLum(1,kStepLum,kStep) = o2fluxLum(j);
            % random order for the organisms in running DMMM
            optOrder(:,kStepLum,kStep) = randperm(nSp)';
            % O2 consumed by the luminal community
            o2consume = 0;
            % update concentration vector due to exchanges by the mucosal community
            C = C + C_changeByMuc * dtLum;
            % increase in luminal microbial biomass
            dXlum = zeros(nSp,1);
            % Run DMMM
            for jSpCt = 1:nSp
                jSp = optOrder(jSpCt,kStepLum,kStep);
                if XlumCur(jSp) > 0
                    % shut down other members
                    modelJ = killSpecies(model, setdiff(1:nSp, jSp));
                    % uptake bounds
                    modelJ.ub(modelJ.indCom.EXcom(:, 1)) = min([10000 * ones(nCom, 1), C / dtLum, Cmuc / dtLum], [], 2);
                    modelJ.ub(model.indCom.EXcom(O2, 1)) = o2fluxLum(j);
                    modelJ.ub(modelJ.indCom.EXcom(modelJ.ub(modelJ.indCom.EXcom(:, 1)) < feasTol, 1)) = 0;
                    [~, resLum(jSp, kStepLum, kStep)] = optimizeCbModelFixBm(modelJ, XlumCur, 1, solverParam);
                    if strcmp(resLum(jSp,kStepLum,kStep).stat,'infeasible')
                        dying(jSp,kStepLum,kStep) = true;
                        % if infeasible, find maximum biomass sustainable.
                        optionsLumSp = optionsLum;
                        optionsLumSp.GR = 0;
                        optionsLumSp.BMcon = sparse(1:nSp,1:nSp,ones(nSp,1),nSp,nSp);
                        optionsLumSp.BMcon(jSp,:) = [];
                        optionsLumSp.BMrhs = zeros(nSp-1,1);
                        optionsLumSp.BMcsense = char('E'*ones(1,nSp-1));
                        [~, resFixGr] = optimizeCbModelFixGr(modelJ, optionsLumSp, solverParam);
                        resFixGr.GRmax = zeros(nSp, 1);
                        resLum(jSp, kStepLum, kStep) = resFixGr;
                        if resLum(jSp, kStepLum, kStep).BM(jSp) > XlumCur(jSp)
                            resLum(jSp, kStepLum, kStep).flux = resLum(jSp, kStepLum, kStep).flux ...
                                * XlumCur(jSp) / resLum(jSp,kStepLum,kStep).BM(jSp);
                            resLum(jSp,kStepLum,kStep).Ut = resLum(jSp,kStepLum,kStep).Ut ...
                                * XlumCur(jSp) / resLum(jSp,kStepLum,kStep).BM(jSp);
                            resLum(jSp,kStepLum,kStep).Ex = resLum(jSp,kStepLum,kStep).Ex ...
                                * XlumCur(jSp) / resLum(jSp,kStepLum,kStep).BM(jSp);
                            resLum(jSp,kStepLum,kStep).BM(jSp) = XlumCur(jSp);
                        end
                        dXlum(jSp) = resLum(jSp,kStepLum,kStep).BM(jSp) - XlumCur(jSp);
                        if dXlum(jSp) > 1e-8
                            error('Sim #%d: Current amount of luminal biomass should be unsustainable.', nSim);
                        end                        
                    else
                        %feasible case
                        dXlum(jSp) = resLum(jSp,kStepLum,kStep).vBM(jSp) * dtLum;
                    end
                    
                    if any(isnan(resLum(jSp,kStepLum,kStep).vBM))
                        error('Sim #%d: Error in predicting growth in the luminal microbiota!', nSim);
                    end
                    % update concentration vector due to exchanges by the luminal microbiota
                    C = C + (resLum(jSp,kStepLum,kStep).Ex - resLum(jSp,kStepLum,kStep).Ut) * dtLum;
                    C(O2) = 10000;
                    % ensure non-negative values with minimum threshold to avoid too
                    % many small steps caused by them
                    C(C < tol) = 0;
                    o2consume = o2consume + resLum(jSp,kStepLum,kStep).Ut(O2) - resLum(jSp,kStepLum,kStep).Ex(O2);
                    GRlum(jSp, kStepLum, kStep) = resLum(jSp,kStepLum,kStep).GRmax(jSp);
                else
                    resLum(jSp,kStepLum,kStep).Ex = zeros(nCom,1);
                    resLum(jSp,kStepLum,kStep).Ut = zeros(nCom,1);
                    resLum(jSp,kStepLum,kStep).GRmax = zeros(nSp,1);
                    resLum(jSp,kStepLum,kStep).BM = XlumCur;
                    resLum(jSp,kStepLum,kStep).vBM = zeros(nSp,1);
                    GRlum(jSp, kStepLum, kStep) = 0;
                end
            end
            % store the actual oxygen uptake
            o2utLum(2,kStepLum,kStep) = o2consume;
            
            % increase in luminal biomass = growth + detechment from mucosal community
            % assuming constant biomass of mucosal community in which detech = growth
            XlumCur = XlumCur + XmucCur * resMuc(kStep).GRmax * dtLum + dXlum;
            XlumCur(XlumCur < 0) = 0;
            
            % update time and step
            if kStepLum >= nStepLum
                break
            end
            tLum = tLum + dtLum;
            kStepLum = kStepLum + 1;
            
            % store time, concentration vector and luminal biomass composition
            timeLum(kStepLum,kStep) = tLum;
            Ct(:,kStepLum,kStep) = C;
            Xlum(:,kStepLum,kStep) = XlumCur;
            
        end
        
        % print after each large step
        t = t + dtMuc;
        fprintf('Sim #%d. Section %d. Step %d. time %.1f\t%04d-%02d-%02d %02d:%02d:%02.0f\n', nSim, j, kTotal, t, clock);
        nextJ = t >= maxT;
        if nextJ && j == nSect
            finish = true;
        end
        if kStep == saveFre || nextJ
            % save if new section or new save file
            [j0, kSave0, kTotal0, kStep0] = deal(j, kSave, kTotal, kStep);
            if saveDetail
                save(sprintf(['%s_sect%dsave%0' num2str(digit) 'd.mat'], saveName, j, kSave), ...
                    'time', 'scSolve', 'scFinish', 'resMuc','fluxMuc', 'Xmuc', 'o2utMuc', ...  % muc-level variables
                    'timeLum', 'Xlum', 'Ct', 'resLum', 'optOrder', 'XlumCur', 'o2utLum', ...  % lum-level variables
                    't', 'nextJ', 'C', 'finish', ...  % variables across sections
                    'j0', 'kSave0', 'kStep0', 'kTotal0'); %counter
            else
                save(sprintf(['%s_sect%dsave%0' num2str(digit) 'd.mat'], saveName, j, kSave), ...
                    'time', 'scSolve', 'scFinish', 'Xmuc', 'o2utMuc', 'GRmuc', ...  % muc-level variables
                    'timeLum', 'Xlum', 'optOrder', 'XlumCur', 'o2utLum', 'GRlum', ...  % lum-level variables
                    't', 'nextJ', 'C', 'finish', ...  % variables across sections
                    'j0', 'kSave0', 'kStep0', 'kTotal0'); %counter
            end
            fprintf('Sim #%d. Section #%d. Save #%d.\t%04d-%02d-%02d %02d:%02d:%02.0f\n', nSim, j, kSave, clock);
            kSave = kSave + 1;
            kStep = 0;
        end
        
        kTotal = kTotal + 1;
        if nextJ
            kStepPrev = kStep;
            break
        end
        kStep = kStep + 1;
        time(kStep) = t;
    end
end
toc(tReal);
end

function [sol, result] = optimizeCbModelFixBm(modelCom, BMfix, minNorm, solverParam)
% [sol, result] = optimizeCbModelFixBm(modelCom, BMfix, minNorm, solverParam)
% maximize the sum of indiviudal growth rates given the biomass (or relative abundance) of each individual organism.
%
% INPUT:
%    modelCom:       A community COBRA model structure with the following fields (created using `createMultipleSpeciesModel`)
%                    (the first 5 fields are required, at least one of the last two is needed. Can be obtained using `getMultiSpecisModelId`):
%
%                      * S - Stoichiometric matrix
%                      * b - Right hand side
%                      * c - Objective coefficients
%                      * lb - Lower bounds
%                      * ub - Upper bounds
%                      * infoCom - structure containing community reaction info
%                      * indCom - the index structure corresponding to `infoCom`
%
%    BMfix           biomass amount/relative abundance for each organism
%    minNorm         minimize the sum of absolute flux or not, default 0
%
% OUTPUT
%    sol:            solveCobraLP solution structure
%    result:         structure with the following fields:
%                      * GRmax - maximum specific growth rate found (/h)
%                      * vBM - biomass formation rate (gdw/h)
%                      * BM - Biomass vector at GRmax (gdw)
%                      * Ut - uptake fluxes (mmol/h)
%                      * Ex - export fluxes (mmol/h)
%                      * flux - flux distribution for the original model
%    stat:           status at the termination of the algorithm
%                      * 'infeasible' - infeasible model (e.g., not satisfying maintenance requirement)
%                      * 'optimal' - optimal growth rate found

%% Initialization
nSp = numel(modelCom.infoCom.spAbbr); %number of species
if ~exist('BMfix', 'var')
    BMfix = ones(nSp, 1);
end
if ~exist('minNorm', 'var')
    minNorm = 0;
end
if ~exist('solverParam', 'var')
    solverParam = struct();
end
[m, n] = size(modelCom.S);

%% Construct LP 

[feasTol, ~] = getCobraSolverParams('LP',{'feasTol'; 'optTol'}, solverParam);

if ~isfield(modelCom, 'csense')
    cs = char('E' * ones(1, m));
else 
    cs = modelCom.csense(:)';
end
LP = struct();
LP.A = modelCom.S;
LP.b = modelCom.b;
LP.c = modelCom.c;
lbBM = modelCom.lb;
ubBM = modelCom.ub;
for jSp = 1:nSp
    % LB^k_j * X^k <= V^k_j
    lbBM(modelCom.indCom.rxnSps == jSp) = lbBM(modelCom.indCom.rxnSps == jSp) * BMfix(jSp);
    % UB^k_j * X^k >= V^k_j
    ubBM(modelCom.indCom.rxnSps == jSp) = ubBM(modelCom.indCom.rxnSps == jSp) * BMfix(jSp);
end
LP.lb = lbBM;
LP.ub = ubBM;
LP.osense = -1;
LP.csense = cs;

sol = solveCobraLP(LP, solverParam);

result = struct();
[result.GRmax, result.vBM, result.BM, result.Ut, result.Ex, result.flux, ...
    result.iter0, result.iter, result.stat] = deal([]);

% check the feasibility of the solution manually
if sol.stat ~= 1
    result.stat = 'infeasible';
    return
end
GRmax = sol.full(modelCom.indCom.spBm);
GRmax(BMfix > 0) = GRmax(BMfix > 0) ./ BMfix(BMfix > 0);
result.GRmax = GRmax;
result.vBM = sol.full(modelCom.indCom.spBm);
result.BM = BMfix;
result.Ut = sol.full(modelCom.indCom.EXcom(:,1));
result.Ex = sol.full(modelCom.indCom.EXcom(:,2));
result.flux = sol.full(1:n);
result.iter0 = [];
result.iter = [];
result.stat = 'optimal';
  
if isscalar(minNorm) && minNorm == 1
    sol.full(modelCom.indCom.spBm(sol.full(modelCom.indCom.spBm) < 0)) = 0;
    fval = LP.c(:)' * sol.full;
    LP.A = [LP.A,                      sparse(m, n); ...  % Sv = 0
            LP.c(:)',                  sparse(1, n); ...  % c'v >= fval
            sparse(1:n, 1:n, 1, n, n), sparse(1:n, 1:n, -1, n, n); ...  % v <= |v|
            sparse(1:n, 1:n, -1, n, n), sparse(1:n, 1:n, -1, n, n)];  % -v <= |v|
    LP.b = [LP.b; fval * (1 - feasTol * 10); zeros(2 * n, 1)];
    LP.c = zeros(n * 2, 1);
    LP.lb = [LP.lb; zeros(n, 1)];
    LP.ub = [LP.ub; 1000 * ones(n ,1)];
    LP.csense = [LP.csense(:)', 'G', repmat('L', 1, 2 * n)];
    LP.osense = 1;
    sol = solveCobraLP(LP, solverParam);
    if sol.stat ~= 1
        result.stat = 'optimal but unable to minimize total absolute flux';
    else
        result.GRmax = GRmax;
        result.vBM = sol.full(modelCom.indCom.spBm);
        result.BM = BMfix;
        result.Ut = sol.full(modelCom.indCom.EXcom(:,1));
        result.Ex = sol.full(modelCom.indCom.EXcom(:,2));
        result.flux = sol.full(1:n);
    end
end
end

function [sol, result] = optimizeCbModelFixGr(modelCom, options, solverParam)
% [sol, result] = optimizeCbModelFixGr(modelCom, options, solverParam)
% Find maximum sum of biomass of the community given a steady-state growth rate
%
% INPUT:
%    modelCom:       A community COBRA model structure with the following fields (created using `createMultipleSpeciesModel`)
%                    (the first 5 fields are required, at least one of the last two is needed. Can be obtained using `getMultiSpecisModelId`):
%
%                      * S - Stoichiometric matrix
%                      * b - Right hand side
%                      * c - Objective coefficients
%                      * lb - Lower bounds
%                      * ub - Upper bounds
%                      * infoCom - structure containing community reaction info
%                      * indCom - the index structure corresponding to `infoCom`
%
% options (optional)   struct with the following possible fields:
%     GR              Fixed growth rate for the community
%     BMcon           Biomass constraint matrix (sum(a_ij * X_j) </=/> b_i)
%                      (given as K x N_species matrix for K constraints)
%                      e.g. [0 1 1 0] for X_2 + X_3 in a 4-species model
%     BMrhs           RHS for BMcon, K x 1 vector for K constraints
%     BMcsense        Sense of the constraint, 'L', 'E', 'G' for <=, =, >=
%     minNorm         Same as minNorm in optimizeCbModel, default 0
%
% OUTPUT
%    sol:            solveCobraLP solution structure
%    result:         structure with the following fields:
%                      * GRmax - maximum specific growth rate found (/h)
%                      * vBM - biomass formation rate (gdw/h)
%                      * BM - Biomass vector at GRmax (gdw)
%                      * Ut - uptake fluxes (mmol/h)
%                      * Ex - export fluxes (mmol/h)
%                      * flux - flux distribution for the original model
%    stat:           status at the termination of the algorithm
%                      * 'optimal' - optimal growth rate found
%                      * 'infeasible' - infeasible model (should not happen)

%% Initialization

%get paramters
if ~exist('options', 'var')
    options = struct();
end
if ~exist('solverParam', 'var')
    solverParam = struct();
end
[GR, BMcon, BMrhs, BMcsense, BMgdw, minNorm] = getSpatialGutParams(...
    {'GR', 'BMcon', 'BMrhs', 'BMcsense', 'BMgdw', 'minNorm'}, options, modelCom);

[m, n] = size(modelCom.S);
nRxnSp = sum(modelCom.indCom.rxnSps > 0); %number of species-specific rxns
nSp = numel(modelCom.infoCom.spAbbr); %number of species

%% Construct LP 
% upper bound matrix
S_ub = sparse([1:nRxnSp 1:nRxnSp]', [(1:nRxnSp)'; n + modelCom.indCom.rxnSps(1:nRxnSp)],...
    [ones(nRxnSp,1); -modelCom.ub(1:nRxnSp)], nRxnSp, n + nSp);
% lower bound matrix
S_lb = sparse([1:nRxnSp 1:nRxnSp]', [(1:nRxnSp)'; n + modelCom.indCom.rxnSps(1:nRxnSp)],...
    [-ones(nRxnSp,1); modelCom.lb(1:nRxnSp)], nRxnSp, n + nSp);
%growth rate and biomass link matrix
grSp = GR * ones(nSp, 1);
S_gr = sparse([1:nSp 1:nSp]', [modelCom.indCom.spBm(:) (n + 1:n + nSp)'],...
    [BMgdw(:); -grSp], nSp, n + nSp);
if isempty(BMcon)
    A = [modelCom.S sparse([],[],[], m, nSp); S_ub; S_lb; S_gr];
else
    A = [modelCom.S sparse([],[],[], m, nSp); S_ub; S_lb; S_gr;...
        sparse([],[],[],size(BMcon, 1), n) BMcon];
end
LP = struct();
LP.A = A;
LP.b = [modelCom.b; zeros(2 * nRxnSp + nSp, 1)]; 
if ~isempty(BMcon)
    LP.b = [LP.b; BMrhs];
end
LP.c = zeros(n + nSp, 1);
LP.c(modelCom.indCom.spBm) = 1;
% organism-specific fluxes bounded by biomass variable but not by constant
LP.lb = -inf(nRxnSp, 1);
LP.lb(modelCom.lb(1:nRxnSp) >= 0) = 0;
LP.lb = [LP.lb; modelCom.lb((nRxnSp + 1): n); zeros(nSp, 1)];
% biomass upper bound should also be arbitrarily large, but set as 1000 here
LP.ub = inf(nRxnSp, 1);
LP.ub(modelCom.ub(1:nRxnSp)<=0) = 0;
LP.ub = [LP.ub; modelCom.ub((nRxnSp + 1): n); 1000 * ones(nSp, 1)];
%handle constraint sense
if ~isfield(modelCom, 'csense')
    LP.csense = char(['E' * ones(1, m) 'L' * ones(1, 2 * nRxnSp) 'E' * ones(1, nSp) BMcsense(:)']);
else
    LP.csense = [modelCom.csense(:)' char(['L' * ones(1, 2 * nRxnSp) 'E' * ones(1, nSp) BMcsense(:)'])];
end
LP.osense = -1;

[feasTol, ~] = getCobraSolverParams('LP',{'feasTol'; 'optTol'}, solverParam);

result = struct();
[result.GRmax, result.vBM, result.BM, result.Ut, result.Ex, ...
    result.flux, result.iter0, result.iter] = deal([]);

sol = solveCobraLP(LP, solverParam);
if sol.stat ~= 1
    result.stat = 'infeasible';
    return
end
result.GRmax = GR;
result.vBM = sol.full(modelCom.indCom.spBm);
result.BM = sol.full((n + 1):(n + nSp));
result.Ut = sol.full(modelCom.indCom.EXcom(:,1));
result.Ex = sol.full(modelCom.indCom.EXcom(:,2));
result.flux = sol.full(1:n);
result.iter0 = [];
result.iter = [];
result.stat = 'optimal';
  
if isscalar(minNorm) && minNorm == 1
    sol.full(modelCom.indCom.spBm(sol.full(modelCom.indCom.spBm) < 0)) = 0;
    fval = LP.c(:)' * sol.full;
    LP.A = [LP.A,                             sparse(m, n); ...  % original constraints
            LP.c(:)',                         sparse(1, n); ...  % c'v >= fval
            sparse(1:n, 1:n, 1, n, n + nSp),  sparse(1:n, 1:n, -1, n, n); ...  % v <= |v|
            sparse(1:n, 1:n, -1, n, n + nSp), sparse(1:n, 1:n, -1, n, n)];  % -v <= |v|
    LP.b = [LP.b; fval * (1 - feasTol * 10); zeros(2 * n, 1)];
    LP.c(:) = 0;
    LP.lb = [LP.lb; zeros(n, 1)];
    LP.ub = [LP.ub; 1000 * ones(n ,1)];
    LP.csense = [LP.csense(:)', 'G', repmat('L', 1, 2 * n)];
    LP.osense = 1;
    sol = solveCobraLP(LP, solverParam);
    if sol.stat ~= 1
        result.stat = 'optimal but unable to minimize total absolute flux';
    else
        result.GRmax = GRmax;
        result.vBM = sol.full(modelCom.indCom.spBm);
        result.BM = BMfix;
        result.Ut = sol.full(modelCom.indCom.EXcom(:,1));
        result.Ex = sol.full(modelCom.indCom.EXcom(:,2));
        result.flux = sol.full(1:n);
    end
end

end

function varargout = getSpatialGutParams(param2get, options, modelCom)
% get the required default parameters
% [param_1, ..., param_N] = getSpatialGutParams({'param_1',...,'param_N'},options,modelCom)
%
% INPUT:
%   'param_1',...,'param_N': parameter names
%   options: option structure
%            If the required parameter is a field in options, take from
%            options. Otherwise, return the default value.
%   modelCom: the community model for which parameters are constructed.
   
if nargin < 3
    modelCom = struct('rxns',[]);
    modelCom.infoCom.spAbbr = {};
    modelCom.infoCom.rxnSps = {};
end
if nargin < 2 || isempty(options)
    options = struct();
end
if ischar(param2get)
    param2get = {param2get};
end

varargout = cell(numel(param2get), 1);
for j = 1:numel(param2get)
    if isfield(options, param2get{j})
        %if provided in the call
        varargout{j} = options.(param2get{j});
    else
        %use default if default exist and not provided
        %return empty if no default
        varargout{j} = paramDefault(param2get{j}, modelCom);
    end
end
end

function param = paramDefault(paramName,modelCom)
% get default parameter values
switch paramName

    % parameters for spatialGut
    case 'resultTmp',   param = struct('GRmax',[],'vBM',[],'BM',[],'Ut',[],...
                                'Ex',[],'flux',[],'iter0',[],'iter',[],'stat','');  % result template
    case 'dtMuc',       param = 0.5;
    case 'dtLum',       param = 1 / 12;
    case 'saveName',    param = ['spatialGutSim' filesep ' test'];
    case 'saveFre',     param = 21;
    case 'o2Id',        param = 'o2[u]';
    case 'nSim',        param = 1;
    case 'saveDetail',  param = false;
        
    %parameters for optimizeCbModelFixGr
    case 'GR',          param = 0; % growth rate fixed
    case 'BMgdw',       param = ones(numel(modelCom.infoCom.spAbbr), 1);  % relative molecular weight of biomass. For scaling the relative abundance
    case 'minNorm',     param = 0;
    
    otherwise,          param = [];
end
end

