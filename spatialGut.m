function finish = spatialGut(model, options, solverParam)
% finish = spatialGut(model, options, solverParam)
%***Simplify the diffusion process into simple oxygen flux proportional to
%surface area and diffusion coefficients
%Steady state for mucosal community and meanwhile dynamic change for
%luminal community. (SteadyCom for mucosal community, dSteadyCom for 
%luminal community and oxygen diffusion have all different time scales)
%Input:
%  model:    COBRA community model (template for all communities)
%  options: structure with the following fields:
%     C0:  mCom x 2 x nCond array. c_ijk is the amount of metabolites (mmol) 
%         for the i-th community metabolites in the mucsoal (j=1) or
%         luminal (j=2) community under the k-th condition representing 
%         a section of the intestines.
%     X:  nCond x 1 vector of mucosal community biomass density (gdw/L) for various
%         sections of the intestines. (e.g. 0.2 gdw/L)
%     D:  nCond x 1 vector of the dilution rate of the mucosal community
%         for various sections of the intestines (should proportional to
%         the flow rate of luminal content/chyme). Either X or D must be
%         supplied. Use X whenever X is suppled.
%     T:  nCond x 1 vector of transit time for various sections of the intestines
%     dt: small time step for dynamic simulation (1//60);
%    Rmuc: thinkness of the mucosal community (for calculating
%           cross-section area) (mm)
%     R: radius of the intestines (for calculating volume) (mm)
%     Len: length of the intestines (for calculating volume) (mm)
%     L:   number of discretization step for calculating oxygen diffusion
%        (oxygen available to each community calculated by summing the concentration
%         multiplied by the volume occupied by the community over each region
%         in the discretization scheme)
%     saveDsteady: filename for saving
%     saveFre: steps per savefile
%     O2Id: oxygen community metabolite ID (o2[u])
%     O2d: nCond x 1 vector of oxygen diffusion coefficients between mucosal and
%          luminal community (mm^2/hr)
%     pO2: partial pressure of oxygen on the mucosa (which determines
%          oxygen uptake rate) (atm)
%     Kpc: partial pressure / molar concentration for O2 (619 atm / (mol/L))
%     p2c: concentration of O2 per partial pressure (mmol/L / atm)
%          (default 39.3)
%     O2ut: oxygen uptake rate per partial pressure (assume propotional)
%          (nSp x 1 vector)
%     other fields for optimizeCbModelComCplexFixGr (optional)
%
%  solverParam: Cplex solver parameter structure
%
%Given the uptake rate in the k-th section, the mucosal community with
%maximum biomass under a constant flow rate F(k) (as the diltuion rate) is
%determined by SteadyCom. The mucosal community is thus modeled at a community 
%steady state. During the time period T(k). The microbial mass from the mucosal
%community is constantly dilute and becomes the mass in the luminal
%community. The microbes in the luminal community deviate from the
%community steady state since the uptake rate (oxygen availability) has
%changed and evolve to a new steady state over time. The final luminal 
%community in the k-th section becomes the initial community in the k+1-th
%sections. This simulation assumes that the intestines is stable (like in
%mature adults) such that colonization by new microbes into the mucosal
%community is negligible.

%% Check inputs
debug = 0;
finish = true;
tReal = tic;
tol = 1e-10;        
param2get = {'resultTmp','dtMuc','dtLum','dtOxy','C0','T','D','X','O2mucLv','Len',...
    'O2Id','Rmuc','RmucCom','R','L','O2s','O2d','saveName','saveFre', ...
    'O2ut', 'O2total','O2fluxMuc','O2fluxLum',...
    'Kpc','pO2','p2c'};
eval(sprintf('[%s] = getCobraComParams(param2get, options, model);', ...
            strjoin(param2get, ',')...
            )...
    );
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
field = {'C0','T','X'};
for j = 1:numel(field)
    if ~isfield(options,field{j})
        error('%s must be provided in the option structure.',field{j});
    end
end
% number of community metabolites
nCom = size(model.infoCom.Mcom, 1);
if numel(C0) ~= nCom
    error('Size of options.C0 (%d) not equal to the number of community metabolites (%d).', numel(C0), nCom);
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
O2 = find(strcmp(model.infoCom.Mcom, O2Id));
% number of organisms in the community
nSp = numel(model.sps);
% solution structure for mucosal and luminal microbiota
resMuc = repmat(resultTmp,1,saveFre);
resLum = repmat(resultTmp,nSp,nStepLum,saveFre);
% O2 flux bounds and values
o2utMuc = zeros(2,saveFre);
o2utLum = zeros(2,nStepLum,saveFre);
% flux from resMuc may need rescaling after optimization. Record separately.
fluxMuc = zeros(size(model.S,2),saveFre);
% time at which the result was got
[time,o2Muc] = deal(zeros(saveFre,1));
timeLum = zeros(nStepLum,saveFre); 
% metabolite level
[Ct] = deal(zeros(size(C0,1), nStepLum, saveFre));

% biomass (relative abundance)
Xmuc = zeros(nSp, saveFre);
[Xlum,optOrder] = deal(zeros(nSp, nStepLum, saveFre));
%whether the species is dying
dying = false(nSp, nStepLum, saveFre);

%%%%%%%%%%%%%
% Need to recalculate SteadyCom or not
SC = false(saveFre,1);
%Volume of discretization region in each section
Volume = zeros(L, nSect);
%width of each region 
dr = zeros(nSect,1);
%volume ratio occupied by either community in each discretization region
%ratioMuc is the volume ratio in a region with constant oxygen pressure
[ratioMucCom, ratioLumCom, ratioMuc] = deal(zeros(L,nSect));
[O2sMuc,O2sLum] = deal(zeros(nSect,1));


%Variabile for the current step
%nothing in the lumen initially
XlumCur = zeros(nSp,1);

%oxygen concentration profile for calculating diffusion (mmol per unit area)
%assume uniform distribution at the beginning
o2ProfileCur = zeros(L,1);

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

%calculate the physical parameters first
for j = 1:nSect
    %discretization step size for oxygen diffusion
    dr(j) = R(j) / L;
    Volume(:,j) = 2*pi*dr(j)*((dr(j):dr(j):R(j))-(dr(j)/2))';
    %for calculating total amount of oxygen available to each community
    ratioMucCom(L:-1:(L-ceil(RmucCom(j)/dr(j))+1),j) = 1;
    ratioLumCom(1:ceil((R(j)-RmucCom(j))/dr(j)),j) = 1;
    %for identifying the region with constant oxygen pressure
    ratioMuc(dr(j):dr(j):R(j) >= R(j)-Rmuc(j)+dr(j),j) = 1;
    if mod(RmucCom(j),dr(j)) > 0
        areaOut = (((L-ceil(RmucCom(j)/dr(j))+1)*dr(j))^2)*pi - ((R(j)-RmucCom(j))^2)*pi;
        ratioMucCom(L-ceil(RmucCom(j)/dr(j))+1,j) = areaOut/Volume(L-ceil(RmucCom(j)/dr(j))+1,j);
    end
    if mod(R(j)-RmucCom(j),dr(j)) > 0
        areaIn = (R(j)-RmucCom(j))^2*pi - (floor((R(j)-RmucCom(j))/dr(j))*dr(j))^2*pi;
        ratioLumCom(ceil((R(j)-RmucCom(j))/dr(j)),j) = areaIn / Volume(ceil((R(j)-RmucCom(j))/dr(j)),j);
    end
    %for calculating diffusion from epithelial cells to region R(L) 
    O2sMuc(j) = O2s(j) * ratioMucCom(L,j);
    O2sLum(j) = O2s(j) * ratioLumCom(L,j);
end

j0 = 0;
[kTotal0,kStep0] = deal(1);
kSave0 = 0;
while true
    if exist(sprintf(['%s_sect%dsave%0' num2str(digit) 'd.mat'], ...
            saveDsteady, j0, kSave0+1), 'file')
        kSave0 = kSave0 + 1;
    elseif exist(sprintf(['%s_sect%dsave%0' num2str(digit) 'd.mat'], ...
            saveDsteady, j0+1, 1), 'file')
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
                saveDsteady, j0, kSave0),'t','XlumCur','C', ...
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
            j0 = j0 + 1;
            kSave0 = 1;
        else
            kSave0 = kSave0 + 1;
        end
        dataPrev = load(sprintf('%s_pre.mat', saveDsteady),'SpRate','SpRateUt','SpRateEx');
        Id = findRxnIDs(model,dataPrev.SpRate);
        model.lb(Id) = dataPrev.SpRateUt;
        model.ub(Id) = dataPrev.SpRateEx;
    end
end
if j0 == 0
    j0 = 1;
    kSave0 = 1;
    nextJ = true;
    SpRate = model.rxns(model.EXsp(model.EXsp~=0));
    SpRateUt = model.lb(model.EXsp(model.EXsp~=0));
    SpRateEx = model.ub(model.EXsp(model.EXsp~=0));
    save(sprintf('%s_pre.mat', saveDsteady), 'options','solverParam',...
        'SpRate','SpRateUt','SpRateEx','nStepLum','nStepOxy','dtMuc','dtLum',...
        'digit','Nstep','feasTol','tol','exFact', ...
        'Volume','dr','ratioMuc','ratioMucCom','ratioLumCom','O2sMuc','O2sLum','R','Rmuc','L',...
        'O2ut','RmucCom','pO2','Kpc','p2c','Len','O2fluxMuc','O2fluxLum');
    clear SpRate SpRateUt SpRateEx
end
init = true;
fprintf('Start from section %d, save %d...\n',j0,kSave0);
for j = j0:nSect
    %mucosal and luminal o2 level defined by oxygen concentration
    if nextJ
        if 0
            if j == 1
                o2ProfileCur(ratioMuc(:,j) > 0) = O2mucLv(j);
            else
                %O2 from mucosal layer
                o2ProfileCurMuc = zeros(L,1);
                o2ProfileCurMuc(ratioMuc(:,j) > 0) = O2mucLv(j);
                %O2 from lumen.
                %total O2 amount in each region
                weight = o2ProfileCur .* Volume(:,j-1) .* ratioLum(:,j-1);
                %Ensure continuity of total oxygen in the lumen. Assume the
                %same amount in the new corresponding region (same l = 1, 2, ..., L)
                o2ProfileCurLum = weight ./ Volume(:,j);
                %if the boundary between mucosal layer and the lumen lies in
                %a discretization region for the new section. Ignore the
                %contribution of O2 from the lumen for simplicity because that
                %region has been assumed to have constant O2 level
                o2ProfileCurLum(ratioMuc(:,j) > 0) = 0;
                o2ProfileCur = o2ProfileCurMuc + o2ProfileCurLum;
            end
        else
            %simply reset oxygen level
            o2ProfileCur = zeros(L,1);
            o2ProfileCur(ratioMuc(:,j) > 0) = pO2(j);
        end
    end
    %average partial pressure
    o2MucCur = sum(o2ProfileCur .* Volume(:,j) .* ratioMucCom(:,j)) / sum(Volume(:,j) .* ratioMucCom(:,j));
    o2LumCur = sum(o2ProfileCur .* Volume(:,j) .* ratioLumCom(:,j)) / sum(Volume(:,j) .* ratioLumCom(:,j));
    %counter for step and number of save files
    kStep = 1;
    if init
        kSave = kSave0;
        init = false;
    else
        kSave = 1;
    end
    %[kStep,kSave] = deal(1);
    
    %record the initial values
    [time(kStep),o2Muc(kStep)] = deal(t,o2MucCur);

    optionsMuc.BMweight = X(j) * sum(Volume(:,j).*ratioMucCom(:,j)) * Len(j) /(100^3) ;
    %     optionsMuc.BMcon = ones(1,nSp);
    %     optionsMuc.BMrhs = 2 * X(j);
    %     optionsMuc.BMcsense = 'L';
    
    %simulation time for the current section
    maxT = sum(T(1:j)) - dtMuc * 1e-5;
    while true
        %availability profile
%         availMuc = C > 0 & metPossUT;
%         availMuc(O2) = o2MucCur > 1e-5 & metPossUT(O2);
        %%%   unlimited community uptake rate. Uptake rate should be bounded at
        %%%   the organism level
        %Maximum specific uptake rate by organism
        o2ubSp = -abs(O2ut * o2MucCur);
        %maximum community uptake rate
        o2ub = abs(O2total(j) * optionsMuc.BMweight * o2MucCur);
        o2ub = min([o2ub, o2MucCur * p2c * sum(Volume(:,j) .* ratioMucCom(:,j)) * Len(j) / dtMuc / 1e6]);
        %o2ubMuc(:,kStep) = [o2ub; o2ubSp(:)];
        o2utMuc(1,kStep) = O2fluxMuc(j); %* optionsMuc.BMweight;
        ubMucCur = min(10000*ones(nCom,1),C./dtMuc);
        %         ubMucCur(O2) = min([10000, o2MucCur/dtMuc]);
        ubMucCur(O2) = 10000;
        ubMucCur(ubMucCur < feasTol) = 0; %may cause numerical problem for such small ub
        ubDiff = abs(ubMucCur-ubMucPrev);
        if ~(isequal(ubMucCur ~= 0, ubMucPrev ~= 0)) || sum(ubDiff(ubMucPrev ~= 0) ./ ubMucPrev(ubMucPrev ~= 0)) > 0.01
            model.ub(model.EXcom(:,1)) = ubMucCur;
            %model.lb(model.EXsp(O2,model.EXsp(O2,:)~=0)) = o2ubSp(model.EXsp(O2,:)~=0);
            %model.ub(model.EXcom(O2,1)) = o2ub;
            model.ub(model.EXcom(O2,1)) = O2fluxMuc(j);% * optionsMuc.BMweight;
            SC(kStep) = true;
            %perform steadyCom
            if debug
                fprintf('optimizeCbModelComCplexNewMC2:\n')
            end
            [~,resMuc(kStep)] = optimizeCbModelComCplexNewMC2(model, optionsMuc, solverParam);
%             [~,resMuc(kStep)] = SteadyCom(model, optionsMuc, solverParam);
            if strcmp(resMuc(kStep).stat,'time limit exceeded') || strcmp(resMuc(kStep).stat,'infeasible')
                finish = false;
                resMuc(kStep).Ex = zeros(nCom,1);
                resMuc(kStep).Ut = zeros(nCom,1);
                resMuc(kStep).GRmax = zeros(nSp,1);
                resMuc(kStep).BM = zeros(nSp,1);
                resMuc(kStep).vBM = zeros(nSp,1);
                resMuc(kStep).flux = zeros(n,1);
            end
            %scale the optimal solution
            if sum(resMuc(kStep).BM) > 0
                XmucCur = resMuc(kStep).BM * optionsMuc.BMweight  / sum(resMuc(kStep).BM);
                fluxMuc(:,kStep) = resMuc(kStep).flux * optionsMuc.BMweight / sum(resMuc(kStep).BM);
            else
                %this means an infeasible growth for the mucosal community,
                %which should not happen in realistic gut microbiota.
                %The model parameters for the simulation are probably too
                %restrictive.
                XmucCur = zeros(nSp,1);
                fluxMuc(:,kStep) = 0;
            end
             
            %         if isfield(optionsJ,'X')
            %             optionsJ.BMweight = X(j);
            %             [~,resMuc(j),LP] = optimizeCbModelComCplexNewMC2(model, optionsJ, solverParam);
            %         else
            %             optionsJ.GR = D(j);
            %             [~,resMuc(j),LP] = optimizeCbModelComCplexFixGr(model,optionsJ, solverParam);
            %         end
        else
            if kStep == 1
                if kSave == 1
                    %must be a new section
                    fprintf('nCR: %d\.%d\n', options.nCR(1), options.nCR(2));
                    resMuc(kStep) = resMuc(kStepPrev);
                    fluxMuc(:,kStep) = fluxMuc(:,kStepPrev);
                else
                    resMuc(kStep) = resMuc(saveFre);
                    fluxMuc(:,kStep) = fluxMuc(:,saveFre);
                end
            else
                resMuc(kStep) = resMuc(kStep - 1);
                fluxMuc(:,kStep) = fluxMuc(:,kStep-1);
            end
            if sum(resMuc(kStep).BM) > 0
                XmucCur = resMuc(kStep).BM * optionsMuc.BMweight  / sum(resMuc(kStep).BM);
            else
                XmucCur = zeros(nSp,1);
            end
        end
        ubMucPrev = ubMucCur;
        Xmuc(:,kStep) = XmucCur;
        o2utMuc(2,kStep) = fluxMuc(model.EXcom(O2,1),kStep) - fluxMuc(model.EXcom(O2,2),kStep);
        %concentration change vector. Actuate in each dtLum step
        C_changeByMuc = fluxMuc(model.EXcom(:,2),kStep) - fluxMuc(model.EXcom(:,1),kStep);
        C_changeByMuc(O2) = 0;
        %ignore oxygen change in mucosal community. Assume neglegible
        %concentration vector appears to the luminal community (but not the
        %actual recorded concentration
        %Use linear change because the biomass of mucosal community is constant
        Cmuc = C + C_changeByMuc * dtMuc;
        Cmuc(Cmuc < tol) = 0;
        
        kStepLum = 1;
        Xlum(:,kStepLum,kStep) = XlumCur;
        Ct(:,kStepLum,kStep) = C;
        tLum = t;
        timeLum(kStepLum, kStep) = tLum;
        
        while true
            %luminal community
            model.ub(model.EXcom(:,1)) = min([10000*ones(nCom,1),C/dtLum, Cmuc/dtLum],[],2);
            %model.ub(model.EXcom(O2,1)) = min([10000, o2LumCur/dtLum]);
            %model.ub(model.EXcom(O2,1)) = 10000;
            %Maximum specific uptake rate by organism
            o2ubSp = -abs(O2ut * o2LumCur);
            %maximum community uptake rate
            o2ub = abs(O2total(j) * sum(XlumCur) * o2LumCur);
            o2ub = min([o2ub, o2LumCur * p2c * sum(Volume(:,j) .* ratioLumCom(:,j)) * Len(j) / dtLum / 1e6]);
            %o2utLum(:,kStepLum,kStep) = [o2ub; o2ubSp(:)];
            o2utLum(1,kStepLum,kStep) = O2fluxLum(j);% * sum(XlumCur);
            %model.lb(model.EXsp(O2,model.EXsp(O2,:)~=0)) = o2ubSp(model.EXsp(O2,:)~=0);
            %model.ub(model.EXcom(O2,1)) = O2fluxLum(j) * sum(XlumCur);
            %model.ub(model.EXcom(model.ub(model.EXcom(:,1)) < feasTol, 1)) = 0;%may cause numerical problem for such small ub
            optionsLum.BMfix = XlumCur;
            optOrder(:,kStepLum,kStep) = randperm(nSp)';
             %O2 consumed by the luminal community
            o2consume = 0;
            %%%%%%%%
            %%% Use DMMM approach
            %metabolite changed by mucosal community 
            C = C + C_changeByMuc * dtLum;
            dXlum = zeros(nSp,1);
            for jSpCt = 1:nSp
                jSp = optOrder(jSpCt,kStepLum,kStep);
                if XlumCur(jSp) > 0
                    modelJ = killSpecies(model, setdiff(1:nSp,jSp));
                    modelJ.ub(modelJ.EXcom(:,1)) = min([10000*ones(nCom,1),C/dtLum, Cmuc/dtLum],[],2);
                    modelJ.ub(model.EXcom(O2,1)) = O2fluxLum(j);% * sum(XlumCur);
                    modelJ.ub(modelJ.EXcom(modelJ.ub(modelJ.EXcom(:,1)) < feasTol, 1)) = 0;
                    if debug
                        fprintf('optimizeCbModelComCplexFixBMwoSS:\n')
                    end
                    [~, resLum(jSp,kStepLum,kStep)] = optimizeCbModelComCplexFixBMwoSS(modelJ,optionsLum, solverParam);
                    if strcmp(resLum(jSp,kStepLum,kStep).stat,'infeasible')
                        dying(jSp,kStepLum,kStep) = true;
                        %if infeasible, find maximum biomass sustainable.
                        %Other assume dead.
                        optionsLumSp = optionsLum;
                        optionsLumSp.GR = 0;
                        optionsLumSp.BMcon = sparse(1:nSp,1:nSp,ones(nSp,1),nSp,nSp);
                        optionsLumSp.BMcon(jSp,:) = [];
                        optionsLumSp.BMrhs = zeros(nSp-1,1);
                        optionsLumSp.BMcsense = char('E'*ones(1,nSp-1));
                        if isfield(options,'nCR')
                            optionsLumSp.nCR = options.nCR;
                        end
                        if debug
                            fprintf('optimizeCbModelComCplexFixGr:\n')
                        end
                        [~, resFixGr] = optimizeCbModelComCplexFixGr(modelJ,optionsLumSp, solverParam);
                        resFixGr.GRmax = zeros(nSp,1);
                        resLum(jSp,kStepLum,kStep) = resFixGr;
                        if resLum(jSp,kStepLum,kStep).BM(jSp) > XlumCur(jSp)
                            resLum(jSp,kStepLum,kStep).flux = resLum(jSp,kStepLum,kStep).flux ...
                                * XlumCur(jSp) / resLum(jSp,kStepLum,kStep).BM(jSp); 
                            resLum(jSp,kStepLum,kStep).Ut = resLum(jSp,kStepLum,kStep).Ut ...
                                * XlumCur(jSp) / resLum(jSp,kStepLum,kStep).BM(jSp);
                            resLum(jSp,kStepLum,kStep).Ex = resLum(jSp,kStepLum,kStep).Ex ...
                                * XlumCur(jSp) / resLum(jSp,kStepLum,kStep).BM(jSp);
                            resLum(jSp,kStepLum,kStep).BM(jSp) = XlumCur(jSp);
                        end
                        dXlum(jSp) = resLum(jSp,kStepLum,kStep).BM(jSp) - XlumCur(jSp);
                        if dXlum(jSp) > 1e-8
                            error('Current amount of biomass should be unsustainable.');
                        end
                        %resLum(jSp,kStepLum,kStep).Ex = zeros(nCom,1);
                        %resLum(jSp,kStepLum,kStep).Ut = zeros(nCom,1);
                        %resLum(jSp,kStepLum,kStep).GRmax = zeros(nSp,1);
                        %resLum(jSp,kStepLum,kStep).BM = XlumCur;
                        %resLum(jSp,kStepLum,kStep).vBM = zeros(nSp,1);
                    else
                        %feasible case
                        dXlum(jSp) = resLum(jSp,kStepLum,kStep).vBM(jSp) * dtLum;
                    end
                    if any(isnan(resLum(jSp,kStepLum,kStep).vBM))
                        if isfield(options,'nCR')
                            error('c:%d, r:%d\tNaN 1!',options.nCR(1),options.nCR(2));
                        else
                            error('NaN!');
                        end
                    end
                    %update metabolite levels
                    C = C + (resLum(jSp,kStepLum,kStep).Ex - resLum(jSp,kStepLum,kStep).Ut) * dtLum;
                    C(O2) = 10000;
                    %Ensure non-negative values with minimum threshold to avoid too
                    %many small steps caused by them
                    C(C < tol) = 0;
                    o2consume = o2consume + resLum(jSp,kStepLum,kStep).Ut(O2) - resLum(jSp,kStepLum,kStep).Ex(O2);
                else
                    resLum(jSp,kStepLum,kStep).Ex = zeros(nCom,1);
                    resLum(jSp,kStepLum,kStep).Ut = zeros(nCom,1);
                    resLum(jSp,kStepLum,kStep).GRmax = zeros(nSp,1);
                    resLum(jSp,kStepLum,kStep).BM = XlumCur;
                    resLum(jSp,kStepLum,kStep).vBM = zeros(nSp,1);
                end
            end
           
            o2utLum(2,kStepLum,kStep) = o2consume;
            %amount of O2 available to the luminal community in each region
            weight = o2ProfileCur .* Volume(:,j) .* ratioLumCom(:,j);
            if any(weight)
                %as a normalized weight with unit sum
                weight = weight/ sum(weight);
            end
            %to approximate the change in O2 concentration caused by
            %consumption of O2 by the luminal community, assuming a
            %decrease in O2 level proportional to the total O2 level in a
            %region
            o2concChangeByLum = o2consume * weight ./ Volume(:,j) / Len(j) * 1e6;
            
            %increase in luminal biomass = growth + detechment from mucosal
            %community (assume constant biomass of mucosal community,
            %detech = growth)
            if any(isnan(XmucCur * resMuc(kStep).GRmax * dtLum))
                if isfield(options,'nCR')
                    error('c:%d, r:%d\tNaN 2!',options.nCR(1),options.nCR(2));
                else
                    error('NaN!');
                end
            end
            XlumCur = XlumCur + XmucCur * resMuc(kStep).GRmax * dtLum + dXlum;
            %for jSp = 1:nSp
            %    XlumCur = XlumCur + resLum(jSp,kStepLum,kStep).vBM * dtLum ;
            %end
            XlumCur(XlumCur < 0) = 0;
        
            %tOxy = tLum;
            %kStepOxy = 1;
            %timeOxy(kStepOxy,kStepLum,kStep) = tOxy;
            %o2Lum(kStepOxy,kStepLum,kStep) = o2LumCur;
            %only record the end points of oxygen diffusion
            %o2Profile(:,kStepLum,kStep) = o2ProfileCur;
            %o2Profile(:,kStepOxy,kStepLum,kStep) = o2ProfileCur;
            %oxygen diffusion must be calculated at much smaller time steps (as
            %we are using forward Euler, [should implement implicit method])
            
            
            if kStepLum >= nStepLum
                break
            end
            tLum = tLum + dtLum;
            kStepLum = kStepLum + 1;
            
            timeLum(kStepLum,kStep) = tLum;
            Ct(:,kStepLum,kStep) = C;
            Xlum(:,kStepLum,kStep) = XlumCur;
        end
        %print and save
        
        t = t + dtMuc;
        
        fprintf('section %d . step %d . time %.1f\t%04d-%02d-%02d %02d:%02d:%02.0f\n',j,kTotal,t,clock);
        nextJ = t >= maxT;
        if kStep == saveFre || nextJ
            [j0,kSave0,kTotal0,kStep0] = deal(j,kSave,kTotal,kStep);
            save(sprintf(['%s_sect%dsave%0' num2str(digit) 'd.mat'], ...
                saveDsteady, j, kSave), ...
                'time','o2Muc','SC','resMuc','fluxMuc', 'Xmuc', ... %muc level variables
                't', 'nextJ',...
                'timeLum','Xlum','Ct','resLum',... %lum level variables
                'XlumCur','C', ... %lum level current variables
                ... %oxy level variables
                'o2utLum', 'o2utMuc', 'optOrder', ...
                 ... %Oxy level current variables
                'j0','kSave0','kStep0','kTotal0'); %counter
            fprintf('section #%d . save #%d.\t%04d-%02d-%02d %02d:%02d:%02.0f\n',j,kSave,clock);
            kSave = kSave + 1;
            kStep = 0;
        end
        
        kTotal = kTotal + 1;
        if nextJ
            kStepPrev = kStep;
            break
        end
        kStep = kStep + 1;
        [time(kStep),o2Muc(kStep)] = deal(t,o2MucCur);
    end
end
toc(tReal);
end


function o2Flux = oxygenDiffusion(o2Profile, R, D)
%o2Profile = oxygenDiffusion(o2Profile, R, D)
%calculate the O2 diffusion flux given the concentration profile
%from 0 to R, diffusion coefficent D
L = numel(o2Profile);
dr = R/L;
% r = ((1:L)' - 0.5) * dr;
o2Flux = zeros(L,1);
for j = 1:L
    if j == 1
        [c1,c2,c3] = deal(o2Profile(j),o2Profile(j),o2Profile(j+1));
    elseif j == L
        [c1,c2,c3] = deal(o2Profile(j-1),o2Profile(j),o2Profile(j));
    else
        [c1,c2,c3] = deal(o2Profile(j-1),o2Profile(j),o2Profile(j+1));
    end
    %     o2Flux(j) = (D/r(j)/(dr^2)) * ((r(j) - (dr/2)) * c1 - (2 * r(j) * c2) ...
    %         + (r(j) + (dr/2)) * c3);
    o2Flux(j) = (D/(dr^2)) * (c1 - 2*c2 + c3 + (c3-c1)/(2*j-1));
end

end





