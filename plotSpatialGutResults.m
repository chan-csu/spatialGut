function [data, dataMean, dataSD, finish, figHandle] = plotSpatialGutResults(saveName, simIDs, expData, expDataSD, T)
% To retrieve and plot the results of spatialGut
%
% USAGE:
%    [data, dataMean, dataSD, finish, figHandle] = plotSpatialGutResults(saveName, simIDs, expData, expDataSD)
%
% OPTIONAL INPUTS:
%    (all optional, default using the partial data presented in the paper)
%    saveName           the same saveName used in calling spatialGut to retrieve the data
%    simIDs             an index vector for the simulations performed to retrieve and plot
%    expData            structure containing the following experimental data:
%                       * AmucAero - nSect-by-1 vector, relative abundance of the aerobes 
%                                    plus facultative anaerobes on the mucus layer
%                       * AmucAnaero - nSect-by-1 vector, relative abundance of the strict
%                                      anaerobes on the mucus layer
%                       * AlumAero - nSect-by-1 vector, relative abundance of the aerobes 
%                                    plus facultative anaerobes in the lumen
%                       * AlumAnaero - nSect-by-1 vector, relative abundance of the strict
%                                      anaerobes in the lumen
%                       * XlumRel - nSect-by-1 vector, the relative microbial biomass level 
%                                   in the lumen along the intestines
%    expDataSD          struture with standard deviations corresponding to the data in expData
%    T                  the retention time for each intestinal section used
%
% OUTPUTS:
%    data               nSim-by-1 structure, each containing the results of one simuation
%    dataMean           structure containing the meam for the data in 'data'
%    dataSD             structure containing the standard deviation for the data in 'data'
%    finish             nSim-by-1 logical vector, indicating whether each simulation is finished
%    figHandle          structure containing the handles for all the figure-related objects for the plotted figure

if nargin == 0 || isempty(saveName)
    saveName = ['spatialGutExample' filesep 'default_params'];
end
if nargin < 2 || isempty(simIDs)
    simIDs = 1:100;
end
if nargin < 3 || isempty(expData)
    d = load('experimentalData.mat');
    expData = d.expData;
    expDataSD = d.expDataSD;
elseif nargin < 4
    expDataSD = expData;
    fn = fieldnames(expData);
    for j = 1:numel(fn)
        expDataSD.(fn{j})(:) = 0;
    end
end
if nargin < 5 || isempty(T)
    T0 = [2 2 2 3 3 3 3];
    T = zeros(1,numel(T0));
    for j = 1:numel(T)
        T(j) = sum(T0(1:j));
    end
    T = [0 T];
end

% gather the data from saved files
[data, dataMean, dataSD, finish] = getSpatialGutSimData(saveName, simIDs);
dataMean.XlumSum0 = dataMean.XlumSum;
% reference level for calculating the relative change
SIPgm = exp(mean(log(dataMean.XlumSum(dataMean.timeLum(:) < 2 & dataMean.XlumSum(:) > 0))));
dataMean.XlumSum = dataMean.XlumSum / SIPgm;
dataSD.XlumSum = dataSD.XlumSum / SIPgm;

% plot the results
figHandle = plotSimVsExp(dataMean, dataSD, expData, expDataSD, T);

end

function [data,dataMean,dataSD,finish] = getSpatialGutSimData(saveName, simIDs, aerobeIDs)

if nargin < 3
    aerobeIDs = 4:5;
end
field2get = {'time','timeLum','Xmuc','Xlum','resMuc','resLum','C','fluxMuc','GRmuc','GRlum','Alum','Amuc'};
dataTemplate = struct();
for k = 1:numel(field2get)
    dataTemplate.(field2get{k}) = [];
end
data = repmat(dataTemplate, numel(simIDs), 1);
finish = false(numel(simIDs), 1);
% collect the data
for i = 1:numel(simIDs)
  saveNameI = sprintf('%s_sim%d', saveName, simIDs(i));
  dataEachSim = dataTemplate;
  if exist(sprintf('%s_pre.mat', saveNameI), 'file')
      d0 = load(sprintf('%s_pre.mat', saveNameI));
      j = 1;
      k = 1;
      if exist(sprintf(['%s_sect%dsave%0' num2str(d0.digit) 'd.mat'], saveNameI, j, k), 'file')
          d = load(sprintf(['%s_sect%dsave%0' num2str(d0.digit) 'd.mat'], saveNameI,k,j));
          ct = 0;
          while true
              ct = ct + 1;
              dataEachSim.time = [dataEachSim.time; d.time(1:d.kStep0)];
              dataEachSim.timeLum = [dataEachSim.timeLum; reshape(d.timeLum(:, 1:d.kStep0), size(d.timeLum, 1) * d.kStep0, 1)];
              dataEachSim.Xmuc = [dataEachSim.Xmuc d.Xmuc(:, 1:d.kStep0)];
              s = size(d.Xlum);
              dataEachSim.Xlum = cat(2, dataEachSim.Xlum, reshape(d.Xlum(:, :, 1:d.kStep0), s(1), s(2) * d.kStep0));
              if isfield(d, 'resMuc')
                  dataEachSim.resMuc = [dataEachSim.resMuc d.resMuc(1:d.kStep0)];                
                  dataEachSim.GRmuc = [dataEachSim.GRmuc [d.resMuc(1:d.kStep0).GRmax]];
              elseif isfield(d, 'GRmuc')
                  dataEachSim.GRmuc = [dataEachSim.GRmuc, columnVector(d.GRmuc(1:d.kStep0))'];
              end
              
              if isfield(d, 'resLum')
                  s = d.resLum(:, 1:size(d.timeLum,1), 1:d.kStep0);
                  s = reshape(s, size(s, 1),size(s, 2) * size(s, 3));
                  dataEachSim.resLum = [dataEachSim.resLum s];
                  M = zeros(size(s));
                  for j = 1:size(s, 2)
                      for jSp = 1:size(s, 1)
                          if isequal(s(jSp, j).GRmax, 0)
                              s(jSp, j).GRmax = zeros(size(s, 1), 1);
                          end
                      end
                      a = [s(:,j).GRmax];
                      if ~isempty(a)
                          M(:,j) = sum(a,2);
                      end
                  end
                  dataEachSim.GRlum = [dataEachSim.GRlum M];
              elseif isfield(d, 'GRlum')
                  s = d.GRlum(:, 1:size(d.timeLum,1), 1:d.kStep0);
                  s = reshape(s, size(s, 1), size(s, 2) * size(s, 3));
                  dataEachSim.GRlum = [dataEachSim.GRlum, s];
              end
              if isfield(d, 'Ct')
                  s = size(d.Ct);
                  dataEachSim.C = cat(2,dataEachSim.C,reshape(d.Ct(:,:,1:d.kStep0),s(1),s(2)*d.kStep0));
              end
              if isfield(d, 'fluxMuc')
                  dataEachSim.fluxMuc = cat(2,dataEachSim.fluxMuc, d.fluxMuc(:,1:d.kStep0));
              end
              if exist(sprintf(['%s_sect%dsave%0' num2str(d0.digit) 'd.mat'], saveNameI,k,j+1),'file')
                  j = j + 1;
              elseif exist(sprintf(['%s_sect%dsave%0' num2str(d0.digit) 'd.mat'], saveNameI,k+1,1),'file')
                  j = 1;
                  k = k + 1;
              else
                  break
              end
              d = load(sprintf(['%s_sect%dsave%0' num2str(d0.digit) 'd.mat'], saveNameI,k,j));
          end
      end
      sumX = sum(dataEachSim.Xlum,1);
      dataEachSim.Alum = zeros(size(dataEachSim.Xlum));
      dataEachSim.Alum(:,sumX>0) = dataEachSim.Xlum(:,sumX>0) * diag(1./sumX(sumX>0));
      sumX = sum(dataEachSim.Xmuc,1);
      dataEachSim.Amuc = zeros(size(dataEachSim.Xmuc));
      dataEachSim.Amuc(:,sumX>0) = dataEachSim.Xmuc(:,sumX>0) * diag(1./sumX(sumX>0));
      data(i) = dataEachSim;
      if isfield(d, 'finish')
          finish(i) = d.finish;
      else
          finish(i) = true;
      end
  else
      fprintf('data for %s not found.\n', saveNameI);
  end
end

% Get the means and standard deviations
dataMean = data(1);
dataSD = data(1);
field = setdiff(fieldnames(dataMean),{'resMuc','resLum'});
for k = 1:numel(field)
    if size(dataMean.(field{k}),2) == 1
        dataK = [data(finish).(field{k})];
        dataMean.(field{k}) = mean(dataK,2);
        dataSD.(field{k}) = std(dataK,0,2);
    elseif size(dataMean.(field{k}),1) == 1
        dataK = [data(finish).(field{k})];
        dataK = reshape(dataK,size(dataMean.(field{k}),2),sum(finish));
        dataMean.(field{k}) = mean(dataK,2);
        dataSD.(field{k}) = std(dataK,0,2);
    else
        dataK = [];
        for j = find(finish(:))'
            dataK = cat(3,dataK,data(j).(field{k}));
        end
        dataMean.(field{k}) = mean(dataK,3);
        dataSD.(field{k}) = std(dataK,0,3);
        switch field{k}
            case {'Amuc','Alum'}
                nSp = size(dataK, 1);
                anaerobeIDs = setdiff(1:nSp, aerobeIDs);
                dataMean.([field{k} 'Anaero']) = mean(sum(dataK(anaerobeIDs,:,:),1),3);
                dataSD.([field{k} 'Anaero']) = std(sum(dataK(anaerobeIDs,:,:),1),0,3);
                dataMean.([field{k} 'Aero']) = mean(sum(dataK(aerobeIDs,:,:),1),3);
                dataSD.([field{k} 'Aero']) = std(sum(dataK(aerobeIDs,:,:),1),0,3);
            case 'Xlum'
                dataMean.XlumSum = mean(sum(dataK,1),3);
                dataSD.XlumSum = std(sum(dataK,1),0,3);
        end
    end
end
end

function figHandle = plotSimVsExp(data, dataSD, expData, expDataSD, T)
biomassPlotRange = [0.2, 1e4];
biomassPlotTick = 10 .^ (0:4);
if nargin < 2
    dataSD = data;
    dataField = fieldnames(dataSD);
    for j = 1:numel(dataField)
        if isnumeric(dataSD.(dataField{j}))
            dataSD.(dataField{j})(:) = 0;
        end
    end
end
row = 2;
col = 3;
figure('Position',  [1 249 831.5000 333]);
linCol = [0    0.4470    0.7410; 0.8500 0.3250 0.0980];

meanT = (T(1:end-1)+T(2:end))/2;

% plot mucosal microbiota
ax(1,1) = subplot(row,col,1);
ax(1,1).NextPlot = 'add';
lin(1,1,1) = plot(data.time,data.AmucAero(1,:),'LineWidth',2,'Color',linCol(1,:));
lin(1,1,2) = plot(data.time,data.AmucAnaero(1,:),'LineWidth',2,'Color',linCol(2,:));
d = [data.AmucAero(1,:)+dataSD.AmucAero(1,:), fliplr(data.AmucAero(1,:)-dataSD.AmucAero(1,:))];
d(d<0) = 0;
d(d>1) = 1;
pat(1,1,1) = fill([data.time; flipud(data.time)], d, lin(1,1,1).Color);
d = [data.AmucAnaero(1,:)+dataSD.AmucAnaero(1,:), fliplr(data.AmucAnaero(1,:)-dataSD.AmucAnaero(1,:))];
d(d<0) = 0;
d(d>1) = 1;
pat(1,1,2) = fill([data.time; flipud(data.time)], d, lin(1,1,2).Color);
ax(1,1).YLim = [0 1];
ax(1,1).YTick = [0 0.25 0.5 0.75 1];

lg = legend({sprintf('facultative\nanaerobes'), sprintf('strict\nanaerobes')});

% plot luminal microbiota
ax(1,2) = subplot(row,col,2);
ax(1,2).NextPlot = 'add';
lin(1,2,1) = plot(data.timeLum(2:end),data.AlumAero(1,2:end),'LineWidth',2,'Color',linCol(1,:));
lin(1,2,2) = plot(data.timeLum(2:end),data.AlumAnaero(1,2:end),'LineWidth',2,'Color',linCol(2,:));
d = [data.AlumAero(1,2:end)+dataSD.AlumAero(1,2:end), fliplr(data.AlumAero(1,2:end)-dataSD.AlumAero(1,2:end))];
d(d<0) = 0;
d(d>1) = 1;
pat(1,2,1) = fill([data.timeLum(2:end); flipud(data.timeLum(2:end))], d, lin(1,2,1).Color);
d = [data.AlumAnaero(1,2:end)+dataSD.AlumAnaero(1,2:end), fliplr(data.AlumAnaero(1,2:end)-dataSD.AlumAnaero(1,2:end))];
d(d<0) = 0;
d(d>1) = 1;
pat(1,2,2) = fill([data.timeLum(2:end); flipud(data.timeLum(2:end))], d, lin(1,2,2).Color);
ax(1,2).YLim = [0 1];
ax(1,2).YTick = [0 0.25 0.5 0.75 1];

% plot luminal microbial mass
ax(1,3) = subplot(row,col,3);
ax(1,3).NextPlot = 'add';
lin(1,3,2) = plot(data.timeLum(2:end),data.XlumSum(1,2:end),'LineWidth',2,'Color', [0 0 0]);%linCol(2,:));
d = [data.XlumSum(1,2:end)+dataSD.XlumSum(1,2:end), fliplr(data.XlumSum(1,2:end)-dataSD.XlumSum(1,2:end))];
d(d<0) = 0;
pat(1,3,2) = fill([data.timeLum(2:end); flipud(data.timeLum(2:end))], d, lin(1,3,2).Color);
ax(1,3).YScale = 'log';
ax(1,3).YLim = biomassPlotRange; %[1e-6,1e-1];
ax(1,3).YTick = biomassPlotTick; %10.^(-6:2:-1);
ax(1,3).YTickLabel(biomassPlotTick == 1) = {'1'};
tlText = {'Mucosal microbiota', ...
    'Luminal microbiota',...
    'Luminal microbial biomass'};%'Biomass in the lumen',''};
ylab = {'relative abundance','relative abundance','relative increase'};
for r = 1:row
    for c = 1:col
        if r == 1
            if isprop(ax(r,c),'XLim')
                ax(r,c).XLim = [0 T(end)];
                tl(r,c) = title(ax(r,c),tlText{r,c});
                ax(r,c).Color = 'none';
                ax(r,c).XTick = T(:)';
            end
            yl(r,c) = ylabel(ax(r,c),ylab{c});
            xl(r,c) = xlabel(ax(r,c),'time (h)');
            for k = 1:2
                if isprop(pat(r,c,k),'LineStyle')
                    pat(r,c,k).LineStyle = 'none';
                    pat(r,c,k).FaceAlpha = 0.4;
                end
            end
        end
    end
end
pos = {'SI.P','SI.M','SI.D','Cec','LI.P','LI.M','LI.D'};
for k = 1:3
    ax(2,k) = subplot(row,col,3+k);
    hold on
    for j = 1:numel(T)
        plot([T(j) T(j)],[0 1],'-k','LineWidth',1.5);
    end
    ax(2,k).XLim = [0 T(end)];
    ax(2,k).Color = 'none';
    ax(2,k).XColor='none';
    ax(2,k).YColor='none';
    ax(2,k).Position(2)=0.43;
    ax(2,k).Position(4)=0.05;
    for j =1:numel(meanT)
        tx(j,k) = text(meanT(j),-0.5,pos{j},'HorizontalAlignment','center');
        tx(j,k).Rotation = 45;
    end
end
lg.Position = [0.85 0.65 0.1422 0.0585];
lg.FontSize = 13;
lg.Box = 'off';

if nargin > 2 && ~isempty(expData)
    if nargin < 4 || isempty(expDataSD)
        expDataSD = expData;
        fn = fieldnames(expData);
        for j = 1:numel(fn)
            expDataSD.(fn{j})(:) = 0;
        end
    end
    %% Add experimental data points
    dataPoint(1, 1) = errorbar(ax(1, 1), meanT, expData.AmucAero, expDataSD.AmucAero, 'o');
    dataPoint(1, 1).Color = linCol(1,:);
    dataPoint(1, 1).MarkerSize = 6.5;
    
    dataPoint(1, 2) = errorbar(ax(1, 1), meanT, expData.AmucAnaero, expDataSD.AmucAnaero, 'o');
    dataPoint(1, 2).Color = linCol(2,:);
    dataPoint(1, 2).MarkerSize = 6.5;
    
    dataPoint(2, 1) = errorbar(ax(1, 2), meanT, expData.AlumAero, expDataSD.AlumAero, 'o');
    dataPoint(2, 1).Color = linCol(1,:);
    dataPoint(2, 1).MarkerSize = 6.5;
    
    dataPoint(2, 2) = errorbar(ax(1, 2), meanT, expData.AlumAnaero, expDataSD.AlumAnaero, 'o');
    dataPoint(2, 2).Color = linCol(2,:);
    dataPoint(2, 2).MarkerSize = 6.5;
    
    dataPoint(3, 1) = plot(ax(1, 3), meanT, expData.XlumRel, 'o');
    dataPoint(3, 1).Color = [0 0 0];
end

ax(1, 1).Position(1) = 0.07;
ax(2, 1).Position(1) = 0.07;
ax(1, 2).Position(1) = 0.36;
ax(2, 2).Position(1) = 0.36;
lg.Position = [0.5341 0.6130 0.1600 0.2748];
lg.FontSize = 11;
lg.String = lg.String(1:2);
ax(1, 3).Position(1) = 0.75;
ax(2, 3).Position(1) = 0.75;

axBG = axes('Position', [0 0 1 1], 'Visible', 'off');
mark(1) = text(axBG, 0.01, 0.97, 'A', 'FontSize', 14, 'FontWeight', 'bold');
mark(2) = text(axBG, 0.3, 0.97, 'B', 'FontSize', 14, 'FontWeight', 'bold');
mark(3) = text(axBG, 0.69, 0.97, 'C', 'FontSize', 14, 'FontWeight', 'bold');

figHandle = struct();
figHandle.fig = gcf;
figHandle.axes = ax;
figHandle.lines = lin;
figHandle.patches = pat;
figHandle.sections = tx;
figHandle.legend = lg;
figHandle.xlabel = xl;
figHandle.ylabel = yl;
figHandle.axBG = axBG;
figHandle.marks = mark;
if exist('dataPoint', 'var')
    figHandle.dataPoints = dataPoint;
end
if 0
    %%
    f = gcf;
    for j = 1:3
        for k = 1:2
            ax(k,j).Color = 'none';
        end
    end
    f.Color = 'none';
    export_fig -png fig_simSpatialGutExample -transparent -r600
end


end
