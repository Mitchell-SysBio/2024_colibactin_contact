%% user variables
curDate = '20240827'; %for saving files at the end
wellInx = 5; %column with well location from metadata
ctrlNum = 0; %number of control samples you have in your metadata/.fcs files
repNum = 3; %number of replicates you have per sample (used for plotting/taking averages)
%%
load xyCoordinate.mat;
xyCells = xyCoord.Cells;
xySingle = xyCoord.Single;
xyPE = xyCoord.PI;

%% max signal value
maxSigVal = 2^18; %assumed maximum fluorescent signal the machine records

%% read in fcs file names
fcsFiles = dir('*.fcs'); %get all .fcs file names
metadata = readcell('20240827_metadata.xlsx'); %get metadata descriptions
sampleNum = length(metadata)-ctrlNum; %get the number of samples

%% read files 
fcsdat = {}; fcshdr = {};
for iFile = 1:length(fcsFiles)
    [fcsdat{iFile}, fcshdr{iFile}] = fca_readfcs(fcsFiles(iFile).name); %read in all .fcs files
end
allNames = replace({fcsFiles(:).name},'_',' '); %get fcs file names for all samples

%% sort files into same order as metadata
sampInx = [];
for iWell = 1:size(metadata,1)
    curWell = metadata(iWell,wellInx);
    %curPlate = metadata(iWell,plateInx);
    sampInx(iWell,:) = find(contains(allNames,curWell)); %get index of each file in the order of the metadata
end
orderedFcsFiles = fcsdat(sampInx); %reorder the fcs files to match the metadata
description = join(cellfun(@string,metadata)); %combine the descriptors to one string for each sample

%% user defined variables
subset = [1, 5, 10,15, 26];
%% Extract vectors for each parameter (indexes can be found in fcshdr.par)
FSCA = {}; FSCH = {}; SSCA = {}; CFP = {}; YFP = {};
for i = 1:length(subset)
    FSCA{i} = orderedFcsFiles{subset(i)}(:,1);
    FSCH{i} = orderedFcsFiles{subset(i)}(:,2);
    SSCA{i} = orderedFcsFiles{subset(i)}(:,4);
    PE{i} = orderedFcsFiles{subset(i)}(:,7);
end

%% Gate on cells
myColor = jet(length(subset));
fh = figure;
hold on
for i = 1:length(subset)
    plot(FSCA{i},SSCA{i},'o','MarkerSize',2,'MarkerEdgeColor',myColor(i,:),'MarkerFaceColor',myColor(i,:)); %scatter plot of FSC-A x SSC-A
end
xlabel('FSC-A'); ylabel('SSC-A');
set(gca,'Yscale','log','Xscale','log');
legend(description(subset),'location','eastoutside');
cells = drawpolygon(gca); %draw polygon to gate on cells
xyCells = cells.Position; %get coordinates of the polygon points
close(fh);
tfCells = {};
for i = 1:length(subset)
    tfCells{i} = inpolygon(FSCA{i},SSCA{i},xyCells(:,1),xyCells(:,2)); %get indexes of cells that are inside that gate
end

%% Gate on single cells
fh = figure;
hold on
for i = 1:length(subset)
    plot(FSCA{i}(tfCells{i}),FSCH{i}(tfCells{i}),'o','MarkerSize',2,'MarkerEdgeColor',myColor(i,:),'MarkerFaceColor',myColor(i,:)); %scatter plot of FSC-A x FSC-H
end
xlabel('FSC-A'); ylabel('FSC-H');
set(gca,'Yscale','log','Xscale','log');
legend(description(subset),'location','eastoutside');
singleCell = drawpolygon(gca); %draw polygon to gate on single cells
xySingle = singleCell.Position; %get coordinates
close(fh);
tfSingle = {};
for i = 1:length(subset)
    tfSingle{i} = inpolygon(FSCA{i},FSCH{i},xySingle(:,1),xySingle(:,2)); %get indexes of cells inside the gate
end

%% Get indexes of single cells within cells gate
tfSingleCells = {};
for i =1:length(subset)
    tfSingleCells{i} = tfSingle{i} & tfCells{i}; %get indexes of cells within both the cell and single cell gates
end

%% Gate on PI positive cells
fh = figure;
hold on
for i = 1:length(subset)
    plot(FSCA{i}(tfSingleCells{i}),PE{i}(tfSingleCells{i}),'o','MarkerSize',2,'MarkerEdgeColor',myColor(i,:),'MarkerFaceColor',myColor(i,:)) %scatter plot of CFP in single cells
end
xlabel('FSC-A'); ylabel('PE-A');
yline(maxSigVal);
legend(description(subset),'location','eastoutside');
xlim([10 100000]);
set(gca,'Yscale','log','Xscale','log');
grid on;
posPE = drawpolygon(gca); %draw rectangle to gate PE cells
xyPE = posPE.Position; %get coordinates of rectangle
close(fh)
tfPE = {};
for i = 1:length(subset)
    tfPE{i} = inpolygon(FSCA{i},PE{i},xyPE(:,1),xyPE(:,2)); %get indexes of cells within PE pos gate
end

%% sanity check: plot sum of cells in each gate (should be smaller with each subsequent gating step)
figure;
for i = 1:length(subset)
    curData = [length(FSCA{i}),sum(tfCells{i}),sum(tfSingle{i}),sum(tfSingleCells{i}),sum(tfPE{i}&tfSingleCells{i})];
    subplot(3,3,i)
    bar(curData)
    grid on
    title(description(subset(i)));
    set(gca,'xticklabel',{'total events','Cells','Single','Single + Cells','PE + Single + Cells'})
    xtickangle(45)
end

%% Apply gating strategy to all fcs files
%extract cell vectors for each parameter for all files
FSCA = {}; FSCH = {}; SSCA = {}; PE = {};
for i = 1:length(fcsFiles)
    FSCA{i} = orderedFcsFiles{i}(:,1);
    FSCH{i} = orderedFcsFiles{i}(:,2);
    SSCA{i} = orderedFcsFiles{i}(:,4);
    PE{i} = orderedFcsFiles{i}(:,7);
end
%get indexes of cells within each gate
for i = 1:length(fcsFiles)
    tfCells{i} = inpolygon(FSCA{i},SSCA{i},xyCells(:,1),xyCells(:,2)); %get indexes of cells that are inside that gate
    tfSingle{i} = inpolygon(FSCA{i},FSCH{i},xySingle(:,1),xySingle(:,2)); %get indexes of cells inside the gate
    tfPE{i} = inpolygon(FSCA{i},PE{i},xyPE(:,1),xyPE(:,2)); %get indexes of cells within PE pos gate
end
% get indexes of cells that are within both the cell and single cell gates
tfSingleCells = {};
for i =1:length(fcsFiles)
    tfSingleCells{i} = tfSingle{i} & tfCells{i}; %get indexes of cells within both the cell and single cell gates
    tfSingleCellsPE{i} = tfSingleCells{i} & tfPE{i}; %get indexes of single cells in the positive PE gate
end


%% concatenate replicates into a single array
condIdx = [1:3;7:9;13:15;19:21;];
ctrl = []; PI = [];

for iCond = 1:4
    ctrl = PE(condIdx(iCond,:)+3);
    PI = PE(condIdx(iCond,:));
    for i = 1:length(condIdx(iCond,:))
        catCTRLStrain{i} = ctrl{i}(tfSingleCells{condIdx(iCond,i)+3});
        catCtrl{iCond} = vertcat(catCTRLStrain{:});
        catPIStrain{i} = PI{i}(tfSingleCells{condIdx(iCond,i)});
        catPI{iCond} = vertcat(catPIStrain{:});
    end
end


%% histogram of strains/conditions
figure; hold on;
strains = {'BAC-pks','BAC-empty','EcN','EcN clbN'};
posCtrl = PE{end}(tfSingleCells{end});
inx = [1 3];
for i = 1:length(inx)
    subplot(length(inx),1,i); hold on;
    histogram(catPI{inx(i)},10.^[0:0.075:6],'Normalization','probability','EdgeColor','none')
    histogram(catPI{inx(i)+1},10.^[0:0.075:6],'Normalization','probability','EdgeColor','none')
    histogram(posCtrl,10.^[0:0.075:6],'Normalization','probability','EdgeColor','none')
    histogram(catCtrl{inx(i)},10.^[0:0.075:6],'Normalization','probability','EdgeColor','none')
    grid on; box on;
    ylim([0 .14])
    xlim([10^0 10^5.2])
    set(gca,'xscale','log')
    title(strains(inx(i)))
    xline(4.5*10^3)
end
legend({'pks+','pks-','heat killed control','unstained control'});

%% compare % PI positive by histogram visual cutoff
rawPE = []; rawSingle = []; 
histCutoff = 10^4.5;
for i = 1:length(fcsFiles)
    rawSingle = [rawSingle, sum(tfSingleCells{i})]; %sum of single cells
    rawPE = [rawPE, sum(PE{i}(tfSingleCells{i})>histCutoff)]; %sum of single cells above visual PI cutoff
end

percPE = (rawPE./rawSingle)*100; %percent of PI positive cells

%% take average of replicates
meanPE = []; stdPE = [];
for i = [1:repNum:sampleNum-2]
    meanPE = [meanPE, mean([percPE(i),percPE(i+1),percPE(i+2)])]; %mean of replicates
    stdPE = [stdPE, std([percPE(i),percPE(i+1),percPE(i+2)])]; %standard deviation of replicates
end
condKey = description([1:repNum:sampleNum]); %get descriptions for averaged samples


%% plot mean % recA in different strains
x = [1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8];
strain = metadata(1:end-2,1);
figure; hold on;
b = bar(1:length(condKey(1:8)),meanPE);
errorbar(1:length(condKey(1:8)),meanPE,stdPE,'k.')
plot(x, percPE(1:24),'.','MarkerSize',10)

grid on; box on;
set(gca,'xtick',1:length(condKey)-1,'xticklabel',condKey);
xlabel('bacteria strain');
ylabel('% PI positive cells')


%% t-test for strains

%BAC
x = percPE(1:3);
y = percPE(7:9);
[h, p] = ttest2(x,y,'Vartype','unequal');

%EcN
x1 = percPE(13:15);
y1 = percPE(19:21);
[h1, p1] = ttest2(x1,y1,'Vartype','unequal');

%% save % PI
fileStr = [curDate,'_','PI_results.xlsx'];
mat2save = [rawPE', percPE', rawSingle'];
headerStr = {'condition','# PI','% PI','total single cells'};
writecell(headerStr,fileStr,'Range','A1');
writematrix(description,fileStr,'Range','A2');
writematrix(mat2save,fileStr,'Range','B2');

% %% save gating coordinates to use in the future
xyCoord.Cells = xyCells;
xyCoord.Single = xySingle;
xyCoord.PI = xyPE;

save('xyCoordinate.mat','xyCoord');


