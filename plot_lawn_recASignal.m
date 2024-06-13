%% notes on script
% this script is modular. Sections contain repeat lines of code so that
% sections may be skipped depending on the needs of your analysis (i.e.
% plotting all replicates, plotting the average, plotting smoothed
% replicate lines). All sections are useful for sanity checking your
% analysis, but some can be removed or skipped once you want the final
% visuals
%% user defined variables
nTime = 7; %number of time points
times = {'0h','4h','8h','12h','24h','36h','48h'}; % time points for plots
load backgroundSignal.mat; % get background YFP and CFP to subtract from measurements
%% START HERE IF ANALYZING NEW CZI FILES, IF USING A SAVED .MAT FILE, PROCEED TO LINE 43
% read in files
files = dir('*.czi'); % get all microscopy files

for i = 1:length(files)
    fileName{i} = files(i).name;
end

%% segment microscopy images
data = {}; 
pks = {}; empty = {}; EcN = {}; clbN = {};

% run segmenting function and save outputs by pks+/- strain (pulled from
% file name)
for i = 1:length(files)
    if contains(fileName{i},'pks') % BAC-pks
        pks = [pks, segment_lawn_recASignal(fileName{i})];
    elseif contains(fileName{i},'empty') % BAC-empty
        empty = [empty, segment_lawn_recASignal(fileName{i})];
    elseif contains(fileName{i},'clbN') % EcN clbN
        clbN = [clbN, segment_lawn_recASignal(fileName{i})];
    else 
        EcN = [EcN, segment_lawn_recASignal(fileName{i})]; % EcN (WT)
    end
end

%%
data = {pks; empty; EcN; clbN}; % save all outputs in one variable

%% save conds
save data.mat data; %save the data variable 

%% index by time
% set indexes for replicates at each time point
timeInx = [1:15; 16:30; 31:45; 46:60; 61:75; 76:90; 91:105];

%% plot each replicate with and without smoothing (for earlier time points/ variation per line)
% visualize smoothing to find optimal window (especially for earlier time
% points)
conditions = {'BAC-pks','BAC-empty','EcN','EcN-clbN'};
decayMat = {}; mat = []; qw = {}; smoothTraceMat = {};
colorArray = jet(nTime);

figure;
for iCond = 1:length(conditions)
    curCond = data{iCond}; % pull out each pks+/pks- strain one at a time
    decay = [];
    for iTime = 1:nTime
        curTime = curCond(timeInx(iTime,:)); % pull out replicates for each time point
        bgndYFP = backgroundSignal{1}(iTime); % YFP background at that time point
        bgndCFP = backgroundSignal{2}(iTime); % CFP background at that time point
        traces = {};
        for i=1:length(curTime)
            scaleFactor = double(curTime{i}.micron_per_pixel); % microns per pixel from czi file metadata
            % calculate the distance jumps between intensity points
            if ~isempty(curTime{i}.x) % in case a mistake is made drawing the lines
                dx = curTime{i}.x(2)-curTime{i}.x(1); % distance between x points on line
                dy = curTime{i}.y(2)-curTime{i}.y(1); % distance between y points on line
                dist_unit = sqrt(dx^2+dy^2); % linear distance between points on line
                % find index of max YFP signal
                junk = flip(curTime{i}.c3); % flip to be compatible with finding um posiitons
                junkCFP = flip(curTime{i}.c4);
                inx = find(junk == max(junk)); % find the index of max YFP
                decay{i}.yfp = junk(1:inx); % signal from end of line (currently position 1) to the max inx
                decay{i}.pos = [inx:-1:1] * dist_unit * scaleFactor; % create x variable 
                % counting from index of max signal down to 1 scaled by distance 
                % between two points and the microns per pixel
                decay{i}.subYFP = decay{i}.yfp-bgndYFP; %get YFP signal minus background from max YFP signal
                decay{i}.subCFP = junkCFP(1:inx)-bgndCFP; %get CFP signal minus background from max YFP signal
            end
        end
        
        % plot raw traces
        subplot(2,4,iCond); hold on;
        traces = {}; x = {};
        for i = 1:length(curTime)
            x{i} = flip(decay{i}.pos); % flip so positions start at peak YFP
            traces{i} = flip(decay{i}.subYFP); % flip YFP minus background so line starts at peak YFP
            plot(x{i},traces{i},'color',colorArray(iTime,:));
        end
        title([conditions(iCond) ' raw recA decay from peak signal']);
        grid on; box on;
        set(gca,'XLim',[0 900]);
        xlabel('distance (microns)'); ylabel('YFP (background subtracted)');

        % plot smoothed traces
        subplot(2,4,iCond+4); hold on;
        smoothTrace = {};
        for i = 1:length(curTime)
            smoothTrace{i} = smoothdata(traces{i},'movmean',20); % smooth each replicate with a moving window
            plot(x{i},smoothTrace{i},'color',colorArray(iTime,:));
        end

        decayMat{iCond,iTime} = traces; % save all recA signal decay lines (same as sectino above)
        smoothTraceMat{iCond,iTime} = smoothTrace; % save all smoothed recA decay lines

        title([conditions(iCond) ' smoothed recA decay from peak signal']);
        grid on; box on;
        set(gca,'XLim',[0 900]);
        xlabel('distance (microns)'); ylabel('YFP (background subtracted)');
        qw{iTime} = plot(nan,'color',colorArray(iTime,:));
    end
end

sgtitle('recA decay on M9 agar');
legend([qw{:}],times)

% sanity check for smoothing: plot 4 hour time point from BAC-pks condition
% that has low signal with lots of noise on a small y-axis scale
figure; hold on;
for i = 1:length(smoothTraceMat{1,2})
    plot(smoothTraceMat{1,2}{i}); hold on;
end
grid on; box on;
title('smoothed window 20')

%% (raw YFP) average each time point after signal starts (4 hours)
distCutoff = 250; % cutoff for length of decay line (by values, not um)
allMean = {}; allStd = {}; allMed = {}; allSte = {};

figure;
for iCond = 1:length(conditions) % loop through pks+/pks- strains
    curMat = []; meanMat = []; stdMat = []; nCells = []; medianMat = []; qw = {};
    steMat = [];
    for iTime = 1:nTime
        curPos = decayMat{iCond,iTime}; % pull out each pks+/pks- strain at each time point
        curMat = [];
        for  i =1:length(curPos)
            cur = curPos{i};
            if length(cur) > distCutoff % only save measurement if decay line is long enough
                curMat = [curMat;cur(1:distCutoff)]; % take measurements up to cutoff (so all replicates are the same length in matrix)
            end
        end
        meanMat = [meanMat;mean(curMat,1)]; % average each time point
        stdMat = [stdMat; std(curMat,[],1)]; % standard deviation at each time point
        medianMat = [medianMat; median(curMat,1)]; % median of each time point
        nCells(iTime) = size(curMat,1); % number of technical replicates for each time point meeting distance cutoff
        steMat = [steMat; std(curMat,[],1)./sqrt((nCells(iTime)))]; % standard error of each time point
    end
    % save the variables
    allMean{iCond} = meanMat;
    allStd{iCond} = stdMat;
    allMed{iCond} = medianMat;
    allSte{iCond} = steMat;

    colorArray = jet(size(meanMat,1));
    % plot YFP with background subtracted
    subplot(1,4,iCond); hold on;
    for i = 2:size(meanMat,1)
        x = 1:scaleFactor:(distCutoff*scaleFactor);
        y = meanMat(i,:)'; %mean
        e = stdMat(i,:)';
        
        %shaded error bar plot
        shadedErrorBar(x,y,e,'LineProps',{'color',colorArray(i,:),'LineWidth',3})
        text(400,(0.03+0.01*i),[num2str(nCells(i)) ' colonies'],'Color',colorArray(i,:))
        qw{i} = plot(nan,'color',colorArray(i,:),'LineWidth',2);
    end
    grid on; box on;
    xlim([0 distCutoff*scaleFactor]);
    title([conditions{iCond} ' raw YFP'])
    xlabel('distance from peak recA signal (um)');
    ylabel('raw YFP')

end
sgtitle('mean recA signal over time')

legend([qw{:}],times(2:end))

%% RAW YFP - get max signal for each decay line
% plot max signal at each time point (to show when signal starts)
avgMax = []; stdMax = []; nCells = [];
condInx = [1 3]; % only look at pks+ conditions

for iCond = 1:length(condInx)
    curMax = []; % clear variables from past iTime
    for iTime = 1:nTime
        curCond = decayMat{condInx(iCond),iTime}; % pull out all measurements for condition and time point
        curMat = []; 
        for  i =1:length(curCond)
            cur = curCond{i}; % loop through each colony measurement
            if length(cur) > 100
                curMat = [curMat;cur(1:100)]; % only pull out if the measurement is at least 100 units
            end
        end
        for iCell = 1:size(curMat,1)
        % get max YFP signal for each replicate 
            curMax(iTime,iCell) = max(curMat(iCell,:));
            if size(curMat,1) < 15 % 15 = number of technical replicates
                curMax(iTime,(size(curMat,1)+1):15) = nan; % add nan for lines less than 100 measurements
            end
        end
        nCells(iCond,iTime) = size(curMat,1); % number of technical replicates in max measurement
        avgMax(iCond,iTime) = nanmean(curMax(iTime,:),2); % take avg of max
        stdMax(iCond,iTime) = nanstd(curMax(iTime,:),[],2); % take error of max
    end
end

% plot bar plot of max signal over time for each strain condition
figure; hold on;
for iCond = 1:length(condInx)
    x = 1:nTime;
    y = avgMax(iCond,:);
    e = stdMax(iCond,:);
    subplot(2,1,iCond); hold on;
    bar(x,y)
    errorbar(x,y,e,'k.')
    grid on; box on;
    set(gca,'xtick',x, 'xticklabel', times)
    title(conditions(condInx(iCond)));
    xlabel('time (hours)');
    ylabel('max YFP signal')
end

%% SMOOTH YFP - get half max from MEAN traces
% smooth signal so that the baseline signal can be more optimally
% determined

%BAC
BACHalfMax = [];
BAC = allMean{1}; % pull out BAC-pks mean decay curves
smoothBAC = smoothdata(BAC,2,'movmean',20); % smooth with moving window
for iTime = 1:nTime
    curTime = smoothBAC(iTime,:); % pull out one time point
    maxYFP = max(curTime); % get the max at that time point
    baseline = median(curTime(170:end)); % determine the baseline at that time point (from tail of decay trace)
    halfYFP = ((maxYFP - baseline)/2) + baseline; % get the value of 50% of the total YFP response
    halfInx = find(curTime >= halfYFP, 1, 'last'); % find the index of the last value above the 50% response value
    BACHalfMax(iTime) = halfInx*scaleFactor; % scale by the micron per pixel value
end

% EcN
EcNHalfMax = [];
EcN = allMean{3}; % pull out EcN pks+ mean decay curves
smoothEcN = smoothdata(EcN,2,'movmean',20); % smooth with moving window
for iTime = 1:nTime
    curTime = smoothEcN(iTime,:); % pull out one time point
    maxYFP = max(curTime); % get the max at that time point
    baseline = median(curTime(170:end)); % determine the baseline at that time point (from tail of decay trace)
    halfYFP = ((maxYFP - baseline)/2) + baseline; % get the value of 50% of the total YFP response
    halfInx = find(curTime >= halfYFP, 1, 'last'); % find the index of the last value above the 50% response value
    EcNHalfMax(iTime) = halfInx*scaleFactor; % scale by the micron per pixel value
end  

meanHalfMax = [BACHalfMax; EcNHalfMax]; % save values into one variable


%% save important variables
rep1.decay = decayMat; % raw YFP traces
rep1.smoothDecay = smoothTraceMat; % smoothed YFP traces
rep1.avgTraces = allMean; % mean of raw YFP traces
rep1.errorTraces = allStd; % standard deviation of raw YFP traces
rep1.max = avgMax; % max YFP
rep1.meanHalfMax = meanHalfMax; % 50% max response of smoothed curves
rep1.medianTrace = allMed; % median of raw YFP traces
rep1.STE = allSte; % standard error of raw YFP traces

save rep1.mat rep1