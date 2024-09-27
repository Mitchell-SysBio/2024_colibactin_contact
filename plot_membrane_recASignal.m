%% user defined variables
nTime = 1; %number of time points
times = {'24h'};
%% read in files
files = dir('*.czi');

for i = 1:length(files)
    fileName{i} = files(i).name;
end

%% segment microscopy images
pos = {}; neg = {};
for i = 1:length(files)
    if contains(fileName{i},'pks')
        pos = [pos, segment_membrane_recASignal(fileName{i})];
    elseif contains(fileName{i},'empty')
        neg = [neg, segment_membrane_recASignal(fileName{i})];
    end
end

conds = {pos; neg};

%% save conds
save conds.mat conds; %save the entire conds variable

%% extract trace from peak YFP signal of touching and non touching for both pos and neg
% get all traces and means with error in subplots for non touching pks+,
% touching pks+, touching pks-
decayMat = []; mat = [];

for iCond = 1:2
    col = conds{iCond};
    decay = [];
    for i=1:length(col)
        scaleFactor = double(col{i}.micron_per_pixel);
        % calculate the distance jumps between intensity points
        if ~isempty(col{i}.x) & ~isempty(col{i}.edges)
            dx = col{i}.x(2)-col{i}.x(1); dy = col{i}.y(2)-col{i}.y(1);
            dist_unit = sqrt(dx^2+dy^2); 
            % calculate location of colony edges along line
            if size(col{i}.edges,2) == 2
                dxRep = col{i}.x(end) - col{i}.edges(1,1);
                dyRep = col{i}.y(end) - col{i}.edges(1,2);
                distRep = sqrt(dxRep^2 + dyRep^2);
            end
            inx = round(distRep,0);
            decay{i}.pos = [inx:-1:1]*dist_unit*scaleFactor; %get positions in um from edge of membrane to end of roi line
            junk = flip(col{i}.c2);
            decay{i}.yfp = junk(1:inx); % get yfp signal from end of ROI to membrane
            decay{i}.edge = inx; %save index of membrane
        else
            decay{i}.pos = nan;
            decay{i}.yfp = nan;
            decay{i}.edge = nan;
        end
    end
    
    pos = {}; yfp = {}; rep = {};
    TraceMax = {}; Trace = {}; maxYFP = [];
    if iCond == 1 % pks+ 
        for i=1:length(col)
            pos = [pos;decay{i}.pos];
            yfp = [yfp;decay{i}.yfp];
            curColor = [0 0 1];
            if ~isnan(decay{i}.edge)
                junk = decay{i}.yfp;
                Trace{end+1} = flip(junk);
                maxYFP(end+1) = max(Trace{i});
            end
        end
    else % pks- (don't use peak YFP signal, since it's basically flat)
        for i=1:length(col)
            if ~(isnan(decay{i}.yfp))
                pos = [pos;decay{i}.pos];
                yfp = [yfp;decay{i}.yfp];
                curColor = [1 0 0];
                Trace{end+1} = flip(decay{i}.yfp);
                maxYFP(end+1) = max(Trace{end});
            end
        end
    end

    decayMat{iCond} = Trace;
    maxes{iCond} = maxYFP;
end


%% plot
titles = {'pks+','pks-'};
distCutoff = 400; %minimum length of roi lines
ymax = 0.25; %upper y lim

halfMax = [];
figure; hold on;
for iCond = 1:2
    curMat = []; normMat = [];
    subplot(2,2,iCond)
    for i = 1:length(decayMat{iCond})
        cur = decayMat{iCond}{i};
        if length(cur) < distCutoff %find lines shorter than the distance we want to plot
            cur(end+1:distCutoff) = nan; %add nans to end of data
        else
            cur = cur(1:distCutoff); %make sure all are the same length
        end
        baseline = mean(cur((distCutoff-20):distCutoff)); % get average yfp intensity along flat part of decay line
        curMat = [curMat; (cur-baseline)]; %put all rois with nans into matrix
    end
    meanMat = nanmean(curMat,1); %get mean (excluding nans)
    stdMat = nanstd(curMat,[],1);
    halfDecay = (max(meanMat)-nanmean(meanMat((distCutoff-20):distCutoff)))/2; %find what half the max signal is
    halfInx = max(find(meanMat >= (halfDecay+nanmean(meanMat((distCutoff-20):distCutoff))))); %index where half max signal is found (distance)
    x = 1:scaleFactor:(distCutoff*scaleFactor); %set x variable in microns
    plot(x,curMat,'b'); hold on; %plot all traces
    plot(x,meanMat,'k','LineWidth',3) %plot mean of traces
    grid on; box on;
    title(titles{iCond})
    xlim([1 distCutoff*scaleFactor]);
    ylim([0 ymax]);
    text(100,ymax-0.02,['n = ' num2str(size(curMat,1))])

    subplot(2,2,iCond+2) %plot shaded error bars
    shadedErrorBar(x,meanMat,stdMat)
    grid on; box on;
    title(titles{iCond})
    xlim([1 distCutoff*scaleFactor]);
    ylim([-0.01 ymax]);
    xline(halfInx*scaleFactor);
    halfMax(iCond) = halfInx*scaleFactor;
    maxYFP(iCond) = max(meanMat);
end

%% plot max YFP as bar graph
means = [mean(maxes{1}), mean(maxes{2})];
stds = [std(maxes{1}), std(maxes{2})];

figure; hold on;
bar(means)
errorbar(means,stds,'k.');
plot(1,maxes{1},'o','MarkerSize',4,'MarkerFaceColor','b','MarkerEdgeColor','k')
plot(2,maxes{2},'o','MarkerSize',4,'MarkerFaceColor','b','MarkerEdgeColor','k')
grid on; box on;
set(gca,'xtick',[1:2], 'xticklabel',{'pks+','pks-'})

[h p] = ttest2(maxes{1},maxes{2},'Vartype','unequal');
text(2,0.2,['p = ' num2str(p)])
text(2.5,0.2,['n = ',num2str(length(maxes{1}))]);
text(2.5,0.15,['n = ',num2str(length(maxes{2}))]);


%% save variables
save membraneHalfMax1.mat halfMax;
save membraneMax1.mat maxYFP
