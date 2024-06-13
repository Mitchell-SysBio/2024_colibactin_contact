%% user defined variables
nTime = 1; %number of time points
times = {'24h'};
noTouchCutOff = 10; % distance in microns to say colonies are not touching 
touchCutoff = 4; %max distance in microns between colnoies to say they are touching

%% SKIP TO LINE 37 IF YOU HAVE ALREADY ANALYZED FILES AND ARE JUST PLOTTING read in files
files = dir('*.czi');

for i = 1:length(files)
    fileName{i} = files(i).name;
end

%% segment microscopy images and separate by pks+/pks- and strain
BACPos = {}; BACNeg = {}; EcNPos = {}; EcNNeg = {};
for i = 1:length(files)
    if contains(fileName{i},'empty')
        BACNeg = [BACNeg, segment_recASignal(fileName{i})];
    elseif contains(fileName{i},'BAC')
        BACPos = [BACPos, segment_recASignal(fileName{i})];
    elseif contains(fileName{i},'clbN')
        EcNNeg = [EcNNeg, segment_recASignal(fileName{i})];
    elseif contains(fileName{i},'EcN')
        EcNPos = [EcNPos, segment_recASignal(fileName{i})];
    end
end

BAC = {BACPos; BACNeg};
EcN = {EcNPos; EcNNeg};

%% save conds
conds = [BAC, EcN];
save BAC.mat BAC; 
save EcN.mat EcN;
save conds.mat conds; %save the entire conds variable

%% START HERE IF ALREADY ANALYZED IMAGES AND ARE PLOTTING set current strain to plot
cond = BAC; % set the strain you want to plot

%% USING MARKED LINES extract trace from peak YFP signal of touching and non touching for both pos and neg
% get all traces and means with error in subplots for non touching pks+,
% touching pks+, touching pks-
decayMat = []; mat = []; nonTouch = [];
allPosPos = {}; allYFPPos = {}; allRepPos = {};
allPosNeg = {}; allYFPNeg = {}; allRepNeg = {};
allDecay = []; allDist = [];

for iCond = 1:2
    col = cond{iCond};
    decay = [];
    for i=1:length(col)
        scaleFactor = double(col{i}.micron_per_pixel);
        % calculate the distance jumps between intensity points
        if ~isempty(col{i}.x) & ~isempty(col{i}.edges)
            dx = col{i}.x(2)-col{i}.x(1); dy = col{i}.y(2)-col{i}.y(1);
            dist_unit = sqrt(dx^2+dy^2); 
            % calculate location of colony edges along line
            if size(col{i}.edges,1) == 2
                dxRep = col{i}.x(1) - col{i}.edges(1,1);
                dyRep = col{i}.y(1) - col{i}.edges(1,2);
                dxTox = col{i}.x(1) - col{i}.edges(2,1);
                dyTox = col{i}.y(1) - col{i}.edges(2,2);
                distRep = sqrt(dxRep^2 + dyRep^2);
                distTox = sqrt(dxTox^2 + dyTox^2);
            end
            inx = round(distTox,0);
            decay{i}.pos = [inx:-1:1]*dist_unit*scaleFactor; %get positions in um from edge of toxic colony to end of roi line
            decay{i}.rep = decay{i}.pos(1) - distRep*scaleFactor; %get distance between colonies
            decay{i}.yfp = col{i}.c3(1:inx); % get yfp signal from edge of producer colony
            decay{i}.edge = inx; %save index of edge of toxic colony
            decay{i}.repEdge = round(distRep,0); %save index of edge of reporter colony
            decay{i}.dist = (decay{i}.edge-decay{i}.repEdge) * scaleFactor; %get distance between colonies
        else
            decay{i}.pos = nan;
            decay{i}.rep = nan;
            decay{i}.yfp = nan;
            decay{i}.edge = nan;
            decay{i}.dist = nan;
        end
    end

    pos = {}; yfp = {}; rep = {};
    noSignalPos = {}; noSignalYFP = {}; noSignalRep = {};
    noTouchTraceMax = {}; noTouchTrace = {}; touchTrace = {};
    noTouchDist = []; touchDist = []; 
    if iCond == 1 % pks+
        for i=1:length(col)
            % non-touching colonies
            if decay{i}.rep >= noTouchCutOff % use marked edge to get touching v. non-touching by eye 
                pos = [pos;decay{i}.pos];
                yfp = [yfp;decay{i}.yfp];
                rep = [rep; decay{i}.rep];
                curColor = [1 0 0];
                maxYFP = find(decay{i}.yfp == max(decay{i}.yfp)); % find peak signal
                noTouchTraceMax{end+1} = flip(decay{i}.yfp(1:maxYFP)); %decay from peak signal
                junk = decay{i}.yfp(1:decay{i}.edge); %decay from edge of tox colony
                junk(maxYFP:decay{i}.edge) = nan; %replace values between edge of tox and peak YFP with nan
                noTouchTrace{end+1} = flip(junk); % flip for plotting
                noTouchDist(end+1) = decay{i}.rep; % save distance between colonies
            % touching colonies
            elseif decay{i}.rep <= touchCutoff
                pos = [pos;decay{i}.pos];
                yfp = [yfp;decay{i}.yfp];
                rep = [rep; decay{i}.rep];
                curColor = [0 0 1];
                maxYFP = find(decay{i}.yfp == max(decay{i}.yfp)); % find peak YFP signal
                if ~isnan(decay{i}.edge)
                    junk = decay{i}.yfp(1:decay{i}.edge); %decay from edge of tox colony
                    junk(maxYFP:decay{i}.edge) = nan; %replace values from the edge to the max YFP with nan
                    touchTrace{end+1} = flip(junk); % flip for plotting
                    touchDist(end+1) = decay{i}.rep; % save distance between colonies (should be less than cutoff set at beginning of script)
                end
            end
            allPosPos = [allPosPos; decay{i}.pos]; % save into one variable
            allYFPPos = [allYFPPos; decay{i}.yfp];
            allRepPos = [allRepPos; decay{i}.rep];
        end
    else % pks- (don't use peak YFP signal, since it's basically flat)
        for i=1:length(col)
            %non-touching colonies only
            if decay{i}.rep >= noTouchCutOff
                pos = [pos;decay{i}.pos];
                yfp = [yfp;decay{i}.yfp];
                rep = [rep; decay{i}.rep];
                curColor = [1 0 0];
                noTouchTrace{end+1} = flip(decay{i}.yfp(1:decay{i}.edge)); % flip for plotting
                noTouchDist(end+1) = decay{i}.rep; % save distance between colonies
            % touching colonies
            elseif decay{i}.rep <= touchCutoff 
                pos = [pos;decay{i}.pos];
                yfp = [yfp;decay{i}.yfp];
                rep = [rep; decay{i}.rep];
                curColor = [0 0 1];
                if ~isnan(decay{i}.edge)
                    junk = decay{i}.yfp(1:decay{i}.edge); %decay from edge of tox colony
                    touchTrace{end+1} = flip(junk); % flip for plotting
                    touchDist(end+1) = decay{i}.rep; % save distance between colonies
                end
            end
            allPosNeg = [allPosNeg; decay{i}.pos];
            allYFPNeg = [allYFPNeg; decay{i}.yfp];
            allRepNeg = [allRepNeg; decay{i}.rep];
        end
    end


    decayMat{iCond} = {noTouchTrace; touchTrace};
    distMat{iCond} = {noTouchDist; touchDist};
end

%% plot
titles = {'non-touching pks+','touching pks+','non-touching pks-','touching pks-'};
distCutoff = 150; %minimum length of roi lines
ymax = 0.13; %upper y lim by strain

halfMax = []; meanTrace = {};
figure; hold on;
for iCond = 1:2
    for iTouch = 1:2 %1 is no touch, 2 is touch
        if iCond ==1 & iTouch == 1 %no touch pks+
            inx = 1;
        elseif iCond == 1 & iTouch == 2 %touch pks+
            inx = 2;
        elseif iCond ==2 & iTouch ==1 %no touch pks-
            inx = 3;
        else
            inx = 4; %touch pks-
        end
        curMat = []; normMat = [];
        subplot(2,4,inx)
        for i = 1:length(decayMat{iCond}{iTouch})
            cur = decayMat{iCond}{iTouch}{i};
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
        title(titles{inx})
        xlim([1 distCutoff*scaleFactor]);
        ylim([0 ymax]);
        text(100,ymax-0.02,['n = ' num2str(size(curMat,1))])

        subplot(2,4,(inx+4)) %plot shaded error bars
        shadedErrorBar(x,meanMat,stdMat)
        grid on; box on;
        title(titles{inx})
        xlim([1 distCutoff*scaleFactor]);
        ylim([-0.02 ymax]);
        xline(halfInx*scaleFactor); % vertical line at dist50
        halfMax(:,inx) = [halfInx*scaleFactor, halfDecay+nanmean(meanMat((distCutoff-20):distCutoff))]; 
        maxYFP(inx) = max(meanMat);
        meanTrace{iCond}(iTouch,:) = meanMat;
        stdTrace{iCond}(iTouch,:) = stdMat;
    end
end

%% save variables
rep4.halfMax = halfMax;
rep4.max = maxYFP;
rep4.trace = meanTrace;
rep4.stdTrace = stdTrace;
save rep4.mat rep4
