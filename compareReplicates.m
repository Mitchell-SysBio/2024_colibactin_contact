%% load data
load rep1.mat
load rep2.mat
load rep3.mat

reps = {rep1; rep2; rep3}; % save as one variable

%% user definitions
times = {'0h','4h','8h','12h','24h','36h','48h'}; % time points for plots
scaleFactor = 3.632; % um per pixel from czi metadata
nTime = length(times);
%% pull out max signal
BACmax = []; BACnegMax = []; EcNmax = []; EcNNegMax = []; % variable holders

% pull out the mean max for each replicate
for i = 1:3
    cur = reps{i};
    BACmax(i,:) = cur.max(1,:);
    BACnegMax(i,:) = cur.max(2,:);
    EcNmax(i,:) = cur.max(3,:);
    EcNNegMax(i,:) = cur.max(4,:);
end

% average maxes and find error
avgBACMax = mean(BACmax,1);
stdBACMax = std(BACmax,[],1);
avgBnegMax = mean(BACnegMax,1);
stdBnegMax = std(BACnegMax,[],1);
avgEcNMax = mean(EcNmax,1);
stdEcNMax = std(EcNmax,[],1);
avgEnegMax = mean(EcNNegMax,1);
stdEnegMax = std(EcNNegMax,[],1);

% put both pks+ and pks- in one variable
BAC = [avgBACMax; avgBnegMax];
EcN = [avgEcNMax; avgEnegMax];
errBAC = [stdBACMax; stdBnegMax];
errEcN = [stdEcNMax; stdEnegMax];

% t-test for pks+ versus pks- at each time point for both strains
for i = 1:nTime
    [hBAC(i) pBAC(i)] = ttest2(BACmax(:,i),BACnegMax(:,i));
    [hEcN(i) pEcN(i)] = ttest2(EcNmax(:,i), EcNNegMax(:,i));
end

figure; hold on;
% plot average max YFP for BAC
subplot(2,1,1); hold on;
b = bar(BAC');
y = BAC';
err = errBAC';

[ngroups,nbars] = size(y); %holders for number of bars
x = nan(nbars, ngroups); %holders for number of bars
for i = 1:nbars
   x(i,:) = b(i).XEndPoints; %get middle of each bar on x axis
end
hold on
errorbar(x', y, err, 'k.'); %plot error bars
plot(x(1,:),BACmax,'.r','MarkerSize',10); % plot technical replicate values
plot(x(2,:),BACnegMax,'.r','MarkerSize',10); % plot technical replicate values
hold off
title('BAC max signal')
set(gca,'xtick',[1:7],'xticklabel',times);
ylabel('YFP intensity')
grid on; box on;
ylim([0 0.18]);
for i = 1:length(pBAC)
    text([i],0.15,num2str(pBAC(i)));
end

% plot average max YFP for EcN
subplot(2,1,2); hold on;
b = bar(EcN');
y = EcN';
err = errEcN';

[ngroups,nbars] = size(y); %holders for number of bars
x = nan(nbars, ngroups); %holders for number of bars
for i = 1:nbars
   x(i,:) = b(i).XEndPoints; %get middle of each bar on x axis
end
hold on
errorbar(x', y, err, 'k.'); %plot error bars
plot(x(1,:),EcNmax,'.r','MarkerSize',10); % plot technical replicate values
plot(x(2,:),EcNNegMax,'.r','MarkerSize',10); % plot technical replicate values
hold off
title('EcN max signal')
set(gca,'xtick',[1:7],'xticklabel',times);
ylabel('YFP intensity')
grid on; box on;
ylim([0 0.15]);
for i = 1:length(pEcN)
    text([i],0.13,num2str(pEcN(i)));
end


%% plot half max averaged over each replicate (taken from smoothing the median curves)
BAC = []; EcN = []; % variable holders

% pull out Half max from mean of technical replicates (not median)
for i = 1:3
    cur = reps{i}.meanHalfMax;
    BAC(i,:) = cur(1,:);
    EcN(i,:) = cur(2,:);
end

figure; hold on;
% plot average distance to 50% max response for BAC pks+
subplot(2,1,1); hold on;
y = mean(BAC,1);
e = std(BAC,[],1);
bar(y(2:end))
errorbar(y(2:end),e(2:end),'k.');
plot([1:6],BAC(:,2:end),'.r','MarkerSize',10)
set(gca,'xtick',[1:6],'xticklabel',times(2:end));
ylabel('distance to half max signal (um)')
title('BAC-pks half max')
grid on; box on;

% plot average distance to 50% max response for EcN pks+
subplot(2,1,2); hold on;
y = mean(EcN,1);
e = std(EcN,[],1);
bar(y(2:end))
errorbar(y(2:end),e(2:end),'k.');
plot([1:6],EcN(:,2:end),'.r','MarkerSize',10)
set(gca,'xtick',[1:6],'xticklabel',times(2:end));
ylabel('distance to half max signal (um)')
title('EcN half max')
grid on; box on;

allHalfMax = {BAC, EcN}; % save the averages in a new variable to plot below

%% pull out averaged decays by time to plot replicates on same graph
T12 = []; T24 = []; T36 = []; T48 = []; T4 = []; T8 = []; % holders
curTimes = {'4h','8h','12h','24h','36h','48h'}; % times being plotted
conditions = {'BAC-pks','EcN'}; % labels
ylims = [0.2 0.15]; % upper y limit for BAC and EcN respectively

for iStrain = [1 3] % indexes of BAC and EcN pks+
    figure; hold on;
    for i = 1:length(curTimes) % loop through each time point
        iTime = i+1; % only taking from 4 hours and on
        x = 1:scaleFactor:(length(rep1.medianTrace{iStrain}(iTime,:))*scaleFactor); % get x variable in microns
        y = [rep1.avgTraces{iStrain}(iTime,:); rep2.avgTraces{iStrain}(iTime,:); rep3.avgTraces{iStrain}(iTime,:)]; % pull out the mean trace for each biological replicate
        e = [rep1.errorTraces{iStrain}(iTime,:); rep2.errorTraces{iStrain}(iTime,:); rep3.errorTraces{iStrain}(iTime,:)]; % pull out the standard deviation of each trace for each biological replicate
        
        % plot mean plus error bars for each biological replicate with each
        % time in a different subplot
        subplot(3,2,i); hold on;
        shadedErrorBar(x,y(1,:),e(1,:),'LineProps',{'LineWidth',2,'color','r'})
        shadedErrorBar(x,y(2,:),e(2,:),'LineProps',{'LineWidth',2,'color','b'})
        shadedErrorBar(x,y(3,:),e(3,:),'LineProps',{'LineWidth',2,'color','g'})
        grid on; box on;
        title(curTimes(i))
        xlabel('distance from peak YFP (um)');
        ylabel('YFP signal')
        if iStrain == 1
            inx = 1; % index for labeling pks+ strain
        else
            inx = 2; % index for labeling pks+ strain
        end
        % mark vertical lines at distance of 50% max signal for each
        % biological replicate
        xline(allHalfMax{inx}(1,iTime),'r')
        xline(allHalfMax{inx}(2,iTime),'b')
        xline(allHalfMax{inx}(3,iTime),'g')
        ylim([0 ylims(inx)]);
    end
    sgtitle(conditions(inx))
end
