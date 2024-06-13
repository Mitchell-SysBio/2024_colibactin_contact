% read in files from each replicate
load rep1.mat
load rep2.mat
load rep4.mat

repsHalf = [rep1.halfMax(1,:); rep2.halfMax(1,:); rep4.halfMax(1,:)];
repsMax = [rep1.max; rep2.max; rep4.max];
repsTrace = [rep1.trace; rep2.trace; rep4.trace];
scaleFactor = 3.632; % um per pixel from czi metadata

%% get average dist50 for touching and find distance of that signal in non-touching
halfSignal = [rep1.halfMax(2,:); rep2.halfMax(2,:); rep4.halfMax(2,:)]; % pull out signal at dist50 for contacting colonies
avgTouch = mean(halfSignal(:,2)); % average across replicates
for i = 1:3
    curTrace = repsTrace{i,1}(1,:); % pull out non-contacting average decay for each replicate
    dist50Signal(i) = find(curTrace >= avgTouch,1,'last')*scaleFactor; % find the distance where the 50% signal from contacting is observed
end

%% plot example with correct non-touching xlines to overlay on signal decay curves
figure; hold on;
plot([1:scaleFactor:(150*scaleFactor)],curTrace)
xline(dist50Signal)
grid on; box on;


%% plot half max bar graph
meanHalfMax = [mean(dist50Signal), mean(repsHalf(:,2),1)]; % take means
stdHalfMax = [std(dist50Signal), std(repsHalf(:,2),[],1)]; % calculate standard deviation

% plot
figure; hold on;
bar(meanHalfMax)
errorbar(meanHalfMax,stdHalfMax,'k.');
plot([dist50Signal; repsHalf(:,2)'],'.','MarkerSize',10)
grid on; box on;
set(gca,'xtick',[1:2], 'xticklabel',{'non-contacting','contacting'})

% test and overlay p-value on plot
[h p] = ttest2(dist50Signal,repsHalf(:,2));
text(1,120,['p = ' num2str(p)])

%% plot max signal bar graph
meanMax = mean(repsMax(:,1:2),1); % get means
stdMax = std(repsMax(:,1:2),[],1); % calculate standard deviation

%plot
figure; hold on;
bar(meanMax)
errorbar(meanMax,stdMax,'k.');
plot(repsMax(:,1:2)','.','MarkerSize',10)
grid on; box on;
set(gca,'xtick',[1:2], 'xticklabel',{'non-contacting','contacting'})

% t-test and overlay p-values on plot
[h p] = ttest2(repsMax(:,1),repsMax(:,2));
text(1,0.09,['p = ' num2str(p)])