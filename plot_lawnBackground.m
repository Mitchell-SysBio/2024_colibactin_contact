%% user definitions

files = dir('*.czi'); % read in all microscopy files
conds = {};
for i = 1:length(files)
    cur = files(i).name;
    conds{i} = string(extractBefore(cur,'-')); % get conditions from file name
end

conditions = replace([conds{:}],'_',' '); % replace underscore with space to correctly label conditions in plots
n = length(files); % number of files to loop through
times = extractBefore(conditions,'background'); % get time points from file names
%% quantify yfp in colonies
bkgrnd = {};
for iFile = 1:n
    bkgrnd = [bkgrnd, quant_backgroundSignal(files(iFile).name)]; % run background signal function
end

save bkgrnd.mat bkgrnd
%% plot yfp

yfpBgnd = []; cfpBgnd = [];
for iTime = 1:n
    curTime = bkgrnd{iTime};
    yfpBgnd(iTime) = curTime.YFPbackground; % extract YFP
    cfpBgnd(iTime) = curTime.CFPbackground; % extract CFP
    stdYFP(iTime) = std(curTime.medYFP); % calculate standar deviation for YFP
    stdCFP(iTime) = std(curTime.medCFP); % calculate standar deviation for CFP
end

figure; 
% plot YFP with error over time (sanity check)
subplot(2,1,1);
bar(yfpBgnd); hold on;
errorbar(yfpBgnd,stdYFP,'.k');
grid on; box on;
title('YFP backbround')
set(gca,'xtick',1:n,'xticklabel',times)

% plot CFP with error over time (sanity check)
subplot(2,1,2)
bar(cfpBgnd); hold on;
errorbar(cfpBgnd,stdCFP,'.k');
grid on; box on;
title('CFP background')
set(gca,'xtick',1:n,'xticklabel',times)
ylim([0 3.5e-3])


%% save in format for downstream analysis
backgroundSignal = {yfpBgnd, cfpBgnd, times};
save backgroundSignal.mat backgroundSignal
