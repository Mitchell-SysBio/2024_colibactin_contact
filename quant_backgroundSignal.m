function [fluor] = quant_backgroundSignal(fileName,tfDebug)

if(nargin == 2)
    iSampleImage = 1; % image to plot for debug purposes
else
    tfDebug = false;
end

%% open file
data = {}; 
data = bfopen(fileName); % Import the czi file (utility function from bfmatlab)
n = size(data,1);
C1 = {}; C2 = {}; C3 = {}; C4 = {};
for i = 1:n
    C1{i} = double(data{i,1}{1})./2^16; % phase
    C2{i} = double(data{i,1}{2})./2^16; % mcherry (producer cell)
    C3{i} = double(data{i,1}{3})./2^16; % YFP (reporter)
    C4{i} = double(data{i,1}{4})./2^16; % CFP (reporter)
end

%% draw linear roi across the lawn to calculate YFP and CFP signal intensity 
CFPintensity = {}; YFPintensity = {};
xCoords = {}; yCoords = {};

h = figure;
for i = 1:n
    curImg = C1{i};
    imshow(imadjust(curImg)); %open phase image
    title('draw line on the lawn')
    roi = drawline %draw roi line
    x = roi.Position(:,1); y = roi.Position(:,2); % get x and y positions
    distance = sqrt(diff(x).^2 + diff(y).^2); % get length of line
    numPoints = ceil(sum(distance)); % Number of points along the line
    xCoords{i} = linspace(x(1), x(2), numPoints); % get x coordinates
    yCoords{i} = linspace(y(1), y(2), numPoints); % get y coordinates
    pause(1);
end
close(h);

%% get signal intensity along linear roi for YFP and CFP

for i = 1:n 
    YFPintensity{i} = interp2(double(C2{i}), xCoords{i}, yCoords{i}); %interpolate values along the line
    CFPintensity{i} = interp2(double(C3{i}), xCoords{i}, yCoords{i});   
end

%% take median signal per image and average across all images
medYFP = []; medCFP = [];
for i = 1:n
    medYFP(i) = nanmedian(YFPintensity{i}); %take median signal intensity for all pixels per image
    medCFP(i) = nanmedian(CFPintensity{i});
end

YFPbgnd = mean(medYFP); % take mean of all images
CFPbgnd = mean(medCFP);
%% save variables
fluor = {};
fluor.medYFP = medYFP;
fluor.medCFP = medCFP;
fluor.YFPbackground = YFPbgnd;
fluor.CFPbackground = CFPbgnd;

end

