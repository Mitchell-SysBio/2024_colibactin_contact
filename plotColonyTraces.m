%% this is to plot colony traces for paper figures
% set user variables 
fileName = '24h_BAC-pks_1_1600-01.czi';
myInx = 28; % scene that you want to plot signal profiles for

%% get the physical size of a pixel
reader = bfGetReader(fileName);
omeMeta = reader.getMetadataStore();
micron_per_pixel = omeMeta.getPixelsPhysicalSizeX(0).value(); % unitX = omeMeta.getPixelsPhysicalSizeX(0).unit().getSymbol()
reader.close();

%% open file
data = {}; 
data = bfopen(fileName); % Import the czi file (utility function from bfmatlab)
n = size(data,1);

%% separate channels
C1 = {}; C2 = {}; C3 = {}; C4 = {};
nameArray = {};
for i = 1:n
    C1{i} = double(data{i,1}{1})./2^16; % phase
    C2{i} = double(data{i,1}{2})./2^16; % mcherry (BAC cell)
    C3{i} = double(data{i,1}{3})./2^16; % YFP (reporter)
    C4{i} = double(data{i,1}{4})./2^16; % CFP (reporter)
end

%% draw line across reporter and tox colonies
CFPintensity = {}; YFPintensity = {}; mChintensity = {};
xCoords = {}; yCoords = {};

h = figure('color','white');
% only draw roi lines for the scene of interest in the czi file
for i = myInx
    curImg = imfuse(C1{i},C4{i},'falsecolor'); %merge mcherry and YFP channels
    imshow(curImg);
    title(['Draw a line between colony centers (purple to green), mark edges of purple then green; ',num2str(i)])
    roi = drawline %draw line to calculate intensity along (always draw from mCh to CFP)
    x = roi.Position(:,1); y = roi.Position(:,2); % pull out x and y values on the roi line
    distance = sqrt(diff(x).^2 + diff(y).^2); % get distance between points
    line([x(1),x(end)], [y(1),y(end)]);
    roi1 = drawline %mark edge of colony 1
    x1 = roi1.Position(:,1); y1 = roi1.Position(:,2);
    roi2 = drawline %mark edge of colony 2
    x2 = roi2.Position(:,1); y2 = roi2.Position(:,2);
    [xRep{i} yRep{i}] = polyxpoly([x(1),x(2)],[y(1),y(2)],[x1(1),x1(2)],[y1(1),y1(2)]); % find the intersect of the colony edge line with the roi line
    [xTox{i} yTox{i}] = polyxpoly([x(1),x(2)],[y(1),y(2)],[x2(1),x2(2)],[y2(1),y2(2)]); % find the intersect of the colony edge line with the roi line
    numPoints = ceil(sum(distance)); % Number of points along the line
    xCoords{i} = linspace(x(1), x(2), numPoints);
    yCoords{i} = linspace(y(1), y(2), numPoints);
    pause(1);
end
close(h);

%% Calculate the signal intensity along the line for each channel

for i = myInx
    DICintensity{i} = interp2(double(C1{i}), xCoords{i}, yCoords{i});
    mChintensity{i} = interp2(double(C2{i}), xCoords{i}, yCoords{i});
    YFPintensity{i} = interp2(double(C3{i}), xCoords{i}, yCoords{i});
    CFPintensity{i} = interp2(double(C4{i}), xCoords{i}, yCoords{i});   
end

%% compact all output to a single structure
col = {};
for i=myInx
    col{i}.x = xCoords{i};
    col{i}.y = yCoords{i};
    col{i}.c1 = DICintensity{i}; col{i}.c1_norm = normalize(col{i}.c1,'range');
    col{i}.c2 = mChintensity{i}; col{i}.c2_norm = normalize(col{i}.c2,'range');
    col{i}.c3 = YFPintensity{i}; col{i}.c3_norm = normalize(col{i}.c3,'range');
    col{i}.c4 = CFPintensity{i}; col{i}.c4_norm = normalize(col{i}.c4,'range');
    col{i}.micron_per_pixel = micron_per_pixel;
    col{i}.edges = [xRep{i},yRep{i}; xTox{i},yTox{i}];
end

%% get variables for plotting
i = myInx; % index to pull out for plotting
scaleFactor = double(col{i}.micron_per_pixel); % microns per pixel
dx = col{i}.x(2)-col{i}.x(1); dy = col{i}.y(2)-col{i}.y(1);
dist_unit = sqrt(dx^2+dy^2); % distance between points
totalDist = length(col{i}.x); % get length of roi line across both colonies
decayDist = [1:1:totalDist]*dist_unit*scaleFactor; % convert roi line to microns for plotting

figure; hold on;
plot(decayDist,normalize(col{i}.c2,'range'),'r','LineWidth',3);
plot(decayDist,normalize(col{i}.c3,'range'),'Color',[1 0.8 0.1250],'LineWidth',3);
plot(decayDist,normalize(col{i}.c4,'range'),'c','LineWidth',3);
grid on; box on;

title('touching recA decay pks+');
grid on; box on;
xlim([0 decayDist(end)])
xlabel('distance (microns)'); ylabel('normalized signal');
