function [col] = segment_recASignal(fileName,tfDebug)

if(nargin == 2)
    iSampleImage = 1; % image to plot for debug purposes
else
    tfDebug = false;
end

%% definitions
% get the physical size of a pixel
reader = bfGetReader(fileName);
omeMeta = reader.getMetadataStore();
micron_per_pixel = omeMeta.getPixelsPhysicalSizeX(0).value(); % unitX = omeMeta.getPixelsPhysicalSizeX(0).unit().getSymbol()
reader.close();

%% open file
data = {}; 
data = bfopen(fileName); % Import the czi file (utility function from bfmatlab)
n = size(data,1); % number of positions within file

%% separate channels
C1 = {}; C2 = {}; C3 = {}; C4 = {};
nameArray = {};
for i = 1:n
    C1{i} = double(data{i,1}{1})./2^16; % phase
    C2{i} = double(data{i,1}{2})./2^16; % mcherry (BAC cell)
    C3{i} = double(data{i,1}{3})./2^16; % YFP (reporter)
    C4{i} = double(data{i,1}{4})./2^16; % CFP (reporter)
end

if(tfDebug) % look at each channel as a sanity check you know which is which
    figure; hold on;
    subplot(2,2,1); hold on; imshow(imadjust(C1{iSampleImage})); title('channel #1');
    subplot(2,2,2); hold on; imshow(imadjust(C2{iSampleImage})); title('channel #2');
    subplot(2,2,3); hold on; imshow(imadjust(C3{iSampleImage})); title('channel #3');
    subplot(2,2,4); hold on; imshow(imadjust(C4{iSampleImage})); title('channel #4');
    pause(2);
end

%% draw line across reporter and tox colonies
CFPintensity = {}; YFPintensity = {}; mChintensity = {};
xCoords = {}; yCoords = {};

h = figure('color','white');
for i = 1:n
    curImg = imfuse(C1{i},C4{i},'falsecolor'); %merge mcherry and YFP channels
    imshow(curImg);
    title(['Draw a line between colony centers (purple to green), mark edges of purple then green; ',num2str(i)])
    roi = drawline %draw line to calculate intensity along (always draw from mCh to CFP)
    x = roi.Position(:,1); y = roi.Position(:,2); % get x and y to determine space between points on the drawn lines
    distance = sqrt(diff(x).^2 + diff(y).^2); % calculate distance between points
    line([x(1),x(end)], [y(1),y(end)]); % make a line
    roi1 = drawline %mark edge of colony 1
    x1 = roi1.Position(:,1); y1 = roi1.Position(:,2); % get x and y of reporter colony
    roi2 = drawline %mark edge of colony 2
    x2 = roi2.Position(:,1); y2 = roi2.Position(:,2); % get x and y of producer colony
    [xRep{i} yRep{i}] = polyxpoly([x(1),x(2)],[y(1),y(2)],[x1(1),x1(2)],[y1(1),y1(2)]); % find intersect of reporter colony line with cross-section line
    [xTox{i} yTox{i}] = polyxpoly([x(1),x(2)],[y(1),y(2)],[x2(1),x2(2)],[y2(1),y2(2)]); % find intersect of producer colony line with cross-section line
    numPoints = ceil(sum(distance)); % Number of points along the line
    xCoords{i} = linspace(x(1), x(2), numPoints); % save the x coordinates of the cross-section line
    yCoords{i} = linspace(y(1), y(2), numPoints); % save the y coordinates of the cross-section line
    pause(1);
end
close(h);

%% Calculate the signal intensity along the line for each channel

for i = 1:n 
    DICintensity{i} = interp2(double(C1{i}), xCoords{i}, yCoords{i});
    mChintensity{i} = interp2(double(C2{i}), xCoords{i}, yCoords{i});
    YFPintensity{i} = interp2(double(C3{i}), xCoords{i}, yCoords{i});
    CFPintensity{i} = interp2(double(C4{i}), xCoords{i}, yCoords{i});   
end

%% compact all output to a single structure
col = {};
for i=1:n
    col{i}.x = xCoords{i};
    col{i}.y = yCoords{i};
    col{i}.c1 = DICintensity{i}; col{i}.c1_norm = normalize(col{i}.c1,'range');
    col{i}.c2 = mChintensity{i}; col{i}.c2_norm = normalize(col{i}.c2,'range');
    col{i}.c3 = YFPintensity{i}; col{i}.c3_norm = normalize(col{i}.c3,'range');
    col{i}.c4 = CFPintensity{i}; col{i}.c4_norm = normalize(col{i}.c4,'range');
    col{i}.micron_per_pixel = micron_per_pixel;
    col{i}.edges = [xRep{i},yRep{i}; xTox{i},yTox{i}];
end

if(tfDebug) % show image cropped to the roi region and plot the fluorescent profiles
    img1 = imadjust(C1{iSampleImage}); img2 = imadjust(C3{iSampleImage});
    x = col{iSampleImage}.x; y = col{iSampleImage}.y;
    figure; 
    subplot(2,2,1); hold on; 
    imshowpair(img1,img2);
    plot(x,y,'-w','LineWidth',2); 
    title(['Fulll image (C1+C3) - pos #' num2str(iSampleImage)]);

    subplot(2,2,2); hold on;
    img1 = imcrop(img1,[min(x) min(y) range(x) range(y)]);
    img2 = imcrop(img2,[min(x) min(y) range(x) range(y)]);
    imshowpair(img1,img2);
    title(['Cropped image (C1+C3) - pos #' num2str(iSampleImage)]);

    subplot(2,2,[3 4]); hold on;
    plot(col{iSampleImage}.x,col{iSampleImage}.c2_norm,'-r');
    plot(col{iSampleImage}.x,col{iSampleImage}.c3_norm,'-k');
    plot(col{iSampleImage}.x,col{iSampleImage}.c4_norm,'-c');
    legend({'mCherry','YFP','CFP'})
    title(['Normalized intensity profiles - pos #' num2str(iSampleImage)])
end

end
