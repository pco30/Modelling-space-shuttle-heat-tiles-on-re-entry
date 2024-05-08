function [timeData, tempData] = plottempAutoo(name,varargin)
%PLOTTEMPAUTOO Automatically scan temperture data
% 
%   [TIMEDATA, TEMPDATA] = PLOTTEMPAUTOO(NAME,VARARGIN) analyses images of
%   the temperature data NAME and outputs the temperature, TEMPDATA, in
%   kelvin and time, TIMEDATA, in seconds. VARARGIN produces useful plots.

% Read image digitally
img = imread([name '.jpg']);

% Create balck and white image.
bwimg = img(:,:,1) < 205 | img(:,:,2) < 205 | img(:,:,3) < 205;

% Collect axis data:
xaxis = ones(1,2);
yaxis = ones(1,2);

% Determine the row and collumn in image with the most number of white
% cells and store it in the first position of variables 'xaxis' and
% 'yaxis'.
for i = 1:length(bwimg(:,1))
    if sum(bwimg(i,:)) > sum(bwimg(xaxis(1),:))
        xaxis(1) = i;
    end
end

for i = 1:length(bwimg(1,:))
    if sum(bwimg(:,i,1)) > sum(bwimg(:,yaxis(1)))
        yaxis(1) = i;
    end
end

% Find the length of the X and Y axis and store it in the second
% position of variable 'xaxis' and 'yaxis' respectiely.
axisLen = diff(find(diff(bwimg(xaxis(1),:)) ~= 0));
xaxis(2) = max(axisLen(1:2:end));

axisLen = diff(find(diff(bwimg(:,yaxis(1))) ~= 0));
yaxis(2) = max(axisLen(1:2:end));

% Crop the image
cropimg = img(1:xaxis(1),yaxis(1):end,:);

% extract red plot.
data = cropimg(:,:,1) > 140 & cropimg(:,:,2) < 105;

rawtempData = [];
rawtimeData = [];
xscale = 2000;
yscale = 2000;

% Collect temperature data by taking the average of the red plot at each
% pixel along the x axis.
for i = 1:length(data(1,:))
    if any(data(:,i))
        rawtempData = [rawtempData mean(find(data(:,i)))];
        rawtimeData = [rawtimeData i];
    end
end

% Data could be reduced to a fixed number of points for the same time
% range however correct time values are interpolated later on so there is
% no need to reduce number of points.


% Scale the points to match the correct plot values and convert temperature
% from Farenheit to Kelvin
tempData = (((xaxis(1) - rawtempData) * (yscale/yaxis(2))) - 32) * 5/9 + 273.25;

% tempData = ((xaxis(1) - tempData) * (yscale/yaxis(2)));
timeData = ((rawtimeData) * (xscale/xaxis(2)));

% Duplicate the end values so extrapolation produces straight lines.
tempData = [tempData(1) tempData tempData(end)];
timeData = [timeData(1)-1 timeData timeData(end)+1];

% save data to .mat file with same name as image file
save(name, 'timeData', 'tempData')


% Produces plots For troubleshooting
if nargin > 1
    figure
    plot(timeData,tempData, 'r', LineWidth=2)
    title(name)
    xlabel('Time (s)')
    ylabel('Temperature (K)')
    ylim([0 2000])
    % imshow(img(:,:,1) > 140 & img(:,:,2) < 105)
    if cell2mat(varargin) > 0
        % Plot image with automatically scanned points
        figure
        image(img)
        hold on
        plot(rawtimeData + yaxis(1)-1, rawtempData,'o',Color='g')
    end
    if cell2mat(varargin) > 1
        % Plot black and white data
        figure
        imshow(bwimg)
        axis on
    end
end