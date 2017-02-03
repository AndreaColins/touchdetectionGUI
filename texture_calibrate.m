function [C,image] = texture_calibrate(fname,frames,drawing,pausing,roi1,roi2)

% texture_calibrate
%
% Generate calibration matrix C of pixel locations x,y,v,w and X,Y,Z real
% positions in mm (from motor commands and known offset geometry 
% (from measuring the calibration widget and from calibration log file) 
% 
% If point 1 is at 0,0,0mm, point 2 (higher, with both sides pointing out, on the right/more anticlockwise
% in the top view) is at 0,0.6,2mm
%
% Initial look is that: 
% each pixel in top-down view is ~0.067mm,
% each pixel in front-on view is ~0.051mm.
%
% M.Evans 18.11.15. Initially developed by Mathew Evans 03.03.2015 as
% poletracker_3D.m
%
% This script loads a file then find the pole position in each view from
% 3D video data, then finds 2 needle points of calibration widget
%
% Expects certain inputs and outputs
% fname = file name string in current directory, including extension
% frames = frames to track, can be just one or a range
% drawing/pausing = switches for plotting/pausing each frame
% roi1,roi2, top down and mirror view rois
% Typical values: [320 700 10 250],[50 250 100 262]
%
% Outputs:
% C = calibration matrix, where each row is x,y,v,w, X,Y,Z
% image = image array of tracked frame. Useful for plotting

%

%close all

% NOTE: use  cd I:\Sarah Fox\3D_testing
%            fname = 'grid8_23_Y29_20150226_025204 PM.dat'
%            frame = 6800:35000
%
% run the code with
% for i = frame;
%     [C, image] = texture_calibrate(fname,i,0,0);
% end
% save poletracking.mat barPos notchPos

% assumed pole radius:
handles.pole.radius = 40;   % units of pixels
handles.pole.roi = roi1;    % x x, y y . [300 600 50 250];
handles.mirror.roi =  roi2;  % x x , y y . [50 250 225 375]; % [40 210 210 420];

% the video:
handles.fname = fname;

switch handles.fname(end-2:end)
    case 'avi'
        video.type = 'avi';
        video.vObj = VideoReader(handles.fname);
        handles.nframes = video.vObj.NumberOfFrames;
        video.width = video.vObj.Width;
        video.height = video.vObj.Height;
    case 'dat'
        video.type = 'dat';
        video.fname = handles.fname;
        video.fid = fopen(handles.fname,'r');
        video.header = read_mikrotron_datfile_header(video.fid);
        handles.nframes = video.header.nframes;
        video.width = video.header.width;
        video.height = video.header.height;
        video.offset = 8192;
        % set file position to start of first frame
        fseek(video.fid,8192,-1);
        
    otherwise
        fprintf('Please choose a video file to display instead of an analysis output file \n \n')
        msg = ['Unhandled file type: ',handles.fname(end-3:end)];
        error(msg)
        
end
handles.video = video


%% loop to find pole/calibration points
barPos = zeros(2,length(frames));
C = zeros(2*length(frames),7); % 2* numel frames for 2 calibration points
i = 0;
for frameidx = frames; %video.header.nframes
    i = i+1;
    
    frameidx
    
    % Load a frame of data
    handles.frame = load_frame(video,frameidx);
    handles.framemean = mean(handles.frame(:));
    image = double(handles.frame(:,:,1));
    if drawing;
        clf;
        subplot(3,1,[1,2]);
        imagesc(image)
        hold all
        % Draw rois
%         plot([roi1(1),roi1(2),roi1(2),roi1(1),roi1(1)],[roi1(3),roi1(3),roi1(4),roi1(4),roi1(3)]);
%         plot([roi2(1),roi2(2),roi2(2),roi2(1),roi2(1)],[roi2(3),roi2(3),roi2(4),roi2(4),roi2(3)]);
     %   axis image
        colormap gray
        
%             keyboard
        
        % Re-scale figure brightness to 0:0.995 of range
        n = hist(image(:),0:255);
        ncut = find(cumsum(n)/sum(n)>0.995,1);
        caxis([0 ncut])
        clear n ncut
    end
    
    % Find pole centre based on convolution of a pole-sized disk across the
    % image
    % gof = goodness of fit
    [polecentre, pole_gof] = polefit(handles.frame,handles.pole.radius,handles.pole.roi);
    
    % Find calibration points in top-down view. Uses inferred pole location to restrict search
    % Outputs concatenated x,y of both points
    handles.notch_roi = [polecentre(1)-60,polecentre(1)+10, polecentre(2),polecentre(2)+80];
    
    [notchPos] = notchfit(handles.frame,handles.notch_roi);
    

%     notchPos(1) = notchPos(1) + handles.pole.roi(1) - 1;
%     notchPos(2) = notchPos(2) + handles.pole.roi(3) - 1;
%     notchPos(3) = notchPos(3) + handles.pole.roi(1) - 1;
%     notchPos(4) = notchPos(4) + handles.pole.roi(3) - 1;

    
    C(i,1) = notchPos(1);
    C(i,2) = notchPos(2);
    C(i + length(frames),1) = notchPos(3);
    C(i + length(frames),2) = notchPos(4);
    
    % Find calibration points in mirror view. Outputs concatenated x,y of both points
   [notchPos] = notchfit_m(handles.frame,handles.mirror.roi);
    
    C(i,3) = notchPos(1);
    C(i,4) = notchPos(2);
    C(i + length(frames),3) = notchPos(3);
    C(i + length(frames),4) = notchPos(4);

    if drawing
        subplot(3,2,[1:4]);
        imagesc(handles.frame); hold all;
        plot_circle(polecentre,handles.pole.radius);
        plot_calib([C(i,1:4);C(i+length(frames),1:4)]);
        title(sprintf('frame %d: pc = (%d,%d) gof=%.1f', frameidx, polecentre,pole_gof))
        plot([roi1(1),roi1(2),roi1(2),roi1(1),roi1(1)],[roi1(3),roi1(3),roi1(4),roi1(4),roi1(3)]);
        plot([roi2(1),roi2(2),roi2(2),roi2(1),roi2(1)],[roi2(3),roi2(3),roi2(4),roi2(4),roi2(3)]);
        drawnow
        if pausing
            pause
        end
    end
    
    barPos(:,i) = polecentre;
    
    %     clf

%     pause
end

% Set last 3 columns of C to x,y and z
C(:,5) = C(:,1);
C(:,6) = C(:,2);
C(1:length(frames),7) = 0;
C(length(frames)+1:2*length(frames),7) = 29; % 2mm @ 0.0687 mm/pixel


figure
plot3(C(1:length(frames),5),C(1:length(frames),6),C(1:length(frames),3),'.')
hold all;
plot3(C(length(frames)+1:2*length(frames),5),C(length(frames)+1:2*length(frames),6),C(length(frames)+1:2*length(frames),3),'.')
end




%% Subfunctions

function frame = load_frame(videopmtrs,frameidx)
switch videopmtrs.type
    case 'avi'
        video = read(videopmtrs.vObj, frameidx);
        frame = video(:,:,:,1);
    case 'dat'
        videopmtrs.fid = fopen(videopmtrs.fname,'r');
        offset = videopmtrs.header.imagesize * (frameidx-1) + videopmtrs.offset;
        fseek(videopmtrs.fid,offset,-1);
        tmp = fread(videopmtrs.fid,videopmtrs.header.imagesize-24,'uint8=>uint8');
        tmp = reshape([tmp; zeros(24,1)],videopmtrs.width,videopmtrs.height)';
        frame = uint8(zeros(videopmtrs.height,videopmtrs.width,3));
        frame(:,:,1) = tmp;
        frame(:,:,2) = tmp;
        frame(:,:,3) = tmp;
        clear tmp
        fclose('all') ; % Trying to fix 'too many open files' error
    otherwise
        error('Unhandled video file type')
end


end

%

function plot_circle(centre,radius)

dtheta = 2*pi/100;
theta = 0:dtheta:2*pi-dtheta;

x = centre(1) + radius*cos(theta);
y = centre(2) + radius*sin(theta);

plot(x,y,'r:',centre(1),centre(2),'rx')
%plot(notchPos(1),notchPos(2),'c+',notchPos(3),notchPos(4),'c+');
% plot(centre(1),centre(2),'rx')
end

% Plot calibration points. Expects 2 rows of x,y,v,w
function plot_calib(C)
plot(C(1,1),C(1,2),'r.')
plot(C(1,3),C(1,4),'r.')
plot(C(2,1),C(2,2),'b.')
plot(C(2,3),C(2,4),'b.')
end

function [polecentre, gof] = polefit(frame,radius,roi)
% TO DO:
%       find pole of arbitrary size - maybe use a template with a gradiant?
%       otherwise, try a range of sizes and output the best

% ALSO. Take whisker contact into account. Use whisker theta to find
% contact side, also previous pole position to know whether contact is
% expected.
%       Use pacman
% keyboard
frame = double(frame(:,:,1));
frame = frame - mean(frame(:));                         % zero mean
frame = frame(roi(3):roi(4),roi(1):roi(2));         % Restrict search to ROI

% Zero padding to help with edge effects
frame = [zeros(50,size(frame,2)); frame];


[X,Y] = meshgrid(-radius-2:radius+2,-radius-2:radius+2);
kernel = ones(2*radius+5);
% kernel(1:2*radius,1:2*radius) = 0;    % pacman
idx = X.^2 + Y.^2 <= radius^2;
alpha = (sum(kernel(:))-sum(idx(:)))/sum(idx(:));   % pacman
% alpha = (sum(kernel(:))-.75*sum(idx(:)))/sum(idx(:));   % pacman
kernel(idx) = -alpha;

% clear idx
z = conv2(frame,kernel,'valid');
[gof, midx] = max(z(:)/sum(idx(:)));
[X,Y] = meshgrid(1:size(z,2),1:size(z,1));
X = X(:); Y = Y(:);
polecentre(1) = X(midx) + (length(kernel)-1)/2;
polecentre(2) = Y(midx) + (length(kernel)-1)/2;

polecentre(1) = polecentre(1) + roi(1) - 1;
polecentre(2) = polecentre(2) + roi(3) - 1 - 50;

% % To debug uncomment next bit
% clf
% subplot(3,1,[1,2]);
% imagesc(frame)
% colormap gray
% hold all
% plot(X(midx)+ (length(kernel)-1)/2,Y(midx)+ (length(kernel)-1)/2,'r.');
% subplot(3,1,3);
% plot(min(z));
% drawnow;
% keyboard
end

% FIND CALIBRATION WIDGET SPIKES
% IDEA 1 : find peaks in shadow near the pole
% IDEA 2 (untested) find 2 pacman negatives (best and second best) in region near
% identified pole center
function [notchPos] = notchfit(frame,roi)
% V1 find edge of pole silhouette (including calibration widget), then find
% peaks corresponding to two points
% V1 works only at pixel resolution

% keyboard
frame = double(frame(:,:,1));
frame = frame(roi(3):roi(4),roi(1):roi(2));         % Restrict search to ROI
frame = frame - mean(frame(:));                         % zero mean

% Make first row of frame dark to help findpeaks later
frame(1,:) = 0.5*min(frame(:))*ones(1,size(frame,2));

% keyboard
% Smooth image
g = fspecial('gaussian',[5 5],.5); % filter
im = imfilter(frame,g);            % smooth image
% im = frame;

e1 = edge(im); % find edges of pole silhouette
[i,j] = find(e1);

% Work out location of edge contour
[x,jx] = unique(j,'last');
% Find peaks corresponding to calibration points
[p,pi] = findpeaks(i(jx),'NPeaks',2,'MinPeakDistance',6,'MinPeakHeight',50); % Only find 2 please. May need to constrain other ways...
% keyboard
% If 2 peaks are not found by the first criteria, relax the criteria.
if numel(pi)~= 2;
    display('Peaks not found in top-down view, relaxing criteria')
    [p,pi] = findpeaks(i(jx),'NPeaks',2);
end
%% plotting for debugging
clf
subplot(3,2,6)
imagesc(im); colormap gray
hold all;
plot(j(jx),i(jx))
plot(j(jx(pi)),i(jx(pi)),'r.');
% hold off;

%%
notchcentre(1) = j(jx(pi(1))) + roi(1) - 1;
notchcentre(2) = i(jx(pi(1))) + roi(3) - 1;
notchcentre(3) = j(jx(pi(2))) + roi(1) - 1;
notchcentre(4) = i(jx(pi(2))) + roi(3) - 1;

notchPos = notchcentre;

% Uncomment for debugging
% clf
% imagesc(frame)
% colormap gray
% hold all
% plot(X(midx) + (size(notch_temp,1)-1)/2,Y(midx) + (size(notch_temp,2)-1)/2,'b.');
% drawnow;
% keyboard

% reorder such that notch 1 is always on the left
if notchcentre(1) > notchcentre(3);
    notchPos(3) = notchcentre(1);
    notchPos(4) = notchcentre(2);
    notchPos(1) = notchcentre(3);
    notchPos(2) = notchcentre(4);
end
end

% FIND CALIBRATION WIDGET SPIKES in MIRROR
% IDEA 1 : find peaks in shadow near the pole (pole needs to be found here in 
% this view
% IDEA 2 (untested) find 2 pacman negatives (best and second best) in region near
% identified pole center
function [notchPos] = notchfit_m(frame,roi)
% V1 find edge of pole silhouette (including calibration widget), then find
% peaks corresponding to two points
% V1 works only at pixel resolution

frame = double(frame(:,:,1));
% keyboard
[m,midx] = min(sum(frame(roi(3):roi(4),roi(1):roi(2))'));

% recompute roi
% FIRST. Compute gradient of pole object to set second ROI.

roi(3) = roi(3) + midx; %+ 25; % Isolating to only the spikes by shifting roi down a load of steps. 
                           % May be brittle

im = frame(roi(3):roi(4),roi(1):roi(2));         % Restrict search to ROI
im = im - mean(im(:));                           % zero mean

% keyboard
% Smooth image
g = fspecial('gaussian',[5 5],.5); % filter
im = imfilter(im,g);            % smooth image



% Simple V2, find edges, then keep the first value for each x axis location
e1 = edge(im); % find edges of pole silhouette
[i,j] = find(e1);
% Work out location of edge contour
[x,jx] = unique(j,'first'); % enforce first as different matlabs have different defaults
track_contour = i(jx);


% [dx,dy] = gradient(im);          % Find vertical gradient
% % More involved method
% dy(find(dy<0.7*max(dy(:)))) = 0; % Set low values to a single number in
% an attempt to eliminate false peaks
% % Find max of unique(find(dy)) and re-set roi
% [i,j] = find(dy);
% [x,jx] = unique(j);
% roi(3) = roi(3) + max(i(jx))+2;
% % go again
% im = im(max(i(jx))+2:end,1:end);
% [dx,dy] = gradient(im);          % Find vertical gradient
% % dz = dx + dy;
% dy(find(dy<0.7*max(dy(:)))) = 0; % Set low values to a single number

% Simple method skips to here: just finds the max of dy
% [~,track_contour] = max(dy);     % find max in contour to find spikes

% Find peaks corresponding to calibration points
[p,pi] = findpeaks(track_contour,'NPeaks',2,'MinPeakDistance',25,'MinPeakHeight',5); % Only find 2 please. May need to constrain other ways...

if numel(pi)~= 2;
    display('Peaks not found in mirror view, relaxing criteria')
    [p,pi] = findpeaks(i(jx),'NPeaks',2);
end

% keyboard
%% plotting for debugging
subplot(3,2,5)
% imagesc(dy);
imagesc(im); colormap gray
hold all;
plot(j(jx),track_contour)
plot(j(jx(1))+pi-1,p,'r.');
hold off;

%%
notchcentre(1) = j(jx(1)) +pi(1) + roi(1) - 1;
notchcentre(2) = p(1) + roi(3) - 1;
notchcentre(3) = j(jx(1)) + pi(2) + roi(1) - 1;
notchcentre(4) = p(2) + roi(3) - 1;

notchPos = notchcentre;

% Uncomment for debugging
% clf
% imagesc(frame)
% colormap gray
% hold all
% plot(X(midx) + (size(notch_temp,1)-1)/2,Y(midx) + (size(notch_temp,2)-1)/2,'b.');
% drawnow;
% keyboard

% reorder such that notch 1 is always on the right
if notchcentre(3) > notchcentre(1);
    notchPos(3) = notchcentre(1);
    notchPos(4) = notchcentre(2);
    notchPos(1) = notchcentre(3);
    notchPos(2) = notchcentre(4);
end
end


















% TO DO: UPDATE TO FIND CALIBRATION WIDGET SPIKES
function [notchPos] = notchfit2(frame,roi,drawing)
% Find notches using pole gradient, and deviations from this
% V1 find pole gradient, extract 10 pixels above and below, then find dips
% in this array


frame = double(frame(:,:,1));
frame = frame(roi(3):roi(4),roi(1):roi(2));         % Restrict search to ROI
frame = frame - mean(frame(:));                         % zero mean

keyboard
% Smooth image
g = fspecial('gaussian',[5 5],.5); % filter
im = imfilter(frame,g);            % smooth image

e1 = edge(im); % find edges of pole silhouette
[i,j] = find(e1);

% Work out location of edge contour
[x,jx] = unique(j);

% Find peaks corresponding to calibration points
[p,pi] = findpeaks(i(jx),'NPeaks',2); % Only find 2 please. May need to constrain other ways...

%% plotting for debugging
clf;
imagesc(im); colormap gray
hold all;
plot(j(jx),i(jx))
plot(j(jx(pi)),i(jx(pi)),'.');
hold off;

[dx, dy] = gradient(im);
grad = dx+dy;
%%



[~,pole_contour] = min(grad);
pole_contour = conv(pole_contour,gausswin(10)/sum(gausswin(10)),'same');
x_roi = 1:150;
% pole_contour = pole_contour(x_roi); % trim ends off
p = polyfit(x_roi(40:120),pole_contour(40:120),1); % extrapolate from middle segment as it's cleaner
pv = polyval(p,x_roi);
pv = round(pv);

n = 0;

for i = 15:length(x_roi)-5; % watch out for edges
    n = n+1;
% for i = 1:length(x_roi);
    track_roi(:,n) = im(pv(i)-5:pv(i)+5,x_roi(i));
end
% 
% imagesc(track_roi)
[dx, dy] = gradient(track_roi);
[~,track_contour] = min(dx+dy);

track_contour = conv(track_contour,gausswin(4)/sum(gausswin(4)),'same');

% find peaks
[pks, locs] = findpeaks(track_contour,'minpeakdistance',40,'NPeaks',2,'minpeakheight',6.5);

 if drawing;
        hold all;
        subplot(3,1,3);
        imagesc(track_roi)
        colormap gray
        hold all;
        plot(track_contour,'c');
        plot(locs,pks,'r.');
 end

% keyboard
notchcentre(1) = x_roi(1) + locs(1)+15 - 1;
notchcentre(2) = pv(locs(1)+15) + pks(1) - 6;

notchcentre(3) = x_roi(1) + locs(2)+15 - 1;
notchcentre(4) = pv(locs(2)+15) + pks(2) - 6;


notchcentre(1) = notchcentre(1) + roi(1) - 1;
notchcentre(2) = notchcentre(2) + roi(3) - 1;
notchcentre(3) = notchcentre(3) + roi(1) - 1;
notchcentre(4) = notchcentre(4) + roi(3) - 1;


notchPos = notchcentre;

% Uncomment for debugging
% clf
% imagesc(frame)
% colormap gray
% hold all
% plot(X(midx) + (size(notch_temp,1)-1)/2,Y(midx) + (size(notch_temp,2)-1)/2,'b.');
% drawnow;
% keyboard


% reorder such that notch 1 is always on the left
if notchcentre(1) > notchcentre(3);
    notchPos(3) = notchcentre(1);
    notchPos(4) = notchcentre(2);
    notchPos(1) = notchcentre(3);
    notchPos(2) = notchcentre(4);
end

% CODE BELOW WAS NO BETTER. WORSE IN FACT
% % for each point along pv, work out 5 locations perpendicular to sample for
% % generating track_roi (instead of 10 pixels in y axis).
% 
% p_2 = -1/p(1); % slope of perpendicular line
% 
% keyboard
% 
% % extract points 5 steps either side of pv
% n = 0;
% for i = 6:length(x_roi)-5; % watch out for edges
%     n = n+1;
%     xrange = i-5:i+5;
%     yrange = round((p_2*xrange) + (pv(i)-p_2*i)); % work out points along line
%     yrange = round(yrange);
%     indices = [yrange;xrange]';
%     for j = 1:11;
%     track_roi2(j,n) = im(yrange(j),xrange(j));
%     end
% end

end

