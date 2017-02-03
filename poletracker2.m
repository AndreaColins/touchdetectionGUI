function [barPos,image,radius,best_gof] = poletracker2(fname,frames,drawing)

% Updated helper function for extracting pole location, this time estimating 
% pole appearance size too (not taking radius as an input)
% 
% Developed by Mathew Evans 09.09.2015
%
% Original version created by RSP 040914
%
% 
%


% the video:
handles.fname = fname;

% assumed pole radius:
% handles.pole.radius = radius;   % units of pixels
handles.pole.roi = [1 300 1 180];   % x,ys

% get header information from the video file:
video.fid = fopen(handles.fname,'r');
video.header = read_mikrotron_datfile_header(video.fid);
handles.nframes = video.header.nframes;
video.width = video.header.width;
video.height = video.header.height;
video.offset = 8192;

% set file position to start of first frame
fseek(video.fid,8192,-1);
% keyboard
% figure;

%% loop to find pole
barPos = zeros(2,length(frames));
i = 0;
for frameidx = frames; %video.header.nframes
    i = i+1;
    
    %frameidx
    
    % Load a frame of data
    handles.frame = load_frame(video,frameidx);
    handles.framemean = mean(handles.frame(:));
    image = double(handles.frame(:,:,1));
    if drawing;
        clf;
        imagesc(image)
        axis image
        colormap gray
        
        %     keyboard
        
        % Re-scale figure brightness to 0:0.995 of range
        n = hist(image(:),0:255);
        ncut = find(cumsum(n)/sum(n)>0.995,1);
        caxis([0 ncut])
        clear n ncut
    end
    
    % Find pole centre based on convolution of a pole-sized disk across the
    % image
    % gof = goodness of fit
    [polecentre, gof, radius] = polefit(handles.frame,handles.pole.roi);
    
    [best_gof,best_r] = max(gof);
    radius = radius(best_r);
    polecentre = polecentre(best_r,:);
    
    if drawing
        hold on
        
        plot_circle(polecentre,radius)
        title(sprintf('frame %d: pc = (%d,%d) gof = %.1f, radius = %d', frameidx, polecentre,best_gof,radius))
        drawnow
%         pause
    end
    
    barPos(:,i) = polecentre;

    %     clf
end

end

%% Subfunctions

function frame = load_frame(videopmtrs,frameidx)

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
end

%

function plot_circle(centre,radius)

dtheta = 2*pi/100;
theta = 0:dtheta:2*pi-dtheta;

x = centre(1) + radius*cos(theta);
y = centre(2) + radius*sin(theta);

plot(x,y,'r:',centre(1),centre(2),'rx')
% plot(centre(1),centre(2),'rx')
end

%

function [polecentre,gof,radius] = polefit(frame,roi)
% Find pole of arbitrary size -
%       Trying a range of sizes and output the best

% ALSO. Take whisker contact into account. Use whisker theta to find
% contact side, also previous pole position to know whether contact is
% expected.

frame = double(frame(:,:,1));
frame = frame - mean(frame(:));                         % zero mean
frame = frame([roi(3):roi(4)],[roi(1):roi(2)]);         % Restrict search to ROI

radius = 12:15;

for i = 1:4;
[X,Y] = meshgrid(-radius(i)-2:radius(i)+2,-radius(i)-2:radius(i)+2);
kernel = ones(2*radius(i)+5);
% kernel(1:2*radius,1:2*radius) = 0;    % pacman
idx = X.^2 + Y.^2 <= radius(i)^2;
alpha = (sum(kernel(:))-sum(idx(:)))/sum(idx(:));   % pacman
% alpha = (sum(kernel(:))-.75*sum(idx(:)))/sum(idx(:));   % pacman
kernel(idx) = -alpha;

% clear idx
z = conv2(frame,kernel,'valid');
[gof(i), midx] = max(z(:)/sum(idx(:)));
[X,Y] = meshgrid([1:size(z,2)],[1:size(z,1)]);
X = X(:); Y = Y(:);
polecentre(i,1) = X(midx) + (length(kernel)-1)/2;
polecentre(i,2) = Y(midx) + (length(kernel)-1)/2;
end



end