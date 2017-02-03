function [barPos, frame] = poletracker_avi(fname, radius,frames,drawing)

% helper function for extracting pole location
% Created by Mathew Evans 08.09.15
% based on code by R.Petersen 040914
%
% 
%
% AIM: use this script to develop pole tracking code and to find
% parameters, eg for the pole radius
%
% 

%close all

    
% Define some parameters

% assumed pole radius:
handles.pole.radius = radius;   % units of pixels
handles.pole.roi = [1 300 1 180];   % x,ys
    
    
% Speed up poletracking if video object has already been created
if strcmp(class(fname),'VideoReader')
    v = fname;        
else
    
    % the video:
    handles.fname = fname;
    % Set up videoReader object if none supplied
    v = VideoReader([fname,'.avi']);
end


%% loop to find pole
barPos = zeros(2,length(frames));
i = 0;
for frameidx = frames; %video.header.nframes
    i = i+1;
    
    %frameidx
    
    % Load a frame of data
    clear frame
    frame = read(v,frameidx);
    frame = squeeze(frame(:,:,1));
    handles.frame = frame;
    
    if drawing;
        clf;
        imagesc(handles.frame)
        axis image
        colormap gray
        
        %     keyboard
        
        % Re-scale figure brightness to 0:0.995 of range
        n = hist(handles.frame(:),0:255);
        ncut = find(cumsum(n)/sum(n)>0.995,1);
        caxis([0 ncut])
        clear n ncut
    end
    
    % Find pole centre based on convolution of a pole-sized disk across the
    % image
    % gof = goodness of fit
    [polecentre, gof] = polefit(handles.frame,handles.pole.radius,handles.pole.roi);
    
    if drawing
        hold on
        
        plot_circle(polecentre,handles.pole.radius)
        title(sprintf('frame %d: pc = (%d,%d) gof=%.1f', frameidx, polecentre,gof))
        drawnow
%         pause
    end
    
    barPos(:,i) = polecentre;

    %     clf
end

end

%% Subfunctions

function plot_circle(centre,radius)

dtheta = 2*pi/100;
theta = 0:dtheta:2*pi-dtheta;

x = centre(1) + radius*cos(theta);
y = centre(2) + radius*sin(theta);

plot(x,y,'r:',centre(1),centre(2),'rx')
% plot(centre(1),centre(2),'rx')
end

%

function [polecentre, gof] = polefit(frame,radius,roi)
% TO DO:
%       find pole of arbitrary size - maybe use a template with a gradiant?
%       otherwise, try a range of sizes and output the best

% ALSO. Take whisker contact into account. Use whisker theta to find
% contact side, also previous pole position to know whether contact is
% expected.

frame = double(frame(:,:,1));
frame = frame - mean(frame(:));                         % zero mean
frame = frame([roi(3):roi(4)],[roi(1):roi(2)]);         % Restrict search to ROI
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
[X,Y] = meshgrid([1:size(z,2)],[1:size(z,1)]);
X = X(:); Y = Y(:);
polecentre(1) = X(midx) + (length(kernel)-1)/2;
polecentre(2) = Y(midx) + (length(kernel)-1)/2;

end