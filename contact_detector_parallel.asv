function [pole_dist,gof,barPos] = contact_detector_parallel(fname,radius)
% contact_detector_parallel.m
%
% Script to detect whisker-pole contacts for an .avi video that has already
% been tracked using parfor loops to speed things up. This means you cannot
% plot as you go. For this use contact_detector.m
%
% This script finds the pole, then determines the distance between a
% tracked whisker and the pole.
%
% Inputs are a filename (without extension), and a guess at the pole
% radius (usually 14).

time0 = tic;
% Set some parameters
trial = str2double(fname(end-9:end-4)); % Extract trial number

ext_range = 100; % Distance to extrapolate tracked whisker segment, in units of original segment length

% Uncomment/ edit the next section when saving to a local folder
local_dir = 'C:\Users\Dario\Desktop\Mat_touch_gui\sticazzi'; % local, editable folder on your local machine
% local_dir = 'D:evans_work\';
remote_dir = pwd;
disp(['Creating local copy of ',fname])
old_fname = fname;
% UNIX
% fname = [local_dir,'/',old_fname];
% copyfile([remote_dir,'/',old_fname],[fname]);
% copyfile([remote_dir,'/',old_fname(1:end-4),'_clean.mat'],[fname(1:end-4),'_clean.mat']);

% WINDOWS
fname = [local_dir,'\',old_fname];
copyfile([remote_dir,'\',old_fname],[fname]);
copyfile([remote_dir,'\',old_fname(1:end-4),'_clean.mat'],[fname(1:end-4),'_clean.mat']);

disp(['done copying'])

% Set up video object

%% Load old tracking solution for this video. Using r_base for now

% SWITCH FOR .dat or .tr FILES
switch fname(end-2:end)
    case 'avi'
        load([fname(1:end-4),'_clean.mat'],'r_base','best_m');
        
    case 'dat'
        if exist([fname(1:end-4),'.tr'])
            load([fname(1:end-4),'.tr'],'r_all','-mat');
            r_base = r_all;
        else
            load([fname(1:end-4),'_clean.mat'],'r_base','best_m');
        end
        
end


% Initial version of chimera did not fill in dropped frames for r_base, so
% need to do it here sometimes
if size(r_base,1) < size(best_m,2)
    x = find(best_m);
    d = 1:numel(best_m);
    d(x) = 0;
    dropped_frames = find(d);
    
    r_chi = zeros(numel(best_m),2,3);
    r_chi(x,:,:) = r_base;
    
    dropped_frames = find(r_chi(:,1,1)==0);
    last_good_frame = find(r_chi(:,1,1)~=0);
    last_good_frame = last_good_frame(end);
    disp(['Filling in ',num2str(numel(dropped_frames)),' dropped frames']);
    for i = 1:numel(dropped_frames);
        j = dropped_frames(i);
        if j == 1;
            r_chi(j,:,:) = r_chi(last_good_frame,:,:);
        else
            r_chi(j,:,:) = r_chi(j-1,:,:);
        end
    end
    
    r_base = r_chi;
    
    % Save the newly expanded r_base
    save([fname(1:end-4),'_clean.mat'],'r_base','-append');
    copyfile([fname(1:end-4),'_clean.mat'],[remote_dir,'/',old_fname(1:end-4),'_clean.mat']);
    
end

%% Determine frames to track (only when pole is up, specified by excel file)
disp(['Loading excel file metadata'])
x_files = dir('good_trials.x*');
if isempty(x_files)
    error('No good_trials excel file found')
else
    
    if ismac
        if exist('good_trials.xls','file');
            xlsfile = 'good_trials.xls';
        else
            error('xlsx file must be converted to Excel 98 xls file when working on Mac OSX')
        end
    else
        xlsfile = 'good_trials.xlsx'; % Must be saved as .xls Excel 98 format for compatability
        
    end
end
xls_info = xlsread(xlsfile);
file_numbers = xls_info(:,1);
this_index = find(file_numbers==trial); % Find row corresponding to this trial
start_frame = xls_info(this_index,3);
if start_frame == 0
    start_frame = 1;
end
trigger_frame = xls_info(this_index,6);

% Work out frames to track. All but first few before pole comes up (can't
% determine pole down easily at the moment)
Nframes = size(r_base,1);
frames = 1:Nframes;
frames(start_frame:trigger_frame) = [];


pole_dist = zeros(1,Nframes);
barPos = zeros(Nframes,2);
gof = zeros(1,Nframes);
closest_w = zeros(Nframes,2);

%% Make a local copy of the .dat file to speed up pole tracking
% TO DO: ADD LOCAL COPY CREATION AS AN OPTION
% copyfile([fname,'.dat'],[local_dir,'/',fname,'.dat']);


%% FOR A GIVEN FRAME (MAY BE POSSIBLE TO DO ALL FRAMES AT ONCE EVENTUALLY)

time1 = tic;
parfor I = 1:numel(frames); % Can be changed to parfor for speed up (but don't plot anything)
    i = frames(I);
    if mod(i,100)==0
        disp(['Finding pole in frame ',num2str(i)])
    end
    %% Find pole
    
%     [pole_location,~,gof(I)] = poletracker([local_dir,'/',fname,'.dat'],radius,i,0);
    [pole_location,~,gof(I)] = poletracker(fname,radius,i,0);
  

    barPos(I,:) = pole_location;
    
    %% Load whisker tracking solution and compute distance to whisker
    % V1 uses r_base
    temp = WhikerMan_functions({squeeze(1+r_base(i,:,:)),[0:0.01:ext_range]},8); b_tmp = temp{:};
    
    % Compute distance from all whisker points to pole center, to determine
    % closest point
    whisk_dist = sqrt((b_tmp(1,:)-pole_location(1)).^2 + (b_tmp(2,:)-pole_location(2)).^2);
    
    [pole_dist(I),closest_point] = min(whisk_dist);
    
    % Save closest point on whisker to pole
    closest_w(I,:) = b_tmp(:,closest_point)';
    
end


% Sort out indexing after parfor slices
gof(frames) = gof(1:numel(frames));
barPos(frames,:) = barPos(1:numel(frames),:);
pole_dist(frames) = pole_dist(1:numel(frames));
closest_w(frames,:) = closest_w(1:numel(frames),:);

gof(start_frame:trigger_frame) = 0;
barPos(start_frame:trigger_frame,:) = 0;
pole_dist(start_frame:trigger_frame) = 0;
closest_w(start_frame:trigger_frame,:) = 0;

tot_time = toc(time1);
ind_time = tot_time./numel(frames);

disp(['Total tracking time = ',num2str(tot_time),'s for ',num2str(numel(frames)),' frames. ',num2str(ind_time),'s per frame']);

%% Save variables for later contact estimation

% Append only works if the file already exists
if exist(fullfile([fname(1:end-4),'_touch.mat']), 'file') % local_dir,'/',
%     save([local_dir,'/',fname,'_touch.mat'],'start_frame','trigger_frame','pole_dist','gof','radius','barPos','closest_w','-append');
    save([fname(1:end-4),'_touch.mat'],'start_frame','trigger_frame','pole_dist','gof','radius','barPos','closest_w','-append');

else
    disp(['Creating file ',fname(1:end-4),'_touch.mat'])
%     save([local_dir,'/',fname,'_touch.mat'],'start_frame','trigger_frame','pole_dist','gof','radius','barPos','closest_w');
    save([fname(1:end-4),'_touch.mat'],'start_frame','trigger_frame','pole_dist','gof','radius','barPos','closest_w');

end

%% Delete local copy of .dat files

% Uncomment/ edit the next section when saving to a local folder
disp(['Removing local copy of ',fname])
% UNIX
copyfile([fname(1:end-4),'_touch.mat'],[remote_dir,'/',old_fname(1:end-4),'_touch.mat']);
copyfile([fname(1:end-4),'_clean.mat'],[remote_dir,'/',old_fname(1:end-4),'_clean.mat']);
% WINDOWS
% copyfile([fname(1:end-4),'_touch.mat'],[remote_dir,'/',old_fname(1:end-4),'_touch.mat']);
% copyfile([fname(1:end-4),'_clean.mat'],[remote_dir,'/',old_fname(1:end-4),'_clean.mat']);

delete([fname]);
delete([fname(1:end-4),'_touch.mat']);
delete([fname(1:end-4),'_clean.mat']);

disp(['Total time including file transfers =',num2str(toc(time0)),'s'])
