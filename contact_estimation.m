% contact_estimation.m
%
% Script to estimate whisker-pole contact times from tracked whiskers and pole
% location information. Either needs to be run after contact_detector.m has
% finished on a set of videos, or contact_detector needs to be run here.


%% Set inital parameters
path = '/run/user/1000/gvfs/smb-share:server=nasr.man.ac.uk,share=flsrss$/snapped/replicated/Petersen/Dario Campagner/BEHAVIORAL MOVIES/36/050815_36a/2015_08_05';
path = '/media/mathew/Bigger Data/whisker_movies/050815_36a/2015_08_05'
cd (path);


radius = 14;
plotting = 1;
files = dir('*.dat');

dist_dist = [];
ang_dist = [];
touch_ang = [];
pole_pos = [];

%% Optionally run contact_detector/parallel on all videos in this folder
radius = 14; % Trying with 13 on 0408 to see if if helps in gof peak definition
plotting = 0;
files = dir('*.dat');

for i = 191:numel(files); 
    fname = files(i).name;
  %  fname = fname(1:end-4);
    
    contact_detector_parallel(fname,radius);
%     contact_detector(fname,radius,plotting);
end


%% Batch version (across days) of above code

% Dates from threeposition.m
% dates =  [3,7,29,32,35,38,41,45,47,49,50,56]; % 32. Done: 1,2,3,4 (chloe),5,6,7,8,9,10,11,
dates = [30,34,23,24,27,25,28,31,35,54,56]; % 34. Done: 1 (chloe)
% dates = [60,64,67,68,12,14,17,20,21,24,27,63,66,4,6,9,16]; % 36. Done 1,2,3,4,5,6,7,9,10 13 (chloe). Not done - 8 (all ON)
for i = 2:11; %[11:numel(dates)]
    i
    cd /run/user/1000/gvfs/smb-share:server=130.88.94.172',share=test'/Dario/Behavioral_movies/34/ % NAS. 32, 34
%     cd /run/user/1000/gvfs/smb-share:server=nasr.man.ac.uk',share=flsrss$/snapped/replicated/Petersen'/Dario' Campagner'/BEHAVIORAL' MOVIES'/36/ % Isilon. 36
    files = dir;
    cd(files(dates(i)).name);
    files(dates(i)).name
    x = cellstr(ls);
    cd (x{:})
    
    radius = 14; 
    plotting = 0;
    files = dir('*.dat');
   %% 
    for j = 1 : numel(files);
        j
        fname = files(j).name;
        % fname = fname(1:end-4);
        
        %     contact_detector_parallel(fname,radius);
             contact_detector(fname,radius,plotting);
    end
    
    
end

%% Load relevant tracking outputs
local_dir = '/media/mathew/Bigger Data1/janelia_tracker_scratch';

xlsfile = 'good_trials.xlsx'; % Must be saved as .xls Excel 98 format for compatability
xls_info = xlsread(xlsfile);
% Read out pole position (trial type)
polePos = xls_info(1:90,9);

for i = 1:numel(files);
    fname = files(i).name;
    fname = fname(1:end-4);
    
    % Load pole information
    load([local_dir,'/',fname,'_touch.mat'],'start_frame','trigger_frame','pole_dist','gof','radius','barPos','closest_w');
    
    % Load tracked whisker information. Copy file over if it doesn't exist
    % locally
    if ~exist(fullfile([local_dir,'/',fname,'_clean.mat']), 'file')
        disp(['Copying ',fname,'_clean.mat locally']);
        copyfile([fname,'_clean.mat'],[local_dir,'/',fname,'_clean.mat']);
    end
    
    load([local_dir,'/',fname,'_clean.mat'],'r_base','theta_w','kappa_w','best_m');
    
    %% Determine pole up times for a video
    pole_up_window = trigger_frame+200:trigger_frame + 700;
    pole_up_window = mod(pole_up_window,numel(gof)) + 1;
    
    % Set gof threshold to mean(gof)-4.*std(gof)
    gof_peak_thresh = nanmean(gof(pole_up_window))+3.*nanstd(gof(pole_up_window));
    
    % Find 2 peaks in gof corresponding to pole-in-focal-plane times
    [~,pole_updown] = findpeaks(gof,'MinPeakHeight',gof_peak_thresh,'MinPeakDistance',1000,'Npeaks',2);
    
    if numel(pole_updown) == 1;
        disp(['file ',fname,' only has one gof peak']);    
        pole_updown(2) = start_frame;
    elseif numel(pole_updown) > 2;
        disp(['file ',fname,' has more than 2 gof peaks']);
        pole_updown = pole_updown(1:2);
    elseif numel(pole_updown) == 0;
        disp(['No gof peaks found in file ',fname]);
        pole_updown = [trigger_frame,start_frame];
    end
    
    % Work out which order the peaks are (depends on timing of trigger in
    % the video
    if numel(find(pole_updown>=trigger_frame)) == 2||0;
        pole_up_times = pole_updown(1):pole_updown(2);
    else
        pole_up_times = [pole_updown(2):numel(pole_dist),1:pole_updown(1)];
    end
    
    % Work out distance to the pole for all frames during pole-up
    candidate_touch = pole_dist(pole_up_times);
    
    % Whisker_angle during pole up.
    angle_touch = theta_w(pole_up_times);
    
    % Work out angles of contact wrt base and contact wrt pole center 
    % for each frame during pole up, to help  with contact classification
    basePos = r_base(pole_up_times,:,1);
    
    % Angle of touch from whisker base - same frame of reference as theta_w
    touch_angle_Vec = basePos - closest_w(pole_up_times,:);
    touch_angle = atan2(touch_angle_Vec(:,1),touch_angle_Vec(:,2)) *180./pi;
    
    % Angle of touch from pole centre. > 180 = protraction touch
    cont_angle_Vec = barPos(pole_up_times,:) - closest_w(pole_up_times,:);
    cont_angle = atan2(cont_angle_Vec(:,1),cont_angle_Vec(:,2))*180./pi;

    %% Aggregate pole_dist and cont_angle during pole_up_times for all videos in a session to
    %  determine data-driven distance thresholds for touch
    dist_dist = [dist_dist,candidate_touch];
    
    ang_dist = [ang_dist;cont_angle];
    
    touch_ang = [touch_ang; touch_angle];
    
    pole_pos = [pole_pos; polePos(i)*ones(numel(pole_up_times),1)];

    
end

%% Plot ang_dist and dist_dist, coloured by pole position
% load('/Volumes/05/Video_data/post_processing.mat') % Mac

plot(ang_dist(find(pole_pos==5e4)),dist_dist(find(pole_pos==5e4)),'.')
hold all
plot(ang_dist(find(pole_pos==105e3)),dist_dist(find(pole_pos==105e3)),'.')
plot(ang_dist(find(pole_pos==155e3)),dist_dist(find(pole_pos==155e3)),'.')

%% Divide data into 6 classes based on ang_dist sign and pole_pos
pp_1 = find(pole_pos==5e4);
pro_1 = find(ang_dist(pp_1)>=0);
ret_1 = find(ang_dist(pp_1)<0);
touch_class{1} = pp_1(pro_1);
touch_class{2} = pp_1(ret_1);

pp_2 = find(pole_pos==105e3);
pro_2 = find(ang_dist(pp_2)>=0);
ret_2 = find(ang_dist(pp_2)<0);
touch_class{3} = pp_2(pro_2);
touch_class{4} = pp_2(ret_2);

pp_3 = find(pole_pos==155e3);
pro_3 = find(ang_dist(pp_3)>=0);
ret_3 = find(ang_dist(pp_3)<0);
touch_class{5} = pp_3(pro_3);
touch_class{6} = pp_3(ret_3);


%% Plot all 6 clusters in 3D
clf
for i = 1:6;
    touch = touch_class{i};
    plot3(ang_dist(touch),touch_ang(touch),dist_dist(touch),'.');
    drawnow;
    hold all;
end

%% Histogram of pole distance
clf
for i = 1:6;
    touch = touch_class{i};
    subplot(3,2,i);
    hist(dist_dist(touch),50);
end

%% Plot of pole distance
clf
for i = 1:6;
    touch = touch_class{i};
    subplot(3,2,i);
    plot(dist_dist(touch));
end

%% TO DO. Plot pole distance alongside kappa/theta for trials with definite contacts. Plot alongside movies
% Maybe just skip straight to touch viewer gui code?