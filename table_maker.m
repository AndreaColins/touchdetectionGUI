% table_maker.m
%
% A script for making tables of data and pretty graphs
% MORE NOTES TO BE ADDED IN TIME
%

%% First, load and plot some data
load('170415_32a_20150417_133039_clean.mat')
load('170415_32a_20150417_133039_touch.mat')
% figure % new figure
clf % or just clear the exisiting figure
plot(touches)
hold all % Hold the current axes in place for more plots
plot(zscore(kappa_w))
plot(zscore(theta_w))

legend('touches','zscore kappa','theta')

%% Code to plot whether the mouse was protracting or retracting

% cont_angle_Vec = barPos - closest_w;
% cont_angle = atan2(cont_angle_Vec(:,1),cont_angle_Vec(:,2))*180./pi;

% Protraction vs retraction touch
theta_ts = timeseries(theta_w,(1:numel(theta_w))./1000);
bandpass = [6,30];
theta_filt = idealfilter(theta_ts,bandpass,'pass');
H = hilbert(theta_filt.data);

pro = find(angle(squeeze(H))<=0);
ret = find(angle(squeeze(H))>=0);

%% Plot touches on top of angle in protraction/retraction
clf
plot(theta_w)
hold all
plot(pro,theta_w(pro),'g.')
plot(ret,theta_w(ret),'m.')
plot(find(touches),theta_w(find(touches)),'ko')

%% Plot touches on top of curvature in protraction/retraction
clf
plot(kappa_w)
hold all
plot(pro,kappa_w(pro),'g.')
plot(ret,kappa_w(ret),'m.')
plot(find(touches),kappa_w(find(touches)),'ko')
%% Find touch times
touch_times = find(touches);
first_touch = find((touch_times > start_frame),1,'first');

protraction_touch = ismember(first_touch,pro);

%% Now a loop to load all the data (that has been tracked) and computing whether first touch was during protraction or retraction

% First compute filenames from touch_params
clear all;

touchcsv = 'touch_params.csv'; % Name of spreadsheet with file names etc
touch_params = csvread(touchcsv); % Load data from spreadsheet;

clear tracked_files filenames trialtypes choices % clear variables before use

tracked_files = find(touch_params(:,4));
filenames = touch_params(tracked_files,2);
trialtypes = touch_params(tracked_files,3);
choices = touch_params(tracked_files,14);

results.trial = filenames;
results.trialtype = trialtypes;
results.choice = choices;

touches_found = zeros(size(trialtypes));
protraction_touch = zeros(size(trialtypes));

phase_touch = zeros(size(trialtypes));
angle_touch = zeros(size(trialtypes));

for i = 1: numel(filenames)
    clear theta_w start_frame touches
    cleanfile = ['170415_32a_20150417_',num2str(filenames(i)),'_clean.mat'];
    touchfile = ['170415_32a_20150417_',num2str(filenames(i)),'_touch.mat'];
    
    load(cleanfile,'theta_w');
    load(touchfile,'start_frame','touches');
    
    theta = circshift(theta_w',[-start_frame,0]); % circularly permute theta to start at the start_frame
    touches = circshift(touches',[-start_frame,0]); % same for the touches array
    
    
    
    theta_ts = timeseries(theta,(1:numel(theta_w))./1000);
    bandpass = [6,30];
    theta_filt = idealfilter(theta_ts,bandpass,'pass');
    H = hilbert(theta_filt.data);
    phase = angle(squeeze(H)); % whisker phase
    pro = find(phase<=0);
    
    if numel(find(touches)) >= 1
        touches_found(i) = 1;  % Make a note whether there are any touches in this trial
        first_touch = find(touches,1,'first');
    
        ismember(first_touch,pro)
        protraction_touch(i) = 1; % Note whether first touch was protraction or not
        
        phase_touch(i) = phase(first_touch); % Read off phase at first touch
        angle_touch(i) = theta(first_touch); % Read off theta at first touch
    else
       touches_found(i) = 0; 
    end
end


results.touches_found = touches_found;
results.protract = protraction_touch;
results.phase_touch = phase_touch;
results.angle_touch = angle_touch;

% save results 
save results results 

%% try plotting some results
% you can load the results array with 'load results.mat' 

plot(results.trialtype,results.protract,'.')
xlim([0,4]);
ylim([-0.5,1.5]);

%% A better plot with averages
tt1 = find(results.trialtype == 1);
tt2 = find(results.trialtype == 2);
tt3 = find(results.trialtype == 3);

avg_pro1 = mean(results.protract(tt1));
avg_pro2 = mean(results.protract(tt2));
avg_pro3 = mean(results.protract(tt3));

plot([1,2,3],[avg_pro1,avg_pro2,avg_pro3],'.');
xlim([0,4]);
ylim([-0.5,1.5])

%% Homework for Chloe: 
% Try working out how to use the 'errorbar' plot
sem_pro1 = std(results.protract(tt1)) ./sqrt(numel(tt1));
sem_pro2 = std(results.protract(tt2)) ./sqrt(numel(tt2));
sem_pro3 = std(results.protract(tt3)) ./sqrt(numel(tt3));

errorbar([1,2,3],[avg_pro1,avg_pro2,avg_pro3],[sem_pro1,sem_pro2,sem_pro3])


ylim([-0.5,1.5]);

%% Plot angle at first touch for all trials vs trial type
plot(results.trialtype,results.angle_touch,'.')
xlim([0,4])
