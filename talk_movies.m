% talk_movies.m
% A script to generate movies for talks from data

% CD to directory with data
% cd '/media/mathew/Bigger Data/whisker_movies/Campagner_eLife_movie'
% cd '/run/user/1000/gvfs/smb-share:server=130.88.216.49,share=shared/'

% File names for loading/saving 
data.loadfile = files(6).name(1:end-4); %'210514b_20140521_144602';
data.savefile = 'calibrate_movie_2';%'campagner_SI_mp4';
% data.savefile = 'campagner_SI_smearSpks';

% load kappa, theta and tracking solution
% load([data.loadfile,'.tr'],'kappa_all','theta_all','rall','-mat');

% Some plotting params
frames = 1422:4432;
blue = [0 120 255];
orange = [255 120 0];

% % Mikrotron setup
video.type = 'dat';
video.fid = fopen([data.loadfile,'.dat'],'r');
video.header = read_mikrotron_datfile_header(video.fid);
handles.nframes = video.header.nframes;
video.width = video.header.width;
video.height = video.header.height;
video.offset = 8192;
fseek(video.fid,8192,-1);

handles.video = video;

% New movie setup
data.profile = 'Motion JPEG AVI'; % This is the default. Update mp4.Quality below to ensure best conversion quality possible (note setting quality to 100 increases files size from 37MB to 187MB)
% data.profile = 'MPEG-4'; % This is tiny. Good for talks. Cannot be read by Janelia Whisker Tracker, and doesn't work on Linux
mp4obj = VideoWriter(data.savefile,data.profile);
mp4obj.FrameRate = 50;%video.header.framerate/20;
mp4obj.Quality = 100; % Default is 75

open(mp4obj)

% Load spikes
% load('data_with_fixed_contact_neuron5_6_8_10_9_16_04102015','SPIKES');
% spikes = SPIKES{4};
% spikes = [zeros(15,200),spikes];
for i = 1:numel(frames);% 2000:3000;
    frameidx = frames(i);
    offset = video.header.imagesize * (frameidx-1) + video.offset;
    fseek(video.fid,offset,-1);
    tmp = fread(video.fid,video.header.imagesize-24,'uint8=>uint8');
    tmp = reshape([tmp; zeros(24,1)],video.width,video.height)';
    frame = uint8(zeros(video.height,video.width,3));
    frame(:,:,1) = tmp;
    frame(:,:,2) = tmp;
    frame(:,:,3) = tmp;
    clear tmp
    
%     [C,img] = texture_calibrate(files(6).name,i,1,0);
    clf; hold all;
    image(frame);
    set(gca,'Ydir','reverse');
    % Locate bezier curve from rall
    temp = WhikerMan_functions({squeeze(rall(frameidx,:,:)),0:.01:1},8);
    bez = temp{:};
    plot(bez(1,:),bez(2,:),'r','Linewidth',3);
    plot(296- 100.*kappa_all(frameidx-200:frameidx)./max(kappa_all(frames)),'color',blue./255,'Linewidth',3)
    hold all
    plot(20 + 296 - 200.*theta_all(frameidx-200:frameidx)./max(theta_all(frames)),'color',orange./255,'Linewidth',3)
    plot(find(spikes(7,i:i+200)),100.*ones(1,numel(find(spikes(7,i:i+200)))),'w.','Markersize',10)
    axis off;axis image;
    
    % Add axis bar @ k=0.05 mm-1, theta = 40 deg
    plot([390,390],[280-60,280],'color',blue./255,'Linewidth',3)
    plot([390,390],[210-57,210],'color',orange./255,'Linewidth',3)
    
    % Add text 
    text(317,180,'\bf{Angle}','color',orange./255,'Fontsize',14)
    text(300,245,'\bf{Curvature}','color',blue./255,'Fontsize',14)
    text(310,260,'\bf{change}','color',blue./255,'Fontsize',14)


    
%     if spikes(7,i+200);
%         plot(25,25,'w.','Markersize',50);
%     end
    drawnow;

    img = getframe(gcf);
    writeVideo(mp4obj,img);
    
end

fclose('all') ; % Trying to fix 'too many open files' error

fprintf('done\n')
close(mp4obj)


%% Audio from spikes
spks = upsample(spikes(7,1:3010),160);
spks(find(spks>1)) = 1;
s = conv(spks,[ones(1,2),-1*ones(1,2),ones(1,2),-1*ones(1,2)],'same');
wavwrite(s,8000,'spikes_8KHz');

% NOTE: ffmpeg was used to combine the two sources using the command:
% ffmpeg -i campagner_SI_smearSpks.avi -i spikes_8KHz.wav -codec copy -shortest campagner_SI_nodot.avi
% Then changed to mp4 with:
% ffmpeg -i campagner_SI.avi -vf scale=500:-1 -c:v libx264 -pix_fmt yuv420p -b:v 1000k campagner_SI.mp4

%% Kappa theta space
clf;
data.loadfile = '300915_38a_20150930_151506';
load([data.loadfile,'.tr'],'kappa_all','theta_all','rall','-mat');
plot(180-theta_all,kappa_all,'.'); 
hold all;
data.loadfile = '300915_38a_20150930_151630';
load([data.loadfile,'.tr'],'kappa_all','theta_all','rall','-mat');
plot(180-theta_all,kappa_all,'.');
data.loadfile = '300915_38a_20150930_151750';
load([data.loadfile,'.tr'],'kappa_all','theta_all','rall','-mat');
plot(180-theta_all,kappa_all,'.');
legend('Correct L','Correct R','Incorrect R')

