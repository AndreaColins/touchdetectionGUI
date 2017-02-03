% Script to convert a given .dat video file into 5 different compressed
% versions, then run the janelia tracker on each.
%
% This is to assess the effects of video compression on whisker tracking
% quality. All videos will subsequently be tracked with whikerman (ideally
% by more than one person).
%
% The 4 formats are:
% 'Archival': lossless compression in .mj2 format
% 'Motion JPEG AVI', Quality 100
% 'Motion JPEG AVI', Quality 95
% 'Motion JPEG AVI', Quality 85
% 'Motion JPEG AVI', Quality 75 (default)

% State which file you want to convert/track
loadfile = '270415_32a_20150427_150040'; % e.g. in C:\Users\mathew\work\Data\Test_videos\avi_test

% Conversion/saving options to loop over
prof = {'Archival','Motion JPEG AVI','Motion JPEG AVI','Motion JPEG AVI','Motion JPEG AVI'};
suffix = {'Archival','JPEG_Qual_100','JPEG_Qual_95','JPEG_Qual_85','JPEG_Qual_75'};
qual = [100,100,95,85,75];


% Loop to convert and track the files
for i = 2:5;
    
    fid = fopen([loadfile '.dat'],'r');
    
    % Load in header information
    data = read_mikrotron_datfile_header(fid);
    
    data.loadfile = loadfile;
    data.savefile = [loadfile,'_',suffix{i},'.avi'];
    
    % data.profile = 'MPEG-4'; % This is tiny. Good for talks. Cannot be read by Janelia Whisker Tracker
    % data.profile = ; % This is huge. DO NOT USE
    
    data.profile = prof{i}; 'Motion JPEG AVI'; % This is the default
    
    
    %
    % Create mov file to write individual images (img)
    fseek(fid,data.offset,-1);
    
    mp4obj = VideoWriter(data.savefile,data.profile);
    
    mp4obj.FrameRate = data.framerate/20;
    
    if i>=2;
        mp4obj.Quality = qual(i);
    end
    
    
    open(mp4obj)
    
    count = 0;
    stop = 0;
    fprintf('Processing %d frames:\n',data.nframes)
    while(~stop)
        if rem(count,100)==0
            fprintf('%d ', count)
        end
        count = count+1;
        imagecounter = fread(fid,2,'uint8');
        junk = fread(fid,2,'uint8');
        tickcounter = fread(fid,2,'uint32');
        digitalinputs = fread(fid,1,'uint8');
        junk = fread(fid,11,'uint8');
        frame = fread(fid,data.imagesize-24,'uint8=>uint8');
        stop = feof(fid);
        if ~stop
            digitalbit0(count) = bitand(digitalinputs,1);
            frame_res = reshape([zeros(24,1); frame],data.width,data.height)';
            img = repmat(frame_res, [1,1,3]);
            writeVideo(mp4obj,img);
        end
        
    end
    
    fprintf('done\n')
    close(mp4obj);
    
    %% Janelia tracker stuff
        this_loadfile = data.savefile(1:end-4);
      
        disp(['Tracing vid ',this_loadfile])
        str = ['trace ',data.savefile,' ',this_loadfile,'.whiskers'];
        system(str)
        
        disp(['Measuring vid ',this_loadfile])
        str = ['measure --face bottom ',this_loadfile,'.whiskers ',this_loadfile,'.measurements'];
        system(str);
        
        disp(['Classifying vid ',this_loadfile])
        str = ['classify ',this_loadfile,'.measurements ',this_loadfile,'.measurements bottom --px2mm 0.057 -n 1'];
        system(str);
        
        disp(['Re-classifying vid ',this_loadfile])
        str = ['reclassify ',this_loadfile,'.measurements ',this_loadfile,'.measurements -n 1'];
        system(str);
         
end

%% Plot all 5 solutions.
    figure(2);
    clf;
for i = 1:5;
    whiskerfile = [loadfile,'_',suffix{i},'.whiskers'];
    measurementfile = [loadfile,'_',suffix{i},'.measurements'];
    
    w = LoadWhiskers(whiskerfile);
    m = LoadMeasurements(measurementfile);
    

    hold all
    
    dt = 0.001;
    time = double([m(:).fid]).*dt;
    angle = [m(:).angle];
    curvature = [m(:).curvature];
    %     colour = varycolor(3);
    
    whisker_ID = 0;%:max([m(:).label]);
    mask = [m(:).label] == whisker_ID;
    theta = angle(mask);
    kappa = curvature(mask);
    theta = theta - mean(theta);
    kappa = kappa - mean(kappa);
    theta(abs(theta)>5*std(theta)) = 0;
    kappa(abs(kappa)>2*std(kappa)) = 0;
    
    ax(1) = subplot(2,1,1);
    plot(time(mask),theta);
    hold all
    ylabel('Angle');
    legend(suffix,'Location','EastOutside')

    ax(2) = subplot(2,1,2);
    plot(time(mask),kappa);
    ylabel('Kappa');
    legend(suffix,'Location','EastOutside')
    hold all
    
    xlabel('Time (s)');
    drawnow;
    
end


linkaxes(ax,'x');