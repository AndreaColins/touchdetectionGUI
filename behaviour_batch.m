% behaviour_batch.m
%
% Minimal script for converting files to .avi then running the janelia whisker tracker
%
% Use chimera.m to clean up the data, then extract kappa/angle variables using a mask, quadratic fit and whikerman functions.
%
% Also has instructions for running the tracker at the command line
%
% First videos need to be converted into .avi format with dat2mp4.m
% [~] = dat2mp4(FILENAME,FILENAME);
% Or cd into the video directory and run batch_dat2mp4.m (editted to
% include the appropriate session names.
%
% Second, run (in turn):
% trace FILENAME.avi FILENAME.whiskers
% measure --face bottom FILENAME.whiskers FILENAME.measurements
% classify FILENAME.measurements FILENAME.measurements bottom --px2mm 0.057 -n 1
% reclassify FILENAME.measurements FILENAME.measurements -n 0

%% List of directories to convert/track
% PC
% main_dir = 'P:\Dario Campagner\BEHAVIORAL MOVIES' 
% cd(main_dir)


% UBUNTU
cd /run/user/1000/gvfs/smb-share:server=130.88.94.172',share=test'/Dario/Behavioral_movies/;

% Animal 34 (slashes flipped for Unix
dirs_34 = dir(['34/']);
dirs_34 = dirs_34(3:end);
for d = 1:numel(dirs_34);
    subdir = dir(['34/',dirs_34(d).name,'/']);
    if numel(subdir) == 3;
        subdir = subdir(3).name;
        fulldir{d} = ['34/',dirs_34(d).name,'/',subdir,'/'];
    end
end
num_sofar = numel(fulldir);


% Single whisker sessions only
% 38
% keepers = [27,32,34,36,39,40,47,50,52,56,61,64,2,4,6,10,13,17,18,20,33,35,...
%     37,44,53,54,65,7,8,11,19];

% 34
keepers = [19,21,22,25,2,13,16,23,26,29,33,38,52,54,56,3,6,7,9,10,14,17,20,24,27,34,30,41,44,47,49,50,58];
% 21,23,24,27,4,15,18,25,28,31,35,40,54,56,58,5,8,9,11,12,16,19,22,26,29,36,32,43,46,49,51,52,60

dirs_1_whisker = {};
for d = 1:numel(keepers);
    dirs_1_whisker{d} = fulldir{keepers(d)};
end

% BAD FILES 38: 4|140. 5|5,58,62,84. 7|145.
% 8|52,105.15|1,54.16|1,59,75.17|94(crash).22-end(no janelia).
% 22 no good_trials.23 chimera.
% chimera done: 24.25.26.27.

%% Main conversion loop
for d =  1:numel(dirs_1_whisker); % For all chosen sessions (above)

    %% CD to a video directory, then convert all .dat files to .avi

%     cd(main_dir) % PC
    cd /run/user/1000/gvfs/smb-share:server=130.88.94.172',share=test'/Dario/Behavioral_movies/; % UBUNTU
    dir_name = dirs_1_whisker{d};
    cd(dir_name)

    files = dir('*.dat');
    names = files(1);
    names = names.name;
    session = names(1:10);
    vid = names(12:19);
    
    

    fidx = 0;
    %% Main loop. dat - avi - janelia - chimera - contact_detector
    % Check if .dat has been converted to .avi first
    for fidx = fidx+1:length(files)
        loadfile = files(fidx).name(1:end-4);

%         local_dir = 'D:\evans_work\'; % PC
        local_dir = '~/scratch/'; % UBUTNU
        %         remote_dir = pwd;
        savefile = [local_dir,files(fidx).name(1:end-4),'.avi'];

        % Has this file already been converted to .avi?
        if ~exist([files(fidx).name(1:end-4),'.avi'])
            savefile_nas = [files(fidx).name(1:end-4) '.avi'];


            disp(sprintf('Converting %s to .avi ...',loadfile))
            [~] = dat2mp4(loadfile,savefile);
            disp(sprintf('Saved %s',loadfile));

            % Copy files to NAS. Delete local file. Need clever way of doing
            % this after local file is no longer needed
            disp(sprintf('Copying %s to NAS ...',loadfile))
            copyfile(savefile,savefile_nas);
            %         delete(savefile);
        end
        %    end
        %     clear fidx files loadfile savefile

        %% Check whether there is a local .avi file
        if ~exist(savefile)
            savefile_nas = [files(fidx).name(1:end-4) '.avi'];
            copyfile(savefile_nas,savefile);
        end

        %% Now run Janelia tracker on the videos

        %     files = dir([session '_' vid '_*.avi']);
        if ~exist([files(fidx).name(1:end-4),'.whiskers'])
            disp(['Janelia tracking video ',files(fidx).name(1:end-4)])

            % for fidx = 1:length(files)
            %         loadfile = files(fidx).name(1:end-4);
            %         savefile = [loadfile,'.avi'];

            disp(['Tracing vid ',loadfile])
            str = ['trace ',savefile,' ',loadfile,'.whiskers'];
            system(str);
        end
         
        if ~exist([files(fidx).name(1:end-4),'.measurements'])
            disp(['Measuring vid ',loadfile])
            str = ['measure --face bottom ',loadfile,'.whiskers ',loadfile,'.measurements'];
            system(str);
            
            disp(['Classifying vid ',loadfile])
            str = ['classify ',loadfile,'.measurements ',loadfile,'.measurements bottom --px2mm 0.057 -n 1'];
            system(str);
            
            disp(['Re-classifying vid ',loadfile])
            str = ['reclassify ',loadfile,'.measurements ',loadfile,'.measurements -n 1'];
            system(str);
        end

        % Delete local copy of video
        disp(sprintf('Deleting local copy of %s ',loadfile))
        delete(savefile);

        % Then run chimera.m (currently doing this separately, eventually it
        % would happen here, such that the pipeline is run end-to-end on a
        % given session using only one script

        %% Chimera
        % run from current directory, providing a loadfile
        %   [best_m,best_w,kappa_w,theta_w,r_base,dropped_frames] = chimera(loadfile);
        %  save([loadfile,'_clean.mat'],'kappa_w','theta_w','r_base','dropped_frames','-append');

        % TO DO: add check for .whiskers file

        % Run chimera in script mode to see progress etc
        plotting = 0;
%         chimera
    end
end
        % Copy files to NAS

        % janelia files
        % copyfile();
        % chimera tracker

        % TO DO: add check for _clean.mat
        %% Contact detector
        radius = 13; % Trying with 13 on 0408 to see if if helps in gof peak definition
        plotting = 0;
        %     files = dir('*_clean.mat');

        %     for j = 1:numel(files);
        %         fname = files(j).name;
        %         fname = [fname(1:end-10),'.dat'];
        %         %  fname = fname(1:end-4);
        %
        fname = [loadfile,'.dat'];
        contact_detector_parallel(fname,radius); % NOTE: creates local version of files, so change this option depending on computer being used
        %     contact_detector(fname,radius,plotting);
        %     end
        %     j = 1;
    end
end
