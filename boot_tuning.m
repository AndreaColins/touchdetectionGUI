function boot_tuning(spks,stim,method,numbins,bootreps,colour)
%
% boot_tuning(spks,stim,numbins,bootreps)
%
% Function to produce bootstrapped tuning curves from an array of spikes
% and a stimulus array.
% Expects a 1D binary spikes array, a 1D stimulus array, and optional
% tuning curve bin value (default 20), bootstrap repetitions value (default
% 1000), and a 3 value line colour specification (default blue)
% 
% First verison uses uniform bins along the stimulus axis. An optional 
% method of 'nonuniform' uses non-uniform bins, with equal samples in each
% bin instead of equal stimulus range
%
% TO DO add non-uniform bins

if length(numbins) == 0;
    numbins = 20;
end

if length(bootreps) == 0;
    bootreps = 1000;
end

if length(colour) == 0;
    colour = [0,0,1];
end

% Switch for uniform or nonuniform bins
switch method
    case 'nonuniform'
        
        stim_z = stim;
        
        bin_edges_ind = round(linspace(1,numel(stim_z),numbins+1));
        
        sort_stim = sort(stim_z);
        
        bin_edges = sort_stim(bin_edges_ind)';
        
    otherwise
        
        % z-score stimulus
        stim_z = stim./std(stim);
        
        % Find stim range for bin edges
        min_z = min(stim_z);
        max_z = max(stim_z);
        
        % Define linear range from min to max
        bin_edges = linspace(min_z,max_z,numbins+1);
        
end

clear bootmean_spikes boot_ind
tic

bootmean_spikes = zeros(numbins,bootreps);
disp(['Bootstrapping tuning curves, (',num2str(bootreps),' samples)'])
% Bootstrapped mean firing rate in each bin
for j = 1:bootreps;
    boot_ind = randi(numel(spks),1,numel(spks)); % Random integers across the full range of spks
    %    boot_ind = sort(boot_ind); % Sorting for easier visualisation/debugging
    for i = 1:numel(bin_edges)-1
        ind = find(stim_z(boot_ind)>=bin_edges(i)& stim_z(boot_ind)<bin_edges(i+1));
        bootmean_spikes(i,j) = 1000.*nanmean(spks(boot_ind(ind)));
    end
end
toc

% Shift the spikes bootreps times and re-compute tuning curve
shuffmean_spikes = zeros(numbins,bootreps);

disp(['Shuffling tuning curves, (',num2str(bootreps),' shuffles)'])
rand_i = 1000 + round(4000.*rand(1,bootreps)); % bootreps random index shifts
for j = 1:bootreps;
    %             shuff_spks = spks_non_touch(randperm(numel(spks_non_touch)));
    index = 1:numel(spks);
    sh_index = [index(rand_i(j):end),index(1:rand_i(j)-1)];
    shuff_spks = spks(sh_index);
    
    for i = 1:numel(bin_edges)-1
        ind = find(stim_z>=bin_edges(i)& stim_z<bin_edges(i+1));
        shuffmean_spikes(i,j) = 1000.*nanmean(shuff_spks(ind));
    end
end
toc

% Work out median stim value per bin to set x axis for plotting
for i = 1:numel(bin_edges)-1
    x_axis(i) = median(stim_z(find(stim_z>=bin_edges(i)& stim_z<bin_edges(i+1))));
end

% % Mean firing rate in each bin for whole dataset
% for i = 1:numel(bin_edges)-1
%     ind = find(stim_z>=bin_edges(i)& stim_z<bin_edges(i+1));
%     numel(ind);
%     numel(find(spks(ind)));
%     all_mean(i) = 1000.*nanmean(spks(ind));
% end

% Set x axis ticks to midpoint of each bin
% x_axis = mean([bin_edges(1:numbins);bin_edges(2:numbins+1)]);

% Plot bootstraped tuning curves
% Colour should be a lighter version of main plot colour
samp_colour = colour + 0.8;
samp_colour(samp_colour>1) = 1;

% Plot shuffled tuning curves
plot(x_axis,shuffmean_spikes,'color',[.8,.8,.8]);
hold all;
% Compute 95th percentile of shuffled tuning curves
y = prctile(shuffmean_spikes',[2.5,97.5]);
plot(x_axis,y','k','Linewidth',1)

% Plot bootstrapped tuning curves
plot(x_axis,bootmean_spikes,'color',samp_colour);


% % Compute 95th percentile and median
% y = prctile(bootmean_spikes',[2.5,50,97.5]);

% Plot mean firing rate and std of bootstrapped samples
% plot(x_axis,nanmean(bootmean_spikes,2),'color',colour,Linewidth',2)
myeb(x_axis,nanmean(bootmean_spikes,2),0.5*nanstd(bootmean_spikes,0,2),colour,samp_colour)
hold all

% % Plot mean firing rate for whole dataset
% plot(x_axis,all_mean,'Linewidth',2,'color',colour)

% % Plot median firing rate and CI of bootstrapped samples
% plot(x_axis,y','k') 

% Plot location of bin edges
% stem(bin_edges,-0.2.*(max(nanmean(bootmean_spikes,2))).*ones(1,numel(bin_edges)),'k','marker','none');
stem(x_axis,0.2.*(max(nanmean(bootmean_spikes,2))).*ones(1,numel(x_axis)),'k','marker','none');
