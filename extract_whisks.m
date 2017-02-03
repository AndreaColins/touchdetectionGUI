function [whisk_starts,varargout] = extract_whisks(theta,criterion)
% extract_whisks.m
% Script to isolate individual whisk cycle times from a long timeseries of
% whisker angle. Returns whisk cycle starts by default.
% Optionally returns whisk cycle ends and phi depending on varargin options
%
% Available varargin options: 
% 'all': returns all whisk cycles regardless of properties, also returns
% phi + amplitude to allow post-hoc filtering of whisks
% 'parsed': returns only whisks that pass certain criteria e.g. amplitude.
% Also returns whisk cycle end-times, phi and amplitude.

% M.Evans 8.12.15

switch(criterion)
    case {'all'}
        
        % Extract whisk amplitude using the hilbert transform
        clear ts H
        ts = timeseries(theta,(1:length(theta))./1000);
        bandpass = [6 50];
        theta_filt = idealfilter(ts,bandpass,'pass');
        % Hilbert transform to find amplitude and phase
        H = hilbert(theta_filt.data);
        amp = squeeze(abs(H));
        phi = angle(H);
        dphi  = diff(phi);
        whisk_starts = [1; 1+ find(dphi<-2)];
        
        % Extract all whisks
        whisk_ends = [whisk_starts(2:end); numel(theta_filt.data)];
        whisks = zeros(300,numel(whisk_starts));
        phis = zeros(300,numel(whisk_starts));
        amps = zeros(300,numel(whisk_starts));
        for i = 1:numel(whisk_starts); 
            whisks(1:1+(whisk_ends(i)-whisk_starts(i)),i) = theta_filt.data(whisk_starts(i):whisk_ends(i));
            phis(1:1+(whisk_ends(i)-whisk_starts(i)),i) = phi(whisk_starts(i):whisk_ends(i));
            amps(1:1+(whisk_ends(i)-whisk_starts(i)),i) = amp(whisk_starts(i):whisk_ends(i));
        end
        
        varargout = {whisk_ends,whisks,phis,amps};

        
        
    case {'parsed'}
        
    otherwise
        disp('Unknown criterion. Try all or parsed')
        
end