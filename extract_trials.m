function returned_trials = extract_trials(Data,unit,session,criterion)
% extract_trials.m
% Function to work out best trials for a given cell and session in a
% dataset, given some criterion, either opto trials, stable ephys or good behaviour. 
% Available options are:
% 'opto','good','stable','goodstable','goodopto','stableopto','goodstableopto'
% This function supersedes some copy-and-paste code I was using before
% NOTE 15.12.15: session 3 is a special case where the opto mirror val doesn't match
% the other sessions. This session has no opto trials so this is set
% manually.
% 
% M.Evans 07.09.15

switch(criterion)
    case{'good'}
    % Trials with good behaviour but no opto
    good_trials = Data(session).selTrials{unit};
    opto = find(Data(session).events.mirror.val<0.9);
    if session == 3;
        opto = [];
    end
    good_trials(ismember(good_trials,opto)) = [];
    
    returned_trials = good_trials;

    case{'stable'}
    % Trials with stable ephys, but no opto
        stable_trials = Data(session).stableTrials{unit};
        opto = find(Data(session).events.mirror.val<0.9);
        if session == 3;
            opto = [];
        end
        stable_trials(ismember(stable_trials,opto)) = [];
        
        returned_trials = stable_trials;
        
    case{'opto'}
    % All trials with opto
        opto = find(Data(session).events.mirror.val<0.9);

        if session == 3;
            opto = [];
        end

        returned_trials = opto;

        
    case{'goodstable'}
        % Trials with both good behaviour and stable ephys, but no opto
        good_trials = Data(session).selTrials{unit};
        opto = find(Data(session).events.mirror.val<0.9);
        if session == 3;
            opto = [];
        end
        good_trials(ismember(good_trials,opto)) = [];

        stable_trials = Data(session).stableTrials{unit};
        
        non_stable_trials = 1:size(Data(session).licks,1);
        non_stable_trials(stable_trials) = [];
        good_stable_trials = good_trials;
        good_stable_trials(ismember(good_stable_trials,non_stable_trials)) = [];
                
        returned_trials = good_stable_trials;
        
        
    case{'goodopto'}
        % Trials with both good behaviour and opto
        good_trials = Data(session).selTrials{unit};
        opto = find(Data(session).events.mirror.val<0.9);
        if session == 3;
            opto = [];
        end
        non_opto = 1:size(Data(session).licks,1);
        non_opto(opto) = [];
        good_opto_trials = good_trials;
        good_opto_trials(ismember(good_trials,non_opto)) = [];
    
        returned_trials = good_opto_trials;
    
        
    case{'goodstableopto'}
        % Trials with good behaviour, stable ephys and opto
        good_trials = Data(session).selTrials{unit};
        opto = find(Data(session).events.mirror.val<0.9);
        if session == 3;
            opto = [];
        end
        non_opto = 1:size(Data(session).licks,1);
        non_opto(opto) = [];
        good_trials(ismember(good_trials,non_opto)) = [];

        stable_trials = Data(session).stableTrials{unit};
        non_stable_trials = 1:size(Data(session).licks,1);
        non_stable_trials(stable_trials) = [];
        good_stable_trials = good_trials;
        good_stable_trials(ismember(good_stable_trials,non_stable_trials)) = [];
                
        returned_trials = good_stable_trials;
        
    case{'stableopto'}
        % All trials with stable ephys and opto
        stable_trials = Data(session).stableTrials{unit};
        opto = find(Data(session).events.mirror.val<0.9);
        if session == 3;
            opto = [];
        end
        non_opto = 1:size(Data(session).licks,1);
        non_opto(opto) = [];
        stable_opto_trials = stable_trials;
        stable_opto_trials(ismember(stable_trials,non_opto)) = [];
        
        returned_trials = stable_opto_trials;
        
        
    otherwise
        disp('Unknown criterion. Try good, stable, opto, goodstable, goodopto, goodstableopto, stableopto')


end


if numel(returned_trials) == 0;
    display(['No ',criterion,' trials for cell ',num2str(unit),' session ',num2str(session)])
end