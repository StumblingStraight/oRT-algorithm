%identifyEpochofInterest
%
%This script identifies redirections in the data using the Epoch of
%Interest method. Redirections are identified in acceleration traces by
%points where the acceleration exceeds a certain threshold for a certian
%duration, and results in pre-defined amount of change in X-Y heading.
%This script identifies between 1 and 2 redirections in each trial. This
%script also tracks the number of trials that have no redirections and the
%number that have too many redirections (3 or more). When redirections are
%identified, their approximate time is saved in the data structure "RTInd"
%field, and the RTs can be further refined with one of the "findORTs"
%scripts (such as findoRTs_CCS or findoRTs_mGMR).
%
%Author: D. Tanis PhD, 6/2024


%assume data is already loaded into the base workspace

%targ_data is a struct that hold information on which targets are available
%in each trial of each part of the experiment. The structure of targets_exp
%cells is 7 rows for each condition and 2 columns for targets present. The
%first column is the target that the subject is instructed to reach to, and
%the second column contains any additional targets that may be present (for
%split trials).
targ_data = struct;
targets_exp1 = {0,[];1,[];-1,[];3,[];-3,[];1,-3;-1,3};
targets_exp2 = {0,[];1,[];-1,[];3,[];-3,[];-3,1;3,-1};
targ_data.targets = {targets_exp1, targets_exp2};

%define the subjects, experiments, conditions, and trials to look for EOIs
%in:
sn = [1:length(data.subject)];
en = [1:length(data.subject(1).exp)];
cn = 'Perturb';
tn = 'Included';

%defin which trials are the jumps:
jumps = 2:5;

%set parameters that define a valid redirection:
acc_thresh = 42;
heading_thresh = 50;
dur = 60;

ymax = 18;   
r_thresh = 0.75; %radial distance indicating being inside the target

%set the channel to look at:
channel = 'Xacc_nobaseline';

%pull out all the data indicated:
temp = pullData(data, sn, en, cn, 'Included', 'Right_HandX');
tempj = pullData(data, sn, en, jumps, 'Included', 'Right_HandX');
[accs, listing] = pullData(data, sn, en, cn, tn, channel);
head = pullData(data, sn, en, cn, tn, 'Heading');
dhead = pullData(data, sn, en, cn, tn, 'Change in Heading');
xs = pullData(data, sn, en, cn, tn, 'Right_HandX');
ys = pullData(data, sn, en, cn, tn, 'Right_HandY');
ts = pullData(data, sn, en, cn, tn, 'Time_rel_Cue');

%set up data structures:
n_included = length(temp);  %total trials included
n_total = length(accs);     %total trials overall
j_total = length(tempj);    %total jumps overall
n_unidentifiable = 0;       %total number of trials in which a redirection could not be identified
j_unidentifiable = 0;       %same as above, for jumps only
n_overshoot = 0;            %number of trials in which the subject overshot the target
j_overshoot = 0;            %same as above, for jumps only
n_toomany = 0;              %number of trials in which there were > 2 redirections identified
j_toomany = 0;              %same as above, for jump trials
n_samedir = 0;              %number of trials in which there were two redirections oriented in the same direction
j_samedir = 0;              %same as above, for jumps only 


%step through each subject and experiment individually, and then through
%each condition:
sns = unique(listing(:,1))'; %find all subjects that are returned in the listing

%step through all available subjects:
for si = sns
    
    %step through each available experiment for this subject:
    s_match = listing(:,1)== si;
    ens = unique(listing(s_match,2))';
    for ei = ens
        
        %step through each available condition for this subject and
        %experiment:
        e_match = listing(:,2)==ei;
        cns = unique(listing(s_match&e_match,3))';
        for ci = cns
            
            %initialize the data:
            ntrials = length(data.subject(si).exp(ei).condition(ci).IsIncluded);
            data.subject(si).exp(ei).condition(ci).RTInd = cell(ntrials,2);
            data.subject(si).exp(ei).condition(ci).Direction = cell(ntrials,2);
            
            %get the indices for this specific subject/experiment/condition
            %combination:
            c_match = listing(:,3)==ci;
            keep = s_match&e_match&c_match;
            
            %access the acceleration and time profiles for these trials:
            these_accs = accs(keep);
            these_ts = ts(keep);
            these_xs = xs(keep);
            these_head = head(keep);
            these_dhead = dhead(keep);
            these_ys = ys(keep);
            this_list = listing(keep,:);
            
            %Look for redirections to the right:
            redir_inds_r = findRedirs(these_accs, these_ts, these_ys, these_dhead, acc_thresh, heading_thresh, dur);
            
            %Look for redirections to the left:
            redir_inds_l = findRedirs(invertTraces(these_accs), these_ts, these_ys, invertTraces(these_dhead), acc_thresh, heading_thresh, dur);
            
            %combine the right and left rtinds and remove any that dont
            %meet criteria:
            for ii = 1:length(these_accs)
                ti = this_list(ii,4);
                
                %identify the final location for this trial:
                y19 = find(these_ys{ii}>18,1,'first');
                all_targ = targ_data.targets(ei);
                these_targs = cell2mat(all_targ{1}(ci,:));
                end_samp = length(these_xs{ii});
                for jj = 1:length(these_targs)
                    %identify the final target choice based on the first
                    %target that the cursor is within range
                   new_end = find((sqrt((these_xs{ii}-these_targs(jj)).^2 + (these_ys{ii}-20).^2) <= r_thresh) == 1, 1, 'first');
                   if new_end < end_samp
                       end_samp = new_end;
                   end
                end
                
                %if no final target was achieved, broaden the search
                %radius:
                if end_samp == length(these_xs{ii})
                    for jj = 1:length(these_targs)
                        new_end = find((sqrt((these_xs{ii}-these_targs(jj)).^2 + (these_ys{ii}-20).^2) <= 1) == 1, 1, 'first');
                        if new_end < end_samp
                            end_samp = new_end;
                        end
                    end
                end
                
                %if still no target was achieved, set the end sample to the
                %time that the cursor crossed 95% of the way to the target:
                if end_samp == length(these_xs{ii})
%                     warning('No target acheived');
                     end_samp = y19;
                end
                
                %determine the final direction (left vs right) based on the
                %cursor position at the final sample:
                loc = these_xs{ii}(end_samp);
                if loc > 0
                    data.subject(si).exp(ei).condition(ci).Direction(ti,2) = {'R'};
                else
                    data.subject(si).exp(ei).condition(ci).Direction(ti,2) = {'L'};
                end
                
                
                redir_inds = [redir_inds_r{ii}, redir_inds_l{ii}];
                is_right = [ones(size(redir_inds_r{ii})), zeros(size(redir_inds_l{ii}))];
                
                %if no redirections were identified, add this to the tally
                %of unidentified trials and go to the next trial:
                if isempty(redir_inds)
                    n_unidentifiable = n_unidentifiable+1;
                    if any(ci == jumps)
                        j_unidentifiable = j_unidentifiable + 1;
                    end
                    continue
                end
                
               %sort redirections in temporal order:
               [redir_inds,I] = sort(redir_inds);
               is_right = is_right(I);
               
               %remove any inds that occur after outside the search
               %window:
               [~,t0] = min(abs(these_ts{ii}));
               yend = find(these_ys{ii}>ymax,1,'first');
               if isempty(yend)
                   yend = length(these_ys{ii});
               end
               keep = redir_inds<=yend & redir_inds >= t0+90 & redir_inds<=t0+500;
               redir_inds = redir_inds(keep);
               is_right = is_right(keep);
               
               %if there are no redirections left, add this to the tally
                %of unidentified trials and go to the next trial:
               if isempty(redir_inds)
                   n_unidentifiable = n_unidentifiable+1;
                   if any(ci == jumps)
                       j_unidentifiable = j_unidentifiable + 1;
                   end
                   continue
               end
               
               %verify that the last redirection to occur is congruent with
               %the final location of the trial. If it is not, remove it,
               %since this is an unwanted terminal redirection. (Some
               %subjects "sweep through" the targets and return to the
               %starting position rather than stopping in the target.)
               if length(redir_inds) > 1
               if logical(is_right(end)) ~= (loc > 0)
                   redir_inds = redir_inds(1:end-1);
                   is_right = is_right(1:end-1);
               end
               end
               
               %if there are no redirections left, add this to the tally
                %of unidentified trials and go to the next trial:
               if isempty(redir_inds)
                   n_unidentifiable = n_unidentifiable+1;
                   if any(ci == jumps)
                       j_unidentifiable = j_unidentifiable + 1;
                   end
                   continue
               end
               
               
               %exclude trials with more than 2 redirections:
               if length(redir_inds) > 2
                   data.subject(si).exp(ei).condition(ci).IsIncluded{ti} = false;
                   n_toomany = n_toomany +1;
                   if any(ci == jumps)
                       j_toomany = j_toomany + 1;
                   end
                   continue
               end
               
               %convert directions to a string:
               opts = {'L','R'};
               dirs_str = opts(is_right+1);
               
               %exclude trials with 2 redirections in the same direction:
               if length(redir_inds) == 2 && is_right(1) == is_right(2)
% % %                    data.subject(si).exp(ei).condition(ci).IsIncluded{ti} = false;
                    data.subject(si).exp(ei).condition(ci).RTInd{ti,1} = NaN;
                    data.subject(si).exp(ei).condition(ci).RTInd{ti,2} = NaN;
                    data.subject(si).exp(ei).condition(ci).Direction(ti,2) = dirs_str(2);
                    n_samedir = n_samedir +1;
                    if any(ci == jumps)
                       j_samedir = j_samedir + 1;
                   end
                   continue
               end
               
               %Because of the max y and max time cutoffs, there can be
               %redirections to targets that are unidentified because they
               %occured after this window. Thus a trajectory can end up at
               %a rightward target without an identifiable rightward
               %redirection. If a leftward redirection was identified in
               %this trial, we need to make a placeholder for the final
               %unidentifed redirection. This way we will know the final
               %direction of the trial, but we may not know the timing of
               %that final redirection. These trials are identified as
               %trials with a single redirection in the opposite direction
               %of where they end up (such as a trial that ends up at +3cm,
               %but has only a leftward redir detected). An issue arises with
               %trials that overshoot their target and then come back to
               %it. In those trials, the final location is also incongruent
               %with this "terminal redirection". (such as a trial that
               %ends up at +3cm, but there was an overshoot to +5cm, so
               %there was a leftward redirection to the +3cm). We don't
               %want to detect such overshoots. The following code
               %identifies the trials that we want as explained above,
               %based on the projected target that the trajectory will
               %achieve. If the trajectory will achieve a point further
               %away than the final target, it is considered an overshoot,
               %otherwise, it is considered a true redirection
               %get the target info
               all_targ = targ_data.targets(ei);
               these_targs = cell2mat(all_targ{1}(ci,:));
               targ_ind = find((these_targs>0)==(loc>0));
               if ~isempty(targ_ind)
                   targ_val = these_targs(targ_ind);
               else
                   targ_val = NaN;
               end
               
               %calculate the projected target:
               y_to_go = 20-these_ys{ii}(redir_inds(1));
               x_st = these_xs{ii}(redir_inds(1));
               proj_x = y_to_go*tand(these_head{ii}(redir_inds(1)))+x_st;
               
               if length(redir_inds) == 1 && (logical(is_right(1)) ~= (loc>0)) && ~isnan(targ_val)
                   if ((targ_val > 0) && (proj_x < targ_val)) || ...
                           ((targ_val < 0) && (proj_x > targ_val))
                       %the projected target is closer than the final endin
                       %location: this is a redir
                       redir_inds = [redir_inds, NaN];
                       is_right = [is_right, loc>0];
                       
                   else
                       %The projected target is further than the final
                       %target, thus this is an overshoot. Remove it:
                       redir_inds = [];
                       is_right = [];
%                        plot(these_xs{ii},these_ys{ii})
%                        ab = 2;
                   end
               end
               
               %if there are no redirections left, that is becuse an
               %overshoot redir was removed. Tally this, and move to the
               %next trial:
               if isempty(redir_inds)
                   n_unidentifiable = n_unidentifiable+1;
                   n_overshoot = n_overshoot+1;
                   if any(ci == jumps)
                       j_unidentifiable = j_unidentifiable + 1;
                       j_overshoot = j_overshoot+1;
                   end
                   continue
               end
               
               %save the identified redirections:
               opts = {'L','R'};
               dirs_str = opts(is_right+1);
               if length(redir_inds) == 2
                   data.subject(si).exp(ei).condition(ci).Direction(ti,:) = dirs_str;
                   data.subject(si).exp(ei).condition(ci).RTInd(ti,:) = num2cell(redir_inds);
               else
                   data.subject(si).exp(ei).condition(ci).Direction(ti,2) = dirs_str;
                   data.subject(si).exp(ei).condition(ci).RTInd{ti,2} = redir_inds;                   
               end
                
            end

        end
    end
end
%

function newtraces = invertTraces(traces)
%This function inverts each given traces. This is inversion is done so that
%the findRedirs function can also identify redirections to the left.

newtraces = traces;

if ~iscell(traces)
    newtraces = -traces;
else
    for ii = 1:length(traces)
        newtraces{ii} = -traces{ii};
    end
end

end

function redir_inds = findRedirs(accs, ts, ys, change_headings, a_thresh, head_thresh, dur)
%this function identifies redirections in acceleration data as points where
%the acceleration crosses above a certain threshold, remains above that
%threshold for a certain duration, and results in a certain velocity
%change.

redir_inds = cell(size(accs));
for ii = 1:length(accs)
    
    above = accs{ii}(1:end)>a_thresh | change_headings{ii}(1:end)>head_thresh;
    cross = [0; diff(above)];
    cross_up_inds = find(cross == 1);
    if isempty(cross_up_inds)
        continue
    end
    cross_down_inds = find(cross(cross_up_inds(1):end) == -1)+cross_up_inds(1)-1;
    cross_up_inds = cross_up_inds(1:length(cross_down_inds));
    
    count = 1;
    %step through each threshold crossing:
    for jj = 1:length(cross_up_inds)
        
        %skip this crossing if the duration above threshold is less than
        %dur:
        if cross_down_inds(jj) - cross_up_inds(jj) < dur
            continue
        end
        
        %if we get to this point, keep this point as a redirection:
        redir_inds{ii}(count) = cross_up_inds(jj);
        count = count+1;
    end
    
end

end