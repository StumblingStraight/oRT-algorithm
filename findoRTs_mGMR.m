%findoRTs_mGMR
%
%This script identifies individual oRTs in trials, after redirections and
%approximate oRTs have been identified by the identifyEpochofInterest
%script. This script uses the modified Grand Mean Regression (mGMR) method.
%
%Author: D. Tanis PhD, 6/2024

%assume data is in the base workspace

%set the subjects, experiments, conditions, and trials to look at:
sn = [1:length(data.subject)];
en = [1:length(data.subject(1).exp)];
cn = 'Perturb';
tn = 'All';

channel = 'Xacc_nobaseline';
search_window = [-50,+50];  %window around the estimated RT to search

%pull the data indicated:
[accs, listing] = pullData(data, sn, en, cn, tn, channel);
xs = pullData(data, sn, en, cn, tn, 'Right_HandX');
ys = pullData(data, sn, en, cn, tn, 'Right_HandY');
ts = pullData(data, sn, en, cn, tn, 'Time_rel_Cue');
rtinds = pullData(data, sn, en, cn, tn, 'RTInd');
dirs = pullData(data, sn, en, cn, tn, 'Direction');

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
            
            %get the indices for this specific subject/experiment/condition
            %combination:
            c_match = listing(:,3)==ci;
            keep = s_match&e_match&c_match;
            
            %access the acceleration and time profiles for these trials:
            these_accs = accs(keep);
            these_ts = ts(keep);
            these_xs = xs(keep);
            these_ys = ys(keep);
            these_rtinds = rtinds(keep,:);
            these_dirs = dirs(keep,:);
            this_list = listing(keep,:);
            
            %sort the RTinds into left vs right redirections:
            all_accs_r = cell(2*length(these_accs),1);
            all_accs_l = all_accs_r;
            all_ts_r = all_accs_r;
            all_ts_l = all_accs_r;
            rtinds_l = nan(2*length(these_accs),1);
            rtinds_r = rtinds_l;
            trial_inds_r = rtinds_l;
            trial_inds_l = rtinds_l;
            r_count = 1;
            l_count = 1;
            for ii = 1:length(these_rtinds)
                for jj = 1:2
                    if ~isempty(these_dirs{ii,jj})
                        if strcmp(these_dirs{ii,jj},'R') && ~isempty(these_rtinds{ii,jj}) && ~isnan(these_rtinds{ii,jj})
                            rtinds_r(r_count,1) = these_rtinds{ii,jj};
                            all_accs_r(r_count) = these_accs(ii);
                            all_ts_r(r_count) = these_ts(ii);
                            trial_inds_r(r_count) = ii;
                            r_count = r_count+1;
                        elseif strcmp(these_dirs{ii,jj},'L') && ~isempty(these_rtinds{ii,jj}) && ~isnan(these_rtinds{ii,jj})
                            rtinds_l(l_count,1) = these_rtinds{ii,jj};
                            all_accs_l(l_count) = these_accs(ii);
                            all_ts_l(l_count) = these_ts(ii);
                            trial_inds_l(l_count) = ii;
                            l_count = l_count+1;
                        end
                    end
                end
            end
            all_accs_r = all_accs_r(1:r_count-1);
            all_ts_r = all_ts_r(1:r_count-1);
            rtinds_r = rtinds_r(1:r_count-1,:);
            trial_inds_r = trial_inds_r(1:r_count-1);
            all_accs_l = all_accs_l(1:l_count-1);
            all_ts_l = all_ts_l(1:l_count-1);
            rtinds_l = rtinds_l(1:l_count-1,:);
            trial_inds_l = trial_inds_l(1:l_count-1);
            
            %first run a threshold identification on the trials:
            %do the right side first
            rtinds_r3 = [];
            if ~isempty(all_accs_r)
                rtinds_r2 = individualoRT(all_accs_r, all_ts_r, rtinds_r+search_window, 12);
                rtinds_r2(isnan(rtinds_r2)) = rtinds_r(isnan(rtinds_r2));
                
                %generate a time-shifted grand mean trace:
                acc_clip = nan(length(rtinds_r2),1000);
                for ii = 1:length(rtinds_r2)
                    this_data = all_accs_r{ii};
                    
                    %clip the data to the same size:
                    this_data = this_data(rtinds_r2(ii)-100:end);
                    if length(this_data) > 1000
                        this_data = this_data(1:1000);
                    end
                    acc_clip(ii,1:length(this_data)) = this_data;
                end
                %nanmean all the trials:
                grand_mean_r = {nanmean(acc_clip,1)'};
                gmean_t = [-100:1:900-1];
                
                %run a regression on the mean acceleration:
                window = [50,300];
                regRTind_r = individualRegression(grand_mean_r, {gmean_t}, window, 15, 100);
                
                if isnan(regRTind_r) || regRTind_r > window(2)
                    [~,pkind] = max(grand_mean_r{1}(1:window(2)));
                    regRTind_r = find(grand_mean_r{1}(1:pkind)<=0, 1, 'last');
                    if isempty(regRTind_r)
                        regRTind_r = 1;
                    end
                end
                
                RToffset_r = regRTind_r-100-1;
                rtinds_r3 = rtinds_r+RToffset_r;
                
                %remove any oRTs identified before 90 ms after the cue:
                for ii = 1:length(rtinds_r3)
                    [~,t90] = min(abs(all_ts_r{ii}-0.09));
                    if rtinds_r3(ii) < t90
                        rtinds_r3(ii) = NaN;
                    end
                end
                
            
            end
            
            %do the left side:
            rtinds_l3 = [];
            if ~isempty(all_accs_l)
                rtinds_l2 = individualoRT_search(invertTraces(all_accs_l), all_ts_l, rtinds_l+search_window, 10);
                rtinds_l2(isnan(rtinds_l2)) = rtinds_l(isnan(rtinds_l2));

                %generate a time-shifted grand mean trace to the left:
                acc_clip = nan(length(rtinds_l2),1000);
                for ii = 1:length(rtinds_l2)
                    this_data = all_accs_l{ii};
                    
                    this_data = this_data(rtinds_l2(ii)-100:end);
                    if length(this_data) > 1000
                        this_data = this_data(1:1000);
                    end
                    acc_clip(ii,1:length(this_data)) = this_data;
                end
                grand_mean_l = invertTraces({nanmean(acc_clip,1)'});
                gmean_t = [-100:1:900-1];
                
                regRTind_l = individualRegression_search(grand_mean_l, {gmean_t}, window, 100, 15);
                
                if isnan(regRTind_l) || regRTind_l > window(2)
                    [~,pkind] = max(grand_mean_l{1}(1:window(2)));
                    regRTind_l = find(grand_mean_l{1}(1:pkind)<=0, 1, 'last');
                    if isempty(regRTind_l)
                        regRTind_l = 1;
                    end
                end
                
                RToffset_l = regRTind_l-100-1;
                rtinds_l3 = rtinds_l+RToffset_l;
                
                for ii = 1:length(rtinds_l3)
                    [~,t90] = min(abs(all_ts_l{ii}-0.09));
                    if rtinds_l3(ii) < t90
                        rtinds_l3(ii) = NaN;
                    end
                end
            end
            
            %now update the data structure with the new rts:
            rcount = 1;
            lcount = 1;
            for ii = 1:length(these_accs)
                for jj = 1:2
                    if ~isempty(these_dirs{ii,jj})
                        if strcmp(these_dirs{ii,jj},'R') && ~isempty(these_rtinds{ii,jj}) && ~isnan(these_rtinds{ii,jj})
                            if ~isempty(rtinds_r3(rcount))% && ~isnan(rtinds_r3(rcount))
                                data.subject(si).exp(ei).condition(ci).RTInd{ii,jj} = rtinds_r3(rcount);
%                                 temp = all_accs_r{rcount}(rtinds_r3(rcount));
                                rcount = rcount + 1;
                            end
                        elseif strcmp(these_dirs{ii,jj},'L') && ~isempty(these_rtinds{ii,jj}) && ~isnan(these_rtinds{ii,jj})
                            if ~isempty(rtinds_l3(lcount))% && ~isnan(rtinds_l3(lcount))
                                data.subject(si).exp(ei).condition(ci).RTInd{ii,jj} = rtinds_l3(lcount);
%                                 temp = all_accs_l{lcount}(rtinds_l3(lcount));
                                lcount = lcount + 1;
                            end
                        end
                    end
                end
            end
            
            ab = 2;
        end
    end
end


function newtraces = invertTraces(traces)
%This function inverts each given trace.
newtraces = traces;

if ~iscell(traces)
    newtraces = -traces;
else
    for ii = 1:length(traces)
        newtraces{ii} = -traces{ii};
    end
end

end