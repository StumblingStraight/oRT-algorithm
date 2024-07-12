%findoRTs_CCS
%
%This script identifies individual oRTs in trials, after redirections and
%approximate oRTs have been identified by the identifyEpochofInterest
%script. This script uses the Canonical Correction Search (CCS) method.
%
%Author: D. Tanis PhD, 6/2024

%assume data is in the base workspace

%set the subjects, experiments, conditions, and trials to look at:
sn = [1:length(data.subject)];
en = [1:length(data.subject(1).exp)];
% en = [2:3];
cn = 'Perturb';
tn = 'All';

%structure of which trials to pull template redirections from:
template_cns_left = [NaN, NaN, 3, NaN, 5, 5, 3];
template_cns_right = [NaN, 2, NaN, 4, NaN, 2, 4];

channel = 'Xacc_nobaseline';
search_window = [-50,+100];  %window around the estimated RT to search

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

%setup a waitbar GUI:
h = waitbar(0,'Estimated remaining time:');
loops = 0;
tic

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
            
            %sort the RTinds into left vs right rts:
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
                        if strcmp(these_dirs{ii,jj},'R')&& ~isempty(these_rtinds{ii,jj}) && ~isnan(these_rtinds{ii,jj})
                            rtinds_r(r_count,1) = these_rtinds{ii,jj};
                            all_accs_r(r_count) = these_accs(ii);
                            all_ts_r(r_count) = these_ts(ii);
                            trial_inds_r(r_count) = ii;
                            r_count = r_count+1;
                        elseif strcmp(these_dirs{ii,jj},'L')&& ~isempty(these_rtinds{ii,jj}) && ~isnan(these_rtinds{ii,jj})
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
            
            %first identify oRTs in redirections to the right:
            rtinds_r2 = [];
            if ~isempty(all_accs_r)
                templ_accs = [];
                templ_ts = [];
                %generate a template by pulling the corresponding jump
                %trial redirections:
                if ~isnan(template_cns_right(ci))
                    %pull the acceleration traces for the template:
                    templ_accs = pullData(data, si, ei, template_cns_right(ci), 'Included',channel);
                    templ_ts = pullData(data, si, ei, template_cns_right(ci), 'Included','Time_rel_Cue');
                    templ_rtinds = pullData(data, si, ei, template_cns_right(ci), 'Included','RTInd');
                    templ_dirs = pullData(data, si, ei, template_cns_right(ci), 'Included','Direction');
                    matches = cellfind(templ_dirs(:,2),'R') & ~cellempty(templ_rtinds(:,2));
                    templ_rtinds = cell2mat(templ_rtinds(matches,2));
                    templ_ts = templ_ts(matches);
                    templ_accs = templ_accs(matches);
                end
                
                %run the CCS algorithm using this template:
                [rtinds_r2] = CCS(all_accs_r, all_ts_r, templ_accs, templ_ts, templ_rtinds, rtinds_r+search_window);
                rtinds_r2(isnan(rtinds_r2)) = rtinds_r(isnan(rtinds_r2));
            
                %remove any rts that are too early:
                for ii = 1:length(rtinds_r2)
                    if ~isnan(rtinds_r2(ii)) && all_ts_r{ii}(rtinds_r2(ii)) < 0.09
                        rtinds_r2(ii) = nan;
                    end
                end
            end
            
            %now identify oRTs to the left:
            rtinds_l2 = [];
            if ~isempty(all_accs_l)
                templ_accs = [];
                templ_ts = [];
                if ~isnan(template_cns_left(ci))
                    %pull the acceleration traces for the template:
                    templ_accs = pullData(data, si, ei, template_cns_left(ci), 'Included',channel);
                    templ_ts = pullData(data, si, ei, template_cns_left(ci), 'Included','Time_rel_Cue');
                    templ_rtinds = pullData(data, si, ei, template_cns_left(ci), 'Included','RTInd');
                    templ_dirs = pullData(data, si, ei, template_cns_left(ci), 'Included','Direction');
                    matches = cellfind(templ_dirs(:,2),'L') & ~cellempty(templ_rtinds(:,2));
                    templ_rtinds = cell2mat(templ_rtinds(matches,2));
                    templ_ts = templ_ts(matches);
                    templ_accs = templ_accs(matches);
                end
                
                [rtinds_l2] = CCS(invertTraces(all_accs_l), all_ts_l, invertTraces(templ_accs), templ_ts, templ_rtinds, rtinds_l+search_window);
                rtinds_l2(isnan(rtinds_l2)) = rtinds_l(isnan(rtinds_l2));
                
                %remove any rts that are too early:
                for ii = 1:length(rtinds_l2)
                    if ~isnan(rtinds_l2(ii)) && all_ts_l{ii}(rtinds_l2(ii)) < 0.09
                        rtinds_l2(ii) = nan;
                    end
                end
            end
            
            %now update the data structure with the new rts:
            rcount = 1;
            lcount = 1;
            for ii = 1:length(these_accs)
                for jj = 1:2
                    if ~isempty(these_dirs{ii,jj})
                        if strcmp(these_dirs{ii,jj},'R')&& ~isempty(these_rtinds{ii,jj}) && ~isnan(these_rtinds{ii,jj})
                            if ~isempty(rtinds_r2(rcount))
                                data.subject(si).exp(ei).condition(ci).RTInd{ii,jj} = rtinds_r2(rcount);
                                rcount = rcount + 1;
                            end
                        elseif strcmp(these_dirs{ii,jj},'L')&& ~isempty(these_rtinds{ii,jj}) && ~isnan(these_rtinds{ii,jj})
                            if ~isempty(rtinds_l2(lcount))
                                data.subject(si).exp(ei).condition(ci).RTInd{ii,jj} = rtinds_l2(lcount);
                                lcount = lcount + 1;
                            end
                        end
                    end
                end
            end
            
            %update the waitbar GUI:
            loops = loops+1;
            loopsfrac = loops/(length(sns)*length(ens)*length(cns));
            telapsed = toc;
            tremaining = telapsed/loopsfrac - telapsed;
            
            if tremaining > 60
                waitbar(loopsfrac,h, ['Estimated time remaining: ', num2str(round(tremaining/60,1)), ' minutes']);
            else
                waitbar(loopsfrac,h, ['Estimated time remaining: ', num2str(round(tremaining, 0)), ' seconds']);
            end
        end
    end
end

disp(['Time elapsed: ', num2str(round(toc/60,2)), ' min']);
delete(h);

function newtraces = invertTraces(traces)
%This function inverts each given trace
newtraces = traces;

if ~iscell(traces)
    newtraces = -traces;
else
    for ii = 1:length(traces)
        newtraces{ii} = -traces{ii};
    end
end

end

function logical = cellfind(data, query)
%this function allows a single function to be used for finding data within
%a cell. Queries can be both strings and numeric.


if isnumeric(query)
    logical = false(length(data(:,1)),length(query));
    for ii = 1:length(data)
        for jj = 1:length(query)
            if ~isempty(data{ii}) && data{ii} == query(jj)
                logical(ii,jj) = true;
            end
        end
    end
elseif ischar(query)
    logical = false(size(data));
    for ii = 1:length(data)
        if strcmp(data{ii},query)
            logical(ii) = true;
        end
    end
end


end

function output = cellempty(data)
%This function determines if the cells in data are empty or not

output = false(size(data));
for ii = 1:length(data)
   if isempty(data{ii})
       output(ii) = true;
   end
end

end

function [RTinds, varargout] = CCS(accs, ts, templ_accs, templ_ts, templ_rtinds, search_window, varargin)
%This function performs the canonical search oRT identification procedure.
%
%Varargout is an array with a row for each trial and the following columns:
% Col 1: error metric value
% Col 2: peak of the trial acceleration
% Col 3: peak of the adjusted canonical
% Col 4: slope adjustment value

varargout{1} = nan(length(accs),5);

RTinds = nan(size(accs));

do_adj = 1;
if ~isempty(varargin)
    do_adj = varargin{1};
end

win = [-25:75]; %window around the estimated RT to look for the true RT
arg = 1;        %determines which error metric to run

try

    
if isempty(templ_rtinds)
    %get an estimate of the template traces' RTs
    templ_rtinds = individualoRT([],templ_accs, templ_ts,12.5,100);
end

%generate the canonical trace by time-aligning all traces to their rts and
%then averaging them together:
clip_canon = nan(length(templ_accs),801);
for ii = 1:length(templ_accs)
   if isnan(templ_rtinds(ii))
       %skip any traces that do not have an identified rt
       continue
   end
   st = templ_rtinds(ii)-200;
   en = templ_rtinds(ii)+600;
   len = length(templ_accs{ii});
   if en > len
       en = len;
   end
   subdata = templ_accs{ii}(st:en);
   clip_canon(ii,1:length(subdata)) = subdata;
end
canon = nanmean(clip_canon,1);

%now perform a regression on the canonical trace to identify its RT:
canonrtind = individualRegression_search({canon'},{[0.001:0.001:0.801]},[100,500],100,25);

%correct the identified in if needed:
if isnan(canonrtind) || canonrtind == 1 || canonrtind > 500
    %find the peak and backtrack to the zero-crossing:
    [~,pkind] = max(canon(1:500));
    canonrtind = find(canon(1:pkind)<=0, 1,'last');
    if isempty(canonrtind)
        canonrtind = 1;
    end
end

%clip the canon to begin at its oRT:
if ~isnan(canonrtind) && canonrtind > 0
    canon = canon(canonrtind:end);
end
    

%now run the optimization on each trial individually:
for ii = 1:length(accs)
    
    if isempty(search_window)
        %no search window: set one
        [~,t0] = min(abs(ts{ii}-0.09));
        search_inds = t0+[0:510];
    elseif isempty(search_window(ii))
        %skip this trial if it has no estimated rtinds given
        continue
    else
        [~,t90] = min(abs(ts{ii}-0.09));
        startind = search_window(ii,1);
        if startind < t90
            startind = t90;
        end
        search_inds = startind:search_window(ii,2);
    end

    thistrial = accs{ii};
    thistrial = thistrial(~isnan(thistrial));
    
    thisjerk = [0; diff(thistrial)];
    
    %identify the peak of this trial within the given window:
    temp = thistrial;
    temp(abs(thisjerk) > 0.5) = 0;
    [peak, peakind] = max(temp(search_inds));
    peakind = peakind + search_inds(1) -1;
    
    %normalize the canon to be the same size as the peak in question:
    if peak > 20 && do_adj
        thiscanon = canon/max(canon) * peak;
    else
        %do not normalize:
        thiscanon = canon;
    end
    
    %now adjust the canon such that its slope matches the slope of the
    %trial:
    adjust = 1;
    cpeak = 0;
    if do_adj

        ind20 = find(thistrial(1:peakind)< 0.2*peak, 1,'last');
        if isempty(ind20)
            ind20 = find(thistrial(1:peakind)< 0.3*peak, 1,'last');
        end
        if ~isempty(ind20)
            ind80 = find(thistrial(ind20:peakind)< 0.8*peak, 1,'last') + ind20-1;
        end

        [cpeak,cpkind] = max(thiscanon);
        cind20 = find(thiscanon(1:cpkind)< 0.2*cpeak, 1,'last');
        cind80 = find(thiscanon(1:cpkind)< 0.8*cpeak, 1,'last');
        
        if isempty(ind20) || isempty(cind20) || isempty(ind80)
            adjust = 1;
        else
            adjust = (cind80-cind20)/(ind80-ind20);
        end
        
        %set bounds on the slope adjustment:
        if adjust > 1.5
            adjust = 1.5;
        elseif adjust < 0.5
            adjust = 0.5;
        end
        
        %resample the canonical to match the slope of the trial:
        [p,q] = rat(1/adjust);
        canon_orig = thiscanon;
        thiscanon = resample(thiscanon,p,q);
        if 1/adjust >= 1
            thiscanon = thiscanon(1:length(canon_orig));
        else
            thiscanon = [thiscanon, nan(1,length(canon_orig)-length(thiscanon))];
        end
    end
    
    cjerk = [0, diff(thiscanon)];
    cjerk(1) = cjerk(2);
    
    %set the search bounds for the algorithm
    start_ind = search_inds(1);
    end_ind = search_inds(end);
    
    %remove any parts of the canon or trial where the acceleration or the
    %jerk are < 0:
    eval_trial = thistrial;
    [tpk,tpkind] = max(eval_trial);
    eval_trial(eval_trial<0) = 0;
    eval_trial(thisjerk<0) = 0;
    eval_trial(tpkind:end) = tpk;
    
    eval_canon = thiscanon;
    eval_canon(eval_canon<0) = NaN;
    eval_canon(cjerk<0) = NaN;
    
    [trial_levels,canon_levels] = quantization(eval_trial, eval_canon, 100);

    %run a genetic algorithm optimization to match the canonical to the
    %trace
    [rtind,errval] = ga(@(x)eval_indiv(x, eval_trial, eval_canon, trial_levels, canon_levels, arg),1,[],[],[],[],start_ind,end_ind,[],1);

    RTinds(ii) = rtind;

    
    varargout{1}(ii,1) = errval;
    varargout{1}(ii,2) = peak;
    varargout{1}(ii,3) = cpeak;
    varargout{1}(ii,4) = adjust;

end

catch ME
    aa = 2;
end

end

function f = eval_indiv(inds, test_acc, canon, trial_levels, canon_levels, arg)
        %this function assesses the error between the mean acceleration trace
        %and the generated acceleration trace from the canonical correction.
        try
            
        if arg~=1
            test_acc = [test_acc, zeros(1,150)];    %zero pad just in case
            
            generated_acc = nan(size(test_acc));
            generated_acc(inds+[0:length(canon)-1]) = canon;
       
            generated_acc = generated_acc(1:length(test_acc));
        
            if arg ~=6
                %run with a sum absolute error
                err = abs(test_acc-generated_acc);
        
                f = nansum(err);    %error function
                return
            elseif arg == 6
                %run with a mean absolute error
                err = abs(test_acc-generated_acc);

                f = nanmean(err);    %error function
                return
            end

        else
            %do nearest neighbor metric
            err = 0;
            [m,~] = size(canon_levels);
            for ii = 1:m
                %step through each pre-defined level
                if isempty(canon_levels{ii})
                    continue
                end
                n_vals = length(canon_levels{ii});
                for jj = 1:n_vals
                   %step through each canon value in this level
                   opts = trial_levels{ii};
                   thisval = canon_levels{ii}(jj)+inds-1;
                   if ~isempty(opts)
                       %find the closest point in the trial:
                       thiserr = min(abs(opts-thisval));
                       err = err+thiserr;
                   else
                       err = err;
                   end
                end
            end
        end
            

        f = err;
        
        catch ME
            aa = 2;
        end
        
end
  
function [trial_levels, canon_levels] = quantization(trial, canon, nlevels)
%This function quantizes the given trial and canon into n_levels and
%identifies the indices of data located at each level

pk = max(max(canon),max(trial));

trial_levels = cell(nlevels,1);
canon_levels = trial_levels;

vec_t = 1:length(trial);
vec_c = 1:length(canon);

for ii = 1:nlevels
   low = (pk/nlevels)*(ii-1);
   hi = (pk/nlevels)*(ii);
   
   %find matching data for trial
   matches = trial>low & trial <=hi;
   inds = vec_t(matches);
   if ~isempty(inds)
       trial_levels{ii} = inds;
   end
   
   %find matching data for canon
   matches = canon>low & canon<=hi;
   inds = vec_c(matches);
   if ~isempty(inds)
       canon_levels{ii} = inds;
   end
    
end


end