function RTs = canonicalCorrectionSearch(accs, ts, fc, error_metric)
%This function runs the canonicalCorrectionSearch oRT detection algorithm.
%This algorithm identifies an oRT in a given trial by finding the best
%match of an estimated canonical correctin to the given trial. First, a
%common threshold crossing is identified in all trial. Next, all the
%trials are shifted such that these threshold crossings align and an
%average of all the trials is taken. This time-aligned grand average is the
%estimated canonical correction (est-CC). Next, for each trial, the est-CC
%is scaled in amplitude and width to approximately match the given trial.
%Then, the est-CC is moved in time such that the error between the given
%trial and the est-CC is minimized. The point where the error is at a
%minimum is identified as the oRT.
%
%Inputs:
%   - accs: an mx1 cell array of acceleration traces
%   - ts:   an mx1 cell array of time data series
%   - fc:   cuttof frequency for a low-pass filter. If >= 100, no filtering
%   will be performed
%   - error_metric:  this specifies the error metric to use for the algorithm
%
%Outputs:  
%   - RTs: an mx1 array of the oRT times detected
%
%Author: D. Tanis PhD, 6/2024



RTs = nan(size(accs));


%first filter the data to the desired Fc, if indicated:
fs = 1000;
n_poles = 3;
[b,a] = butter(n_poles,fc/(fs/2));

accs_filt = accs;
if fc < 100
for ii = 1:length(accs)
    accs_filt{ii} = filtfilt(b,a,accs{ii});
end
end

%Next, build the estimated canonical correction:
%First, get an estimate of each individual trial rt by the threshold method:
RTs_est = individualoRT(accs,ts,20,100);

%Next, time-align the traces to their threshold-crossings and average them
%together:
clip_canon = nan(length(accs),1000);
for ii = 1:length(accs)
   [~,cue] = min(abs(ts{ii}));
   if isnan(RTs_est(ii))
       continue
   end
   st = round(RTs_est(ii)*1000) + cue;
   subdata = accs_filt{ii}(st-200:st+600);
   clip_canon(ii,1:length(subdata)) = subdata;
end

%do the averaging:
if length(accs) > 1
    canon = nanmean(clip_canon);
else
    canon = clip_canon;
end

%now perform a regression on the estimated canonical trace to identify its RT:
canonrtind = individualRegression({canon},{[0:0.001:1-0.001]},100,20);

%clip the canon to begin at its oRT:
if ~isnan(canonrtind) && canonrtind > 0
    canon = canon(round(canonrtind*1000):end);
end

%Now run the optimization process to match the estimated Canonical to each
%trial and identify the best fit
%step through each trial:
for ii = 1:length(accs)
    
    %clip this trial to a large window around the cue:
    [~,ind0] = min(abs(ts{ii}));
    thistrial = accs{ii}(ind0:ind0+600);
    len = length(thistrial);
    
    thisjerk = [0, diff(thistrial)];
    
    %identify the peak of this trial:
    %Note that index 1 is the time of the cue
    st = 1;
    temp = thistrial;
    temp(abs(thisjerk) > 0.5) = 0;
    [peak, peakind] = max(temp(50:end));
    peakind = peakind + 50-1;
    
    %normalize the estimated canonical to be the same size as the peak in the current trial:
    if peak > 20
        thiscanon = canon/max(canon) * peak;
    end
    
    %also adjust the estimated canonical such that its slope matches the
    %slope of the current trial:
    adjust_w = true;
    if adjust_w
        
        %apply a low pass filter to both trials to reduce error in the
        %slope estimate:
        thistrial2 = thistrial;
        fs = 1000;
        fc = 10;
        n_poles = 3;
        [b,a] = butter(n_poles,fc/(fs/2));
        thistrial2 = filtfilt(b,a,thistrial2);
        thiscanon2 = thiscanon(~isnan(thiscanon));
        thiscanon2 = filtfilt(b,a,thiscanon2);
        
        %find the peak of the trial:
        temp2 = thistrial2;
        thisjerk2 = [0,diff(temp2)];
        temp2(abs(thisjerk2)>0.5) = 0;
        [peak2, peakind2] = max(temp2(50:end));
        peakind2 = peakind2 + 50 - 1;
        
        %find the last points before the peak that the trial crosses 20%
        %and 80% of the peak:
        ind20 = find(thistrial2(1:peakind2)< 0.2*peak2, 1,'last');
        if ~isempty(ind20)
        ind80 = find(thistrial2(ind20:peakind2)< 0.8*peak2, 1,'last') + ind20-1;
        end

        %do the same for the estimated canonical correction:
        [cpeak,cpkind] = max(thiscanon2);
        cind20 = find(thiscanon2(1:cpkind)< 0.2*cpeak, 1,'last');
        cind80 = find(thiscanon2(1:cpkind)< 0.8*cpeak, 1,'last');

        if isempty(ind20) || isempty(cind20) || isempty(ind80)
            %if any of the required data was not found, don't adjust the
            %slope and leave it as it was:
            adjust = 1;
        elseif isnan((cind80-cind20)/(ind80-ind20))
            %if any of the required data was not found, don't adjust the
            %slope and leave it as it was:
            adjust = 1;
        else
            %calculate the slope adjustment:
            adjust = (cind80-cind20)/(ind80-ind20);
        end

        %limit the slope adjustment to +/-50%:
        if adjust > 1.5
            adjust = 1.5;
        elseif adjust < 0.5
            adjust = 0.5;
        end
        
        %perform slope adjustment:
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
    

    [pk,pkind] = max(thiscanon);
    tempind = find(thiscanon(1:pkind)<0.99*pk,1,'last');
    pkind = tempind;
    
    %Set the search bounds for the algorithm. The optimization method will
    %search around the threshold crossing in a window of +/- 75 ms.
    thresh = 20;
    estpt = find(thistrial(st:peakind)<thresh, 1,'last')+st-1;
    can_estpt = find(thiscanon(1:pkind)<thresh,1,'last');
    if isempty(can_estpt)
        can_estpt = 1;
    end
    estpt = estpt-can_estpt;

    if isempty(estpt)
        estpt = peakind-50;
    end
    start_ind = estpt-75;
    end_ind = estpt+75;

    if start_ind < st
        start_ind = st;
    end
    if end_ind < start_ind + 20
        end_ind = start_ind + 100;
    end
    
    %remove any negative values in this trial and in the estimated
    %canonical. Also, remove any values where the slope is negative (ie
    %second half of the acceleration):
    [tpk,tpkind] = max(thistrial);
    eval_trial = thistrial;
    
    eval_trial(eval_trial<0) = 0;
    eval_trial(thisjerk<0) = 0;
    eval_trial(tpkind:end) = tpk;
    
    eval_canon = thiscanon;
    eval_canon(eval_canon<0) = NaN;
    eval_canon(cjerk<0) = NaN;
    
    %to improve processing speed, each trial is quantized before running
    %the optimization:
    [trial_levels,canon_levels] = quantization(eval_trial, eval_canon, 100);

    %run a genetic algorithm optimization to match the estimated canonical
    %to the given trial. The point where there is minimal error defines the
    %location of the oRT:
    rtinds = ga(@(x)eval_indiv(x, eval_trial, eval_canon, len, trial_levels, canon_levels, error_metric),1,[],[],[],[],start_ind,end_ind,[],1);
    
    rts_comb = rtinds';
    RTs(ii) = rts_comb./1000;

end

end

function f = eval_indiv(inds, test_acc, canon, maxlen, trial_levels, canon_levels, error_metric)
%This function is a function called by the ga algorithm, and provides an
%error metric for running the algorithm. Here,the error metric is a measure
%of error between the given trial and the estimated canonical correction
%with the provided amount of time-shift.
%
%Several error metrics were experimented with, and the error metric used
%can be set by the input "error_metric."
%   1: Nearest-neighbor
%   2: Sum absolute error
%   3: Mean absolute error
%   4: Mean squared error


if error_metric~=1
    test_acc = [test_acc, zeros(1,150)];    %zero pad just in case

    generated_acc = nan(size(test_acc));
    generated_acc(inds+[0:length(canon)-1]) = canon;

    generated_acc = generated_acc(1:length(test_acc));

    if error_metric == 2
        
        %run with a sum absolute error
        err = abs(test_acc-generated_acc);

        f = nansum(err);    %error function
        return
        
    elseif error_metric == 3
        
        %run with a mean absolute error
        err = abs(test_acc-generated_acc);

        f = nanmean(err);    %error function
        return
        
    elseif error_metric == 4
        
        %run with the mean squared error:
        err = (test_acc-generated_acc).^2;

        f = nanmean(err);    %error function
        return
        
    end


else
    %run with the nearest neighbor metric:
    err = 0;
    [m,~] = size(canon_levels);
    for ii = 1:m
        %step through each quantized amplitude level
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

end
  
function [trial_levels, canon_levels] = quantization(trial, canon, nlevels)
%This function quantizes the given trial and canon into n_levels of amplitude and
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