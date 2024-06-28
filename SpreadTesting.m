%SpreadTesting
%
%This script tests the performance of the oRT-detection algorithms against
%different levels of spread in the underlying oRTs. This script is designed
%to be run with the "plotResults.m" script which by default displays the
%mean error of the algorithms.
%
%Inputs:
%   - this script assumes that the canonical corrections have already been
%   generated and are in the base workspace (see GenerateCCs.m)
%   - this script assumes that the noise basis has already been gnerated
%   and is in the base workspace (see GetNoise.m)
%   - this script also assums that the experiment data structure is in the
%   base workspace and is labeled "data"
%
%Outputs:
%   - meanrtsout: 170 x 2 x 10 array of the mean oRT estimates for each
%   batch of trials and algorithm. The m dimension corresponds to all batches (10) for
%   each subjects (17). The n dimension corresponds to two levels of SNR.
%   The p dimension corresponds to different algorithms (not all columns
%   used)
%   - indiv_errors: 1700 x 2 x 10 array of individual oRT errors for each
%   batch of trials and for each algorithm. Similar structure as
%   meanrtsout, but each row is an individual trial, not a batch
%
%Author: D. Tanis PhD, 6/2024

%SNRs to test with:
test_sigx = 1;

%mean true oRT:
test_mean = 0.25;

%STDs of the true oRTs:
test_std = [0.01, 0.05];

%number of batches per subject:
n_batches_persub = 10;

%number of trials per batch:
n_per_group = 10;

%number of subjects:
n_subs = 17;

n_groups = n_subs*n_batches_persub;

%initialize data structures:
meanrtsout = nan(n_groups, length(test_std), 10);
SSErtsout = nan(n_groups, length(test_std), 10);
stdrtsout = meanrtsout;
perc_nans = SSErtsout;
peaks = cell(n_groups,length(test_std));
peakinds = peaks;   %note that peak inds are relative to the cue
meanpeaks = nan(n_groups,length(test_std));
meanpeakinds = meanpeaks;

outdat_all = nan(n_per_group*n_groups*length(test_std),9);
indiv_errors = nan(n_per_group*n_groups, length(test_std),10);


%generate a figure for plotting, as well as a wait bar:
f = figure();
h = waitbar(0,'Estimated remaining time:');
loops = 0;
tic;


subjs = repmat([1:n_subs]',n_batches,1);

group_mult = 1;

group_size = n_subs*group_mult;



outdat_count = 1;

%Step through each batch to test:
for jj = [1:n_groups]
    
    this_sub = subjs(jj);
    
    %generate the batch of trials to be tested:
    access_inds = repmat(this_sub.*ones(1,n_per_group)',group_mult,1);
    
    %get the accelerations, times, and rts for this subject:
    acc = acc_all(access_inds);
    t = t_all(access_inds);
    rts = rt_all(access_inds);
    
    %pull a random noise trace within each subject:
    noise_orig = cell(length(access_inds),1);
    for ii = 1:length(access_inds)
        noise_orig{ii} = noise_basis{access_inds(ii)}{randi(length(noise_basis{access_inds(ii)}), 1)};
    end
        
    %now loop through each STD to be tested:   
    for ii = 1:length(test_std)
    
        %first set all traces to begin at the test_mean:
        shift1 = (test_mean-rts);
        rts1 = rts+shift1;
        test_acc = indShift(acc, round(shift1*1000));
        noise = indShift(noise_orig, round(shift1*1000));    %noise traces are shifted as well
        noiseROC = indShift(noise_roc, round(shift1*1000));    %noise traces are shifted as well
        
        %now also add in a std if desired:
        shift2 = randn(length(rts1),1);
        shift2 = shift2-mean(shift2);
        std_mult = test_std(ii)/std(shift2);
        shift2 = shift2*std_mult;
        rts2 = rts1+shift2;
        
        %if there are any rts from the random generator are out of range,
        %don't shift those rts: (rare)
        badrts = rts2 < 0.05 | rts2 > 0.5;
        if any(badrts)
            shift2(badrts) = 0;
        end
        
        %perform final shifts:
        rts2 = rts1+shift2;
        rts2 = rts2-nanmean(rts2)+test_mean;
        test_acc = indShift(test_acc, round(shift2*1000));
        noise = indShift(noise, round(shift2*1000));
        noiseROC = indShift(noiseROC, round(shift2*1000));
        
        %now set the signals with desired amount of signal and noise:
        orig_test_acc1 = test_acc;
        orig_test_acc2 = test_acc;
        
        %step through each trial in the batch:
        for kk = 1:length(acc)
            
           %adjust the individual trial amplitude to replicate the
            %underlying data:
            adj = 1.07;
            test_acc{kk} = test_acc{kk}*adj;
            
            %multiply by the SNR:
            orig_test_acc2{kk} = test_acc{kk}*test_sigx(ii);
            
            %if the SNR is > 1, run this without noise:
            if test_sigx(ii) > 1
                test_acc{kk} = acc{kk}; %test without noise
            else
                test_acc{kk} = test_sigx(ii)*test_acc{kk}+noise{kk,1};  %run with noise:
            end
        end
        
        %now generate velocity profiles (for use with ROC):
        %note that a velocity offset is removed from each trace. This
        %ensures that the raw CC velocities are zero at the cue.
        %If this is not done, some CCs start at a high ROC value
        %because the velocity is initially high, which skews the ROC
        %predictions and makes them perform very poorly.
        test_v = test_acc;
        noise_v = noiseROC;
        for kk = 1:length(test_acc)
            test_v{kk} = cumsum(test_acc{kk})/1000-vel_offsets(this_sub);
            noise_v{kk} = cumsum(noiseROC{kk})/1000;
        end

        %get the peak values and locations:
        peaks{jj,ii} = [];
        peakinds{jj,ii} = [];
        thistraceclip = nan(length(test_acc),601);
        for kk = 1:length(test_acc)
            [~,ind0] = min(abs(t{kk}));
            [peaks{jj,ii}(kk), peakinds{jj,ii}(kk)] = max(test_acc{kk}(ind0+[0:600]));
            thistraceclip(kk,:) = test_acc{kk}(ind0+[0:600]);
        end
        
        thismeantrace = nanmean(thistraceclip);
        [meanpeaks(jj,ii),meanpeakinds(jj,ii)] = max(thismeantrace);
        
        %now test the algorithms:
        
        %Individual Threshold:
        vals = individualoRT(test_v, test_acc, t, 12.5, 100);

        meanrtsout(jj,ii,1) = nanmean(vals);
        stdrtsout(jj,ii,1) = nanstd(vals-rts2);
        SSErtsout(jj,ii,1) = nansae(vals, rts2);
        perc_nans(jj,ii,1) = sum(isnan(vals))/length(vals);
        indiv_errors((jj-1)*n_per_group+1+[0:n_per_group-1],ii,1) = vals-rts2;
        
        %Individual Regression:
        vals2 = individualRegression_testing(test_acc,t,100,20, orig_test_acc2, rts2);

        meanrtsout(jj,ii,2) = nanmean(vals2);
        stdrtsout(jj,ii,2) = nanstd(vals2-rts2);
        SSErtsout(jj,ii,2) = nansae(vals2, rts2);
        perc_nans(jj,ii,2) = sum(isnan(vals2))/length(vals2);
        indiv_errors((jj-1)*n_per_group+1+[0:n_per_group-1],ii,2) = vals2-rts2;

        %Grand Mean Regression:
        [gmean, tmean] = grandMean(test_acc, t);
        gm2 = individualRegression({gmean},{tmean},100, 20);
        meanrtsout(jj,ii,4) = nanmean(gm2);
        
        %Modified Grand Mean Regression (thresh-GMR):
        shift_1 = round(1000*(round(nanmean(vals),3)-vals));    %shift based on the individual threshold result
        shift_1(isnan(shift_1)) = 0;
        shifted_accs = indShift(test_acc, shift_1);
        [tsgmean1, tstmean1] = grandMean(indShift(test_acc, shift_1),t);
        meanrt1 = individualRegression({tsgmean1},{tstmean1},100,25);

        meanrtsout(jj,ii,5) = nanmean(meanrt1+vals-nanmean(vals));
        stdrtsout(jj,ii,5) = nanstd(meanrt1+vals-nanmean(vals)-rts2);
        SSErtsout(jj,ii,5) = nansae(meanrt1+vals-nanmean(vals), rts2);
        perc_nans(jj,ii,5) = sum(isnan(meanrt1+vals-nanmean(vals)))/length(vals);
        indiv_errors((jj-1)*n_per_group+1+[0:n_per_group-1],ii,5) = meanrt1+vals-nanmean(vals)-rts2;

        %ROC
        rocvals = ROC_RT(test_v, t, noise_v);
        meanrtsout(jj,ii,7) = nanmean(rocvals);

        %Canonical Correction Search
        [itvals2, outdat] = canonicalCorrectionSearch(test_acc, t, 100,1);
        outdat_all(outdat_count+[0:n_per_group-1],1:end-1) = outdat;
        outdat_all(outdat_count+[0:n_per_group-1],end) = test_sigx(ii);
       
        meanrtsout(jj,ii,9) = nanmean(itvals2);
        stdrtsout(jj,ii,9) = nanstd(itvals2-rts2);
        SSErtsout(jj,ii,9) = nansae(itvals2, rts2);
        perc_nans(jj,ii,9) = sum(isnan(itvals2))/length(itvals2);
        indiv_errors((jj-1)*10+1+[0:9],ii,9) = itvals2-rts2;
        
        
        outdat_count = outdat_count+10
        
        %plot results every 30 loops:
        if loops > 30
            if mod(loops-1,10) == 0
                plotResults(f);
            end
        end
        
        %calculate time remaining in simulation:
        loops = loops+1;
        loopsfrac = loops/(n_batches*length(test_sigx));
        telapsed = toc;
        tremaining = telapsed/loopsfrac - telapsed;
        
        %update the waitbar GUI:
        if tremaining > 60
            waitbar(loopsfrac,h, ['Estimated time remaining: ', num2str(round(tremaining/60,1)), ' minutes']);
        else
            waitbar(loopsfrac,h, ['Estimated time remaining: ', num2str(round(tremaining, 0)), ' seconds']);
        end
    end
    
    a = 2;
end


delete(h);
disp(' ');
disp(['Total time for simulation: ', num2str(round(toc/60,2)), ' minutes']);

plotResults(f);
 


function accs2 = indShift(accs, indshift)
%This function shifts the given data by the given number of indices from
%the indshift arrray. Empty data is padded with zeros. Positive shifts
%will shift towards later times

accs2 = accs;
for ii = 1:length(accs(:,1))
    for jj = 1:length(accs(1,:))
        pad = zeros(abs(indshift(ii)),1);
        
        if indshift(ii) > 0
            accs2{ii,jj} = [pad; accs{ii,jj}(1:end-indshift(ii))];
        elseif indshift(ii) < 0
            accs2{ii,jj} = [accs{ii,jj}(-indshift(ii)+1:end); pad];
        else
            accs2{ii,jj} = accs{ii,jj};
        end
    end
end


end

function [gmean, tmean] = grandMean(accs, ts)
%This function generates a grand mean of the given data, time-aligned to
%the cue

    starttime = -.1;
    entime = .6;
    acc_clip = nan(length(accs), 701);
    for ii = 1:length(accs)
        [~,st] = min(abs(ts{ii}-starttime));
        [~,en] = min(abs(ts{ii}-entime));
        acc_clip(ii,:) = accs{ii}(st:en);
    end

    gmean = nanmean(acc_clip)';
    tmean = starttime:0.001:entime;

end


function result = nansae(vals1, vals2)

    result = nanmean(abs(vals1-vals2));
 
end
