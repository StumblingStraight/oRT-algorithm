%RemoveBaselineAcc

%This script removes the tailored baseline acceleration from each trial. By
%default, this script performs this as done in Kurtzer 2020, based on the
%lateral acceleration only. However, this script can also include other
%variables, such as lateral position and heading to better match the
%trajectories.
%
%The algorithm steps through each trial that is marked as included. For
%each non-baseline trial, the algorithm ranks the baseline trials based on
%similarity to the trial in question The most similar baseline
%trajectories are averaged and this average is removed from the trial in
%question.
%
%Before running this script, the data should be loaded into the MATLAB
%workspace, in a variable named "data"

%pull something from all included trials to determine total number of
%trials that are included:
n_exp = length(data.subject(1).exp);
n_sub = length(data.subject);
try
    temp = pullData(data, [1:n_sub],[1:n_exp],'Perturb','included','condition');
catch
    temp = pullData(data, [1:12],[1:2],'Perturb','included','condition');
end

header = {'Subject','Experiment','Condition','Trial','Dir 1','Dir 2','RT 1','RT 2'};
outdata = nan(length(temp), length(header));
outdata = num2cell(outdata);
outcount = 1;

tx = 0.0;       %time to look for RT
win = 90;         %window size to look at
thresh = 6;   %threshold for maximum tolerable score

targx_coord = [-5, -3, -1, 1, 3, 5];    %x location of targets
targx_coord = [-10:1:10];    %x location of targets
yy = 20;    %y location of targets
r = 0.75;    %radius of target

for sn = 1:n_sub
    disp(['Processing Subject ', num2str(sn)]);
    for en = 1:length(data.subject(1).exp)
        %assemble a matrix of headings and x postions at +100 ms from
        %baseline trials:
        en_b = en;
        ts_b = pullData(data, sn, en_b, 1, 'included','Time_rel_Cue');
        xs_b = pullData(data, sn, en_b, 1, 'included','Right_HandX');
        heads_b = pullData(data, sn, en_b, 1, 'included','Heading');
        diff_heads_b = pullData(data, sn, en_b, 1, 'included','Change in Heading smoothed');
        xacc_b = pullData(data, sn, en_b, 1, 'included','Right_HandX_Acc');
%         xacc_b = pullData(data, sn, en_b, 1, 'included','Xacc_nobaseline');
        xacc_bns = pullData(data, sn, en_b, 1, 'included','Right_HandX_Acc');
%         xacc_bns = pullData(data, sn, en_b, 1, 'included','Xacc_nobaseline');
        factors = cell(length(ts_b),3);
        norm_factors = zeros(2,3);
        lengths = nan(1,length(xs_b));  %vector containing lengths of each baseline trial
        start_inds = nan(1,length(xs_b)); %vector containing t = 0 points of each baseline trial
        
        %assemble a cell of all the timepoints of baseline trials to
        %consider. While doing this, identify the max and min of each data
        %type for normalization purposes.
        for ii = 1:length(ts_b)
            [~,ind] = min(abs(ts_b{ii}-tx));
            factors{ii,1} = xs_b{ii}(ind+[1:win]);
            factors{ii,2} = heads_b{ii}(ind+[1:win]);
            factors{ii,3} = xacc_b{ii}(ind+[1:win]);
            
            for jj = 1:3
                if max(factors{ii,jj}) > norm_factors(1,jj)
                    norm_factors(1,jj) = max(factors{ii,jj});
                end
                if min(factors{ii,jj}) < norm_factors(2,jj)
                    norm_factors(2,jj) = min(factors{ii,jj});
                end
            end
            
            lengths(ii) = length(xs_b{ii});
            start_inds(ii) = ind-tx*1000-100;   %start ind is taken as 100 ms BEFORE the t = 0
        end
        
        %normalize the data:
        [ntrials, nfacs] = size(factors);
        
        for ii = 1:ntrials
           for jj = 1:nfacs 
               factors{ii,jj} = (factors{ii,jj} - norm_factors(2,jj))./(norm_factors(1,jj)-norm_factors(2,jj));
           end
        end
        
        
        %now pull all the non-baseline trials:
        conds = 2:length(data.subject(sn).exp(en).condition);
%         conds = 1:length(data.subject(sn).exp(en).condition);
%         conds = 1;
        trial_type = 'all';
        [ts, listing] = pullData(data, sn, en, conds, trial_type,'Time_rel_Cue');
        xs = pullData(data, sn, en, conds, trial_type,'Right_HandX');
        ys = pullData(data, sn, en, conds, trial_type,'Right_HandY');
%         heads = pullData(data, sn, en, conds, trial_type,'Heading');
        heads = pullData(data, sn, en, conds, trial_type,'Right_HandX_Vel');
        diffheads = pullData(data, sn, en, conds, trial_type,'Change in Heading smoothed');
        xaccs = pullData(data, sn, en, conds, trial_type,'Right_HandX_Acc_smoothed');
        xaccsns = pullData(data, sn, en, conds, trial_type,'Right_HandX_Acc');
        
        %copy over the subject, exp, cond, trial info:
        outdata(outcount+[1:length(listing)],1:4) = num2cell(listing);
        
        
        
        %step thru each trial and remove the average baseline acceleration:
        for ii = 1:length(ts)
            
            sn = listing(ii,1);
            en = listing(ii,2);
            cn = listing(ii,3);
            tn = listing(ii,4);
            
            %initialize data:
%             data.subject(sn).exp(en).condition(cn).Xacc_nobaseline{tn} = [];
%             data.subject(sn).exp(en).condition(cn).RTInd(tn,1:2) = {NaN,NaN};
            

            [~,ind] = min(abs(ts{ii}-tx)); %time = 100 ms
            a = xs{ii}(ind+[1:win]);
            b = heads{ii}(ind+[1:win]);
            c = xaccs{ii}(ind+[1:win]);
            
            %normalize this data the same as baseline data:
            a = (a - norm_factors(2,1))/(norm_factors(1,1)-norm_factors(2,1));
            b = (b - norm_factors(2,2))/(norm_factors(1,2)-norm_factors(2,2));
            c = (c - norm_factors(2,3))/(norm_factors(1,3)-norm_factors(2,3));
            
            %calculate "match" scores:
            scores = nan(ntrials,1);
            for jj = 1:ntrials
%                 scores(jj) = sum(sqrt((1.*(a-factors{jj,1})).^2 + 1.*(b-factors{jj,2}).^2 + 1*(c-factors{jj,3}).^2)); %This caluclation incorporates all factors (x position, heading, and x acc) 
                scores(jj) = sum(sqrt((mean(c)-mean(factors{jj,3})).^2));
            end
            
            %sort scores in ascending order:
            [scores_sorted, sort_ind] = sort(scores);
            scores_keep = zeros(size(scores_sorted));
            
            
            %Two ways to identify best matches:
            %Method 1 (default): keep top 10 best scores:
            if cn == 1
                scores_keep(2:11) = 1;  %skip the first one, as that is the trial in question (don't want to use the same trial for the average)
            else
                scores_keep(1:10) = 1;
            end
            %Method 2: keep top scores under a certain threshold: (not
            %currently used)
%             scores_keep = scores_sorted < thresh;
            
            
            %If there are less than 3 similar trials, skip this trial and exclude it:
            if sum(scores_keep) < 3
                disp(['Not enough matching baseline S',num2str(sn), ' E', num2str(en), ' C', num2str(cn), ' T', num2str(tn)])
                data.subject(sn).exp(en).condition(cn).IsIncluded{tn} = false;
                data.subject(sn).exp(en).condition(cn).RTInd(tn,:) = {[];[]};
                continue
            elseif sum(scores_keep) > 30
                %if there are more than 30 similar trials, take only the 30
                %best matching trials: (this is only for METHOD 2 above,
                %and will do nothing in Method 1, as by definition the top
                %10 matches are included)
                scores_keep = scores_keep(1:30);
            end
            
            
            avg_data = nan(sum(scores_keep),max(lengths),4);
            
            ys_b = pullData(data, sn, en_b, 1, 'included','Right_HandY');
            
            for jj = 1:sum(scores_keep)
                
                thisind = sort_ind(jj);
%                 plot(ax1,xs_b{thisind},ys_b{thisind},'k');
                
                
%                 plot(ax2,ts_b{thisind},xs_b{thisind},'k');
                temp = xs_b{thisind}(start_inds(thisind):end);  %data for this baseline trial to include in avg
                avg_data(jj,1:length(temp),2) = temp;
                
%                 plot(ax3,ts_b{thisind},xacc_b{thisind},'k');
                temp = xacc_bns{thisind}(start_inds(thisind):end);  %data for this baseline trial to include in avg
                avg_data(jj,1:length(temp),3) = temp;
                
%                 plot(ax4,ts_b{thisind},diff_heads_b{thisind},'k');
                temp = diff_heads_b{thisind}(start_inds(thisind):end);  %data for this baseline trial to include in avg
                avg_data(jj,1:length(temp),4) = temp;
                
%                 plot(ax5, ts_b{thisind}, ys_b{thisind},'k');
                
                
            end
            
            avg_data_old = avg_data;
            std_data = nanstd(avg_data);
            avg_data =  nanmean(avg_data);
            avg_data = reshape(avg_data,max(lengths),4);
            std_data = reshape(std_data,max(lengths),4);
            
            %trim data to same size:
            lens = [length(avg_data(:,3)),length(xaccs{ii}(ind-tx*1000+1:end))];
            enpt = min(lens);
            
            %subtract the average acceleration from this trial's
            %acceleration:
            newacc = xaccsns{ii};
            newacc(ind-tx*1000-100+[1:enpt]) = newacc(ind-tx*1000-100+[1:enpt]) - avg_data(1:enpt,3);
            
            if cn ~= 1
            %also demean the trial acceleration in this window:
            newacc = newacc - mean(newacc(ind+[0:win]));
            end
            
            %save this new acceleration:
            data.subject(sn).exp(en).condition(cn).Xacc_nobaseline{tn} = newacc;
            
           
            outcount = outcount+1;
            

        end
        
    end 
end

disp('');
disp('Processing Complete');


