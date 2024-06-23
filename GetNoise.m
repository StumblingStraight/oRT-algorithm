%GetNoise
%
%This script curates a set a noise profiles for each subject. These noise 
%traces can then be added to the subject's canonical correction to generate
%hybrid trials.
%
%Inputs:
%   - this script assumes the data is already loaded into the base
%   workspace labeled "data"
%   - this script assumes that the canonical corrections have already been
%   generated and are in the base workspace (see GenerateCCs.m)
%
%Outputs:
%   - this script outputs a noise_basis cell array and a noise_basis_t cell
%   array to work with the cell arrays of the canonical corrections
%
%Author: D. Tanis PhD, 6/2024

subjs = 1:length(acc_all);

noise_basis = cell(length(acc_all),1);
noise_basis_t = noise_basis;

total_included = 0;

%step through each subject:
for sn = subjs
    
    %pull baseline data for this subject:
    t = pullData(data, sn, [1:2],1,'Included','Time_rel_Cue');
    x = pullData(data, sn, [1:2],1,'Included','Right_HandX');
    y = pullData(data, sn, [1:2],1,'Included','Right_HandY');
%     acc = pullData(data, sn, [1:2],1,'Included','Right_HandX_Acc');
    acc = pullData(data, sn, [1:2],1,'Included','Xacc_nobaseline');
    
    %check each trial for inclusion. Inclusion criteria are:
    % x within -0.5 and + 0.5
    % xacc < 50% of the peak of the canonical acc
    keepthese = false(size(x));
    canon_pk = max(acc_all{sn});
    count = 1;
    noise_basis{sn} = cell(1);
    
    %step through each trial:
    for ii = 1:length(x)
        
        %access this trial:
        thisacc1 = acc{ii}(1:end-100);
        
        %expand the trial in case the canonical correction is shorter than
        %this trial:
        thisacc = [-flipud(thisacc1); thisacc1; -flipud(thisacc1)];
        
        %align the cue-times of the noise trace and the CC trace and pull
        %the noise trace with that alignment:
        [~,t0] = min(abs(t{ii}));
        lowind = round(min(t_all{sn})*1000);

        access_inds = lowind+t0+[0:length(acc_all{sn})-1]+length(thisacc1);
        
        newnoise = thisacc(access_inds);
        
        %high-pass filter the noise data to remove any overt corrections
        %that might be present:
        fc = 1;
        fs = 1000;
        n_poles = 3;
        [b,a] = butter(n_poles,fc/(fs/2),'high');
        newnoise = filtfilt(b,a,newnoise);
        
        %generate the X position:
        new_x = cumsum(cumsum(newnoise))/1000000;
        if any(abs(new_x(1:1000))>0.5)
            %if this meets exclusion criteria, don't include it
            continue
        end
        
        if any(abs(newnoise) > 0.9*canon_pk)
            %if this meets exclusion criteria, dont include it:
            continue
        end

        %save this noise trace:
        keepthese(ii) = true;
        
        noise_basis{sn}{count} = newnoise;
        noise_basis{sn}{count+1} = -newnoise;   %keep the trial's inverse as well
        count = count + 2;
        

    end
    
    disp([num2str(sum(keepthese)), ' baseline noise trials found for subject ', num2str(sn)]);
    
    total_included = total_included + sum(keepthese);
    
end

function accs2 = indShift(accs, indshift)
%This function shifts the given data by the given number of indices from
%the indshift arrray. Empty data is padded with zeros. Positive shifts
%will shift towards later times.

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