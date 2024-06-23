%GenerateCCs
%
%This script uses the jump trials of the given data to generate canonical
%corrections (CCs).
%Inputs:
%   - this script assumes the data is already loaded into the base
%   workspace labeled "data"
%
%Outputs:
%   - this script outputs the data for the CCs in the following data
%   arrays or cells: acc_all, t_all, x_all, y_all, list_all, rt_all
%   These are all 17x1 cells or arrays (with the exception of the listing)
%   - this script will also disply the total number of trials used to
%   generate all canonicals
%
%Author: D. Tanis PhD, 6/2024


sns = [1:17];   %subjects to include

datatype = 'Xacc_nobaseline';   %data channel to use

%initialize data structures
acc_all = [];
t_all = [];
x_all = [];
y_all = [];
list_all = [];
rt_all = [];
vel_offsets = zeros(17,1);

acc_l = [];
vel_l = [];
acc_nb_l = [];
x_l = [];
y_l = [];
t_l = [];


totaltrials = 0;

%step through each subject
for sn = sns

acc = [];
t = [];
vel = [];
acc_nb = [];
x = [];
y = [];


cond = 4;   %condition 4 is 3cm jumps to the right
ex = [1,2]; %use both experimental conditions

[acc, listing] = pullData(data, sn, ex,cond,'Included',datatype);
vel = pullData(data, sn, ex,cond,'Included','Right_HandX_Vel');
x = pullData(data, sn, ex,cond,'Included','Right_HandX');
y = pullData(data, sn, ex,cond,'Included','Right_HandY');
acc_nb = pullData(data, sn, ex,cond,'Included','Xacc_nobaseline');
% acc_l = pullData(data, sn, ex,5,'Included',datatype);
t = pullData(data, sn, ex,cond,'Included','Time_rel_Cue');

cond = 5;   %condition 5 is 3cm jumps to the left
[acc_l, listing_l] = pullData(data, sn, ex,cond,'Included',datatype);
vel_l = pullData(data, sn, ex,cond,'Included','Right_HandX_Vel');
x_l = pullData(data, sn, ex,cond,'Included','Right_HandX');
y_l = pullData(data, sn, ex,cond,'Included','Right_HandY');
acc_nb_l = pullData(data, sn, ex,cond,'Included','Xacc_nobaseline');
t_l = pullData(data, sn, ex,cond,'Included','Time_rel_Cue');

%invert the left-ward traces:
for ii = 1:length(acc_l)
    acc_l{ii} = acc_l{ii}*-1;
    acc_nb_l{ii} = acc_nb_l{ii}*-1;
    vel_l{ii} = vel_l{ii}*-1;
    x_l{ii} = x_l{ii}*-1;
end

listing = [listing; listing_l];

acc = [acc; acc_l];
vel = [vel;vel_l];
acc_nb = [acc_nb; acc_nb_l];
x = [x;x_l];
y = [y;y_l];
t = [t;t_l];

%run the grandmeanRT script to get the grand mean of the data as well as
%the regressed oRT:
[gmRT, mean_acc, mean_t] = grandmeanRT(acc,t);

%generate the X position trace:
mean_x = cumsum(cumsum(mean_acc))/1000000;
[mean_y,~] = grandMean(y,t);

%get the velocity offsets:
vel = cumsum(mean_acc)/1000;
[~,ind0] = min(abs(mean_t));
vel_offsets(ii) = vel(ind0);

if ~isnan(gmRT)
    
    totaltrials = totaltrials + length(acc);
    acc_all = [acc_all; {mean_acc'}];
    t_all = [t_all; {mean_t}];
    x_all = [x_all; {mean_x'}];
    y_all = [y_all; {mean_y}];
    list_all = [list_all; listing];
    rt_all = [rt_all; gmRT];
end

end

disp(['Total number of trials used is ', num2str(totaltrials)]);



function [mean_data, mean_t] = grandMean(data, ts)
%This function generates a grand mean of the given data, time-aligned to
%the cue

starttime = -0.5;
entime = .8;
data_clip = nan(length(data), round((entime-starttime)*1000+1));
for ii = 1:length(data)
    [~,st] = min(abs(ts{ii}-starttime));
    [~,en] = min(abs(ts{ii}-entime));
    temp = data{ii}(st:en);
    data_clip(ii,1:length(temp)) = temp;
end

mean_data = nanmean(data_clip);
mean_data = [mean_data, zeros(1,800)];
mean_t = starttime:0.001:entime+0.8;
end

function [RT, mean_acc, mean_t] = grandmeanRT(accs, ts)
%This function generates a grand mean of the given trials and then performs
%a regression to identify the oRT

try
    
%time-algin all trials to the cue and average:
starttime = -0.5;
entime = .8;
acc_clip = nan(length(accs), round((entime-starttime)*1000+1));
for ii = 1:length(accs)
    [~,st] = min(abs(ts{ii}-starttime));
    [~,en] = min(abs(ts{ii}-entime));
    temp = accs{ii}(st:en);
    acc_clip(ii,1:length(temp)) = temp;
end

mean_acc = nanmean(acc_clip);
mean_acc = [mean_acc, zeros(1,800)];
mean_t = starttime:0.001:entime+0.8;

xp = 0.25;  %percentage of peak
win = 10;   %number of samples in the window to regress
mint = 0.05;
maxt = 0.5;

%find the regressed RTind. If the value is too low, repeat with a higher
%threshold.
while true
    
%determine the peak and the threshold values:
[peakval,pkind] = max(mean_acc(600+[0:400])); %look for the peak between +100 and +400 ms of the cue
pkind = pkind+600-1;
thresh = xp*peakval;

ind = find(mean_acc(600:pkind)<thresh,1,'last') + 600-1;    %find the threshold crossing

if isempty(ind) && xp == 0.25
    %if no threshold crossing was identified, repeat with a higher
    %threshold
    xp = 0.4;
    continue
elseif isempty(ind) && xp == 0.4
    %if there is still no trheshold crossing identified within these
    %limits, return NaN
    RTind = NaN;
    break
end

%fit a line to the window preceding the threshold crossing:
p = polyfit(mean_t(ind+[-win:0]), mean_acc(ind+[-win:0]),1);
newy = mean_t.*p(1) + p(2);

%The oRT is the timepoint that this line crosses zero:
[~,RTind] = min(abs(newy));


if xp == 0.25 && (mean_t(RTind) < mint || mean_t(RTind) > maxt)
    %if the regression failed to find a value within the limits, repeat with a
    %higher threshold crossing value:
    xp = 0.4;
    continue;
elseif xp == 0.4 && (mean_t(RTind) < mint || mean_t(RTind) > maxt)
    %if the regression failed with the higher value, return NaN:
    RTind = NaN;
    break
end

%if an RT was successfully found, break out of the while loop:
if mean_t(RTind) >= mint
    break
end
end

%find the time of the RT, given its index:
if ~isnan(RTind)
    RT = mean_t(RTind);
else
    RT = nan;
end

catch ME
    ab = 2; 
end

end
