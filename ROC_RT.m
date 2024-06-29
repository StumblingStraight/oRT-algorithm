function [meanRT] = ROC_RT(accs,ts,nonperturb_accs)
%This function runs the ROC oRT detection algorithm. This algorithm
%identifies oRTs using an ROC curve, comparing the perturbed trials with
%non-preturbed trials. This algorithm fits a regression line to the ROC
%trace around the first point that crosses above 75% discrimination. The
%point on this regressed line that crosses 50% discrimination is identified
%as the oRT for that group of trials.
%
%Note that for purposes of oRT algorithm testing, velocity data was given
%as the input to this algorithm, as most ROC methods utilize velocity.
%
%Inputs:
%   - accs: an mx1 cell array of acceleration traces
%   - ts:   an mx1 cell array of time data series
%   - nonperturb_accs: an mx1 cell array of acceleration traces from trials
%   with no perturbation (AKA baseline trials)
%
%Outputs:  
%   - meanRT: a scalar value indicating the mean RT for this batch of
%   trials
%
%Author: D. Tanis PhD, 6/2024

%threshold:
thresh = 0.75;

%window around the threshold crossing to do the regression:
win = 20;


%first, clip all the traces to be the same length, aligned to the cue:
allb = nan(10000,length(accs));
allp = allb;
minlen = 10000;
cues = nan(length(accs),1);
for ii = 1:length(accs)
    
    [~,cue] = min(abs(ts{ii}));
    cues(ii) = cue;
    perturb = accs{ii}(cue:end);
    baseline = nonperturb_accs{ii}(cue:end);
    allp(1:length(perturb),ii) = perturb;
    allb(1:length(baseline),ii) = baseline;
    
    if length(perturb) < minlen
        minlen = length(perturb);
    end
end

allp = allp(1:minlen,:);
allb = allb(1:minlen,:);

%run the ROC:
[rocout, ~,~,~] = running_roc(allb', allp');
roct = 1:length(rocout);

%find the first point that the ROC crosses above 95%
ind95 = find(rocout>0.95, 1,'first');

%Now backtrack to the most recent time that the ROC crossed the threshold:
crossing = find(rocout(1:ind95) < thresh, 1, 'last');

%If the ROC never achieved 95% discrimination, leave this batch as NaN:
if isempty(crossing) || isnan(crossing)
    meanRT = NaN;
    return
end

%If the crossing detection is too early, warn the user and bail out:
if crossing < win+1
   warning("ROC detected a crossing too early for a valid regression window"); 
   meanRT = NaN;
   return
end

%perform the regression:
newt = [-1000:1000];
p = polyfit(roct(crossing+[-win:win]), rocout(crossing+[-win:win]),1);
newy = newt.*p(1) + p(2); 

%The RT is defined as the point in the regressed line that crosses 50%:
[~,RTind] = min(abs(newy-0.5));

RTind = RTind-1001; %correction for the 1000 samples before the cue time

%sometimes the regressed line is too shallow and points to a non-sensical
%time before the cue. To prevent this, a 1000ms window is applied around
%the threshold crossing, and the RT must fall within this window. If it
%doesn't, warn the user and bail out:
if RTind == 1 || RTind == length(newt)
   warning("ROC regression failed to detect a 50% crossing within the designated +/- 1000ms window");
   meanRT = NaN;
   return
end

%convert the RT index to the time:
meanRT = RTind/1000;

end

function [temp outt tp fp] = running_roc(x1,x2,conf);
% function [temp outt tp fp] = running_roc(x1,x2,conf);
% This function performs ROC analysis on a series of data and returns the
% descriminability of one observation over the other.
% Inputs: x1 = matrix of repeated observations, condition 1
%         x2 = as above, condition 2
%         conf{1} = how many sig points in a row specify distinguishibility
%                   (default = 5)
%         conf(2) = significance level (default = 0.75)
%                    provide positive number if > conf(2) or negative 
%                    number if < conf(2).
% Outputs: Confidence with which an ideal observer can descriminate
%           the two signals as different, with direction information

if(nargin <= 2)
    conf = [5,0.75,1];
end
if(length(conf<3))
    conf(3) = 1;
end

for(q=1:max(size(x1(1,:))))
            obs = [(x1(:,q)); (x2(:,q))];
            tru = [1*ones(max(size(x1(:,q))),1); -1*ones(max(size(x2(:,q))),1)];
            [tp(q,:), fp(q,:)] = roc(tru,obs);
            temp(q) = sum([0 diff(tp(q,:))].*fp(q,:));
end

if(conf(2)>0)
    for(i = conf(3):max(size(x1(1,:)))-conf(1))
        if(sum(temp(i:i+(conf(1)-1)) > conf(2)) == conf(1))
            outt = i;
            return
        end
    end
else
   conf(2) = abs(conf(2));
   for(i = conf(3):max(size(x1(1,:)))-conf(1))
        if(sum(temp(i:i+(conf(1)-1)) < conf(2)) == conf(1))
            outt = i;
            return
        end
   end
end
outt = [-1];


end

function [tp, fp] = roc(t, y)
%
% ROC - generate a receiver operating characteristic curve
%
%    [TP,FP] = ROC(T,Y) gives the true-positive rate (TP) and false positive
%    rate (FP), where Y is a column vector giving the score assigned to each
%    pattern and T indicates the true class (a value above zero represents
%    the positive class and anything else represents the negative class).  To
%    plot the ROC curve,
%
%       PLOT(FP,TP);
%       XLABEL('FALSE POSITIVE RATE');
%       YLABEL('TRUE POSITIVE RATE');
%       TITLE('RECEIVER OPERATING CHARACTERISTIC (ROC)');
%
%    See [1] for further information.
%
%    [1] Fawcett, T., "ROC graphs : Notes and practical
%        considerations for researchers", Technical report, HP
%        Laboratories, MS 1143, 1501 Page Mill Road, Palo Alto
%        CA 94304, USA, April 2004.
%
%    See also : ROCCH, AUROC

%
% File        : roc.m
%
% Date        : Friday 9th June 2005
%
% Author      : Dr Gavin C. Cawley
%
% Description : Generate an ROC curve for a two-class classifier.
%
% References  : [1] Fawcett, T., "ROC graphs : Notes and practical
%                   considerations for researchers", Technical report, HP
%                   Laboratories, MS 1143, 1501 Page Mill Road, Palo Alto
%                   CA 94304, USA, April 2004.
%
% History     : 10/11/2004 - v1.00
%               09/06/2005 - v1.10 - minor recoding
%
% Copyright   : (c) G. C. Cawley, June 2005.
%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
%

% process targets

t = t > 0;

% sort by classifier output

[Y,idx] = sort(-y);
t       = t(idx);

% compute true positive and false positive rates

tp = cumsum(t)/sum(t);
fp = cumsum(~t)/sum(~t);

% add trivial end-points

tp = [0 ; tp ; 1];
fp = [0 ; fp ; 1];

end