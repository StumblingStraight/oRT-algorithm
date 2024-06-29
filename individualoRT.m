function RTs = individualoRT(accs, ts, varargin)
%This function runs the individual threshold oRT detection algorithm. This
%algorithm identifies an oRT in an acceleration trace as the last point
%before the peak in which the acceleration crosses above a ceratin
%threshold. By default, this threshold is 15.
%
%Inputs:
%   - accs: an mx1 cell array of acceleration traces
%   - ts:   an mx1 cell array of time data series
%   - varargin:
%       1: the first argument passed in after ts is the acceleration
%       threshold. If no argument is passed in, the default is 15.
%       2: the second argument passed in is the cutoff frequency for a
%       low-pass filter. The default is 100, which will not run any
%       filtering
%
%Outputs:  
%   - RTs: an mx1 array of the oRT times detected
%
%Author: D. Tanis PhD, 6/2024

%set defaults:
thresh_a = 15;
fc = 100;

%interpret varargin inputs:
if ~isempty(varargin)
    thresh_a = varargin{1};
    if length(varargin) > 1
        fc = varargin{2};
    end
end
  
%setup filtering:
fs = 1000;
n_poles = 3;
[b,a] = butter(n_poles,fc/(fs/2),'low');


RTs = nan(length(accs),1);

%step through each acceleration trace:
for ii = 1:length(accs)
    
    %if an fc value is set, run filtering:
    if fc < 100
        filtacc = filtfilt(b,a,accs{ii});
    else
        filtacc = accs{ii};
    end
    
    %identify the index that corresponds to the cue:
    [~,cueind] = min(abs(ts{ii}));
    
    %first find the peak:
    [~,peakind] = max(filtacc(cueind+[50:600]));
    peakind = peakind+cueind+50-1;
    if peakind > cueind+500
        peakind = cueind+500;
    end
    
    %backtrack to the most recent point where acceleration crossed its
    %threshold:
    ai = find(filtacc(cueind:peakind)<thresh_a, 1, 'last')+cueind;

    rt = ai;
    
    if ~isempty(rt)
        %convert the identified RT index into the actual time of the RT:
        RTs(ii) = ts{ii}(rt);
    else
        %if no threshold crossing is identified, leave the RT as NaN
    end
    
end

end