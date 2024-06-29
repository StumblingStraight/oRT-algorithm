function RTs = individualRegression(accs, ts, varargin)
%This function runs the individual regression oRT detection algorithm. This
%algorithm identifies an oRT by a two step process. First, the last point
%before the peak of the acceleration is identified where the acceleration
%crosses above a certain threshold. By default, this threshold is set at
%20. Next, a regression on this point and the preceding 9 points is
%performed draw a line on the rising slope of the acceleration. The point
%where this line crosses zero is defined as the oRT for that trial.
%
%Inputs:
%   - accs: an mx1 cell array of acceleration traces
%   - ts:   an mx1 cell array of time data series
%   - varargin:
%       1: the first argument passed in after ts is the acceleration
%       threshold. If no argument is passed in, the default is 20.
%       2: the second argument passed in is the cutoff frequency for a
%       low-pass filter. The default is 100, which will not run any
%       filtering
%
%Outputs:  
%   - RTs: an mx1 array of the oRT times detected
%
%Author: D. Tanis PhD, 6/2024

%set defaults:
fc = 100;
xthresh = 20;

%interpret varargin inputs:
if ~isempty(varargin)
    if length(varargin) == 1
        xthresh = varargin{1};
    elseif length(varargin) >= 2
        xthresh = varargin{2};
        fc = varargin{2};
    end
end

%setup filtering:
fs = 1000;
n_poles = 3;
[b,a] = butter(n_poles,fc/(fs/2),'low');

RTs = nan(length(accs),1);

filtaccs = accs;

%step through each trial:
for ii = 1:length(accs)
    
    %identify a window to identify the oRT as a window from the cue to +500
    %ms
    [~,st] = min(abs(ts{ii}));
    [~,en] = min(abs(ts{ii}-.5));
    
    %set any nan values to zero:
    accs{ii}(isnan(accs{ii})) = 0;
    
    %do filtering if indicated:
    if fc < 100
        thisacc = filtfilt(b,a,accs{ii});
        filtaccs{ii} = thisacc;
    else
        thisacc = accs{ii};
    end
    
    thist = ts{ii};

    %Window in which to do the regression:
    win = 20;
    
    %find the peak of the acceleration first:
    [~,pkind] = max(thisacc(st+50:en+100));
    pkind = pkind+st+50-1;
    if pkind > en
        pkind = en;
    end
    
    %Find the last point before the peak that acceleration is less than the
    %threshold:
    ind = find(thisacc(st:pkind)<xthresh,1,'last') + st -1+1;
    
    %perform a linear regression on the acceleration in a window around the
    %threshold crossing:
    p = polyfit(thist(ind+[-win:win]), thisacc(ind+[-win:win]),1);
    newy = thist.*p(1) + p(2);
    
    %The RT is the point that this regression line crosses zero:
    [~,RTind] = min(abs(newy));
    

    if ~isnan(RTind)
        %Convert the RT index to the actual time of the RT:
        RTs(ii) = thist(RTind);
    end
    
end

end