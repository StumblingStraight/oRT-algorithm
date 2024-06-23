function varargout = pullData(data, subjN, expN, condN, trialN, dataType, varargin)
%This function pulls desired data from the given data struct for use with
%data analysis scripts.
%
%Inputs:
%   - data: the struct containing the experiment data
%   - subjN: the subject whose data should be pulled. Can be specified as a
%   scalar or a range (ie [1:5])
%   - expN: the experiment(s) from which to pull data. Can also be
%   specified as a scalar or a range (ie [1:2])
%   - condN: the condition(s) from which to pull data. Can be specified as
%   a scalar, a range (ie [1:5]), or with the text 'All' for all conditions
%   or with the text 'Perturb' for all conditions except the first (which
%   is generally the baseline condition)
%   - trialN: the trial(s) from which to pull data. Can be specified as a
%   scalar, a rand (ie [1:10]), or with the text 'All' for all trials, or
%   with the text 'Included' for all trials for which the IsIncluded field
%   is 1.
%   - dataType: the data channel to pull. This can be specified as the text
%   of the field name within the data struct (such as 'Right_HandX' for X
%   position of the right hand). There are also a number of derived
%   channels that can be requested, such as 'Right_HandX_Acc', which is the
%   double derivative of 'Right_HandX'. The following is a list of all
%   available channels:
%       ShoAng, ShoAcc, ElbAng, ElbAng_wrtExtension, OpticalTrigger,
%       PostDelt, PectMaj, Brachio, TriLat, Right_HandX, Right_HandX_Vel,
%       Right_HandX_Acc, Right_HandX_Acc_smoothed, Right_HandX_Jerk,
%       Xacc_nobaseline_smoothed, Xacc_nobaseline, Derived_X, Right_HandY,
%       Right_HandY_Vel, Right_HandY_Acc, Right_HandY_Acc_smoothed,
%       Heading, Change in Heading, Right_ShoTorIM, Right_ElbTorIM,
%       Time, Time_rel_Cue
%   Further information about each channel can be seen in its corresponding
%   section below.
%   - varargin: currently unused
%
%Outputs
%   The requested data is passed out as the first output, as in
%   [output] = pullData(inputs);
%   If multiple data series are requested (such as multiple trials or
%   channels), the output will be as follows:
%   [outputdata, listing] = pullData(inputs), where:
%       outputdata is an m x 1 cell array, where m corresponds to all the
%       available trials of data
%       listing is an m x 4 array, in which the columns correspond to subject
%       number, experiment number, condition number, and trial number for
%       the given trial in question.
%
%Author: D. Tanis PhD, 6/2024
%


if any([length(subjN), length(expN), length(condN), length(trialN)] > 1)
    [outdata, outheader] = pullDataMulti(data, subjN, expN, condN, trialN, dataType, varargin);
    varargout = {outdata, outheader};
    return
end

samprate = 1000;

if strcmp(dataType,'Right_HandX_Vel')
    %if requested data is right hand X velocity, pull the position data and
    %differentiate to get the velocity
    
    thisdata = data.subject(subjN).exp(expN).condition(condN).('Right_HandX'){trialN};
    
    outdata = [0; diff(thisdata)];    %note that this is zero-padded at the beginning to yeild the same length array
    outdata = outdata*samprate;      %correct for samprate (ie velocity = change in X/ change in time)
    
elseif strcmp(dataType,'Right_HandX_Acc')
    %if requested data is right hand X acceleration, pull the position data and
    %differentiate twice to get the acceleration
    
    thisdata = data.subject(subjN).exp(expN).condition(condN).('Right_HandX'){trialN};
    
    outdata = [0; diff([0; diff(thisdata)])];    %note that this is zero-padded at the beginning to yeild the same length array
    outdata = outdata*samprate*samprate;      %correct for samprate (ie velocity = change in X/ change in time)
    
    %add NaNs to the end to eliminate spikes that sometimes show up in data
    outdata(1:5) = 0;
    outdata(end-5:end) = 0;
    
elseif strcmp(dataType, 'Right_HandX_Acc_smoothed')
    
    tempdata = pullData(data, subjN, expN, condN, trialN, 'Right_HandX_Acc');
    
    fc1 = 1;
    fc = 12; %cutoff frequency
    fs = samprate;
    n_poles = 6;
    
%     [b,a] = butter(1,[fc1/(fs/2)],'high');
%     tempdata = filtfilt(b,a,tempdata);
    
    [b,a] = butter(n_poles,[fc/(fs/2)],'low');
    outdata = filtfilt(b,a, tempdata);
    
    
elseif strcmp(dataType, 'Xacc_nobaseline_smoothed')
    
    tempdata = pullData(data, subjN, expN, condN, trialN, 'Xacc_nobaseline');
    
    fc1 = 0.5;
    fc = 5; %cutoff frequency
    fs = samprate;
    n_poles = 6;
    
    nanflag = 0;
    if any(isnan(tempdata))
        nanind = find(isnan(tempdata),1,'first');
        endind = length(tempdata);
        tempdata = tempdata(1:nanind-1);
        nanflag = 1;
    end
    
%     [b,a] = butter(1,[fc1/(fs/2)],'high');
%     tempdata = filtfilt(b,a,tempdata);
%     
    [b,a] = butter(n_poles,[fc/(fs/2)],'low');
    outdata = filtfilt(b,a, tempdata);
    
    %add nans back in
    if nanflag
    buffer = nan(endind-nanind+1,1);
    outdata = [outdata;buffer];
    end
    
    
    
elseif strcmp(dataType, 'Right_HandX_Jerk')
    
    tempdata = pullData(data, subjN, expN, condN, trialN, 'Right_HandX_Acc');
    
    fc = 8; %cutoff frequency
    fs = samprate;
    n_poles = 6;
    
    [b,a] = butter(n_poles,fc/(fs/2),'low');
    tempdata = filtfilt(b,a, tempdata);
    
    outdata = [0; diff(tempdata)];
    outdata = outdata*samprate;
    
elseif strcmp(dataType,'Right_HandY_Vel')
    %if requested data is right hand Y velocity, pull the position data and
    %diff to get the velocity
    
    thisdata = data.subject(subjN).exp(expN).condition(condN).('Right_HandY'){trialN};
    
    outdata = [0; diff(thisdata)];    %note that this is zero-padded at the beginning to yeild the same length array
    outdata = outdata*samprate;      %correct for samprate (ie velocity = change in X/ change in time)
    
elseif strcmp(dataType,'Right_HandY_Acc')
    %if requested data is right hand Y acceleration, pull the position data and
    %diff twice to get the acceleration
    
    thisdata = data.subject(subjN).exp(expN).condition(condN).('Right_HandY'){trialN};
    
    outdata = [0; diff([0; diff(thisdata)])];    %note that this is zero-padded at the beginning to yeild the same length array
    outdata = outdata*samprate*samprate;      %correct for samprate (ie velocity = change in X/ change in time)
    
elseif strcmp(dataType, 'Right_HandY_Acc_smoothed')
    
    tempdata = pullData(data, subjN, expN, condN, trialN, 'Right_HandY_Acc');
    fc1 = 1;
    fc = 7; %cutoff frequency
    fs = samprate;
    n_poles = 6;
    
    [b,a] = butter(n_poles,[fc1/(fs/2),fc/(fs/2)]);
    outdata = filtfilt(b,a, tempdata);

elseif strcmp(dataType, 'Heading')
    
    x = data.subject(subjN).exp(expN).condition(condN).('Right_HandX'){trialN};
    y = data.subject(subjN).exp(expN).condition(condN).('Right_HandY'){trialN};
    [~, st] = min(abs(y-2));
    [~, en] = min(abs(y-19));
%     en = length(y);
    
    ydiff = diff(y);
    xdiff = diff(x);
    ydiff2 = [0; diff(ydiff)];
    heading = [0; atand(diff(x)./ydiff)];  %this gives the heading of the cursor in degrees. 0 deg = straight ahead
    prevydir = 1;
    offset = 0;
    for ii = st:en-1
        %if the y direction is now negative, correct the heading:
        thisydir = ydiff(ii) > 0;
        if thisydir ~= prevydir
            if ydiff(ii) < 0
                if xdiff(ii) >0
                    offset = offset +180;
                else
                    offset = offset -180;
                end
            else
                if xdiff(ii) >0
                    offset = offset -180;
                else
                    offset = offset +180;
                end
            end
        end
        prevydir = thisydir;
        heading(ii+1) = heading(ii+1) + offset;
    end
    
    %get rid of data beofre movement to avoid weird discontinuities in
    %heading data
    heading(1:round(st)) = NaN;
    heading(round(en):end) = NaN;
    
    outdata = heading;
    
elseif strcmp(dataType, 'Change in Heading')
    
    thisdata = pullData(data, subjN, expN, condN, trialN, 'Heading');
    outdata = diff([NaN; thisdata]);
    outdata = outdata*samprate;      %correct for samprate (ie velocity = change in X/ change in time)

elseif strcmp(dataType, 'Change in Heading smoothed')
    
    thisdata = pullData(data, subjN, expN, condN, trialN, 'Change in Heading');
    
    %clip data to not have nans
    st = find(~isnan(thisdata),1);
    subdata = thisdata(~isnan(thisdata));
    
    fc = 8; %cutoff frequency
    fs = samprate;
    n_poles = 6;
    
    [b,a] = butter(n_poles,fc/(fs/2),'low');
    subdata = filtfilt(b,a, subdata);
    outdata = thisdata;
    outdata([1:length(subdata)]+st) = subdata;
    
elseif strcmp(dataType, 'TangentialVel')
    
    xv = pullData(data, subjN, expN, condN, trialN, 'Right_HandX_Vel');
    yv = pullData(data, subjN, expN, condN, trialN, 'Right_HandY_Vel');
    
    outdata = sqrt(xv.^2+yv.^2);
        
elseif strcmp(dataType, 'ElbAng_wrtExtension')
    
    thisdata = pullData(data, subjN, expN, condN, trialN, 'ElbAng');
    outdata = 180-thisdata;
    
elseif strcmp(dataType,'ShoAcc')
    %if requested data is Shoulder acceleration, pull the data and
    %diff twice to get the acceleration
    
    thisdata = data.subject(subjN).exp(expN).condition(condN).('ShoAng'){trialN};
    
    %angle data must be smoothed to yeild good data:
    fc = 20; %cutoff frequency
    fs = samprate;
    n_poles = 6;
    
    [b,a] = butter(n_poles,fc/(fs/2),'low');
    thisdata = filtfilt(b,a, thisdata);
    
    outdata = [0; diff([0; diff(thisdata)])];    %note that this is zero-padded at the beginning to yeild the same length array
    outdata = outdata*samprate*samprate;
    
elseif strcmp(dataType, 'Time')
    %requested data is a time vector:
    
    y = pullData(data, subjN, expN, condN, trialN, 'Right_HandY');
    t = [1:length(y)]/samprate;
    
    outdata = t;
    
elseif strcmp(dataType, 'Time_rel_Cue')
    %requested data is a time vector relative to cue presentation:
    %in the color experiment, Subj 9, exp 2, cond 1, trial 14 has a faulty
    %trigger
        
    y = pullData(data, subjN, expN, condN, trialN, 'Right_HandY');
    t = [1:length(y)]/samprate;
    
    
    trig = abs(pullData(data, subjN, expN, condN, trialN, 'OpticalTrigger'));
    
    stdev = std(trig(1:200));
    thresh = mean(trig(1:200))+2.5*stdev;
    thresh2 = mean(trig(1:200))+3*stdev;
    trig(1:200) = mean(trig(1:200));
    
    %normalize:
%     trig = trig-min(trig);
%     trig = trig/max(trig);
%     trigidx = find(trig>0.3, 1);
    
    crossing = trig>thresh2;
    trigidx = find([diff(crossing)] > 0, 1, 'last');
    if trigidx > 0.9*length(trig)
       %The old optical trigger does not stay above threshold for the rest of the time
       %In this case, take the first crossing above threshold:
       trigidx = find([diff(crossing)] > 0, 1, 'first');
    end
%     trigidx = find(trig>thresh2, 1,'last');
%     if isempty(trigidx)
%         trigidx = find(trig>thresh2, 1);
%     end
    
    %if the trigger index is too early or was not found, the optical trigger is not working.
    %Pull the event code instead:
    if isempty(trigidx) || trigidx < 100
        warning(['Faulty optical trigger. Subject ', num2str(subjN), ' Exp ', num2str(expN), ' Cond ', num2str(condN), ' trial ', num2str(trialN)]);
        before_samps = 200;
        delay_samps = 50;
        reach_cue = data.subject(subjN).exp(expN).condition(condN).events{trialN}.TIMES(2);
        zerotime = reach_cue-before_samps/samprate;
        trigtime = data.subject(subjN).exp(expN).condition(condN).events{trialN}.TIMES(4);
        trigtime = trigtime-zerotime;
        trigidx = round(trigtime*samprate) + delay_samps;
        a = 2;
    end
    
    outdata = t-t(trigidx);
    
    
    
elseif strcmp(dataType, 'Dev_in_Heading')
    
    x0 = 0;
    y0 = 20;
    
    x = pullData(data, subjN, expN, condN, trialN, 'Right_HandX');
    y = pullData(data, subjN, expN, condN, trialN, 'Right_HandY');
    heading = pullData(data, expN, subjN, condN, trialN, 'Heading');
    
    xdif = x0-x;
    ydif = y0-y;
    
    [~, st] = min(abs(y-1));
    [~, en] = min(abs(y-18));
    
    headingtarg = atand(xdif./ydif);

    %get rid of data before movement to avoid weird discontinuities in
    %heading data
    headingtarg(1:round(st)) = NaN;
    headingtarg(round(en):end) = NaN;
    
    devinheading = abs(heading-headingtarg);
    
    outdata = devinheading;
    
elseif strcmp(dataType, 'RTs')
    %pull the reaction times, relative to the cue
    %Note that this returns a cell of RTs
    
    t = pullData(data, subjN, expN, condN, trialN, 'Time_rel_Cue');
    rtinds = data.subject(subjN).exp(expN).condition(condN).RTInd(trialN,:);
    
    if ~isempty(rtinds{1}) && ~isnan(rtinds{1})
        if rtinds{1} > length(t)
            temp = ['S', num2str(subjN), ' E', num2str(expN), ' C', num2str(condN), ' T', num2str(trialN)];
            error(['RT index (', num2str(rtinds{1}), ') exceeds length of data (', num2str(length(t)), '), for ', temp]);
        end
    end
    if ~isempty(rtinds{2}) && ~isnan(rtinds{2})
        if rtinds{2} > length(t)
            temp = ['S', num2str(subjN), ' E', num2str(expN), ' C', num2str(condN), ' T', num2str(trialN)];
            error(['RT index (', num2str(rtinds{2}), ') exceeds length of data (', num2str(length(t)), '), for ', temp]);
        end
    end
    
    if ~isempty(rtinds{1}) && ~isnan(rtinds{1})
        RT1 = t(rtinds{1});
    else
        RT1 = [];
    end
    
    if ~isempty(rtinds{2}) && ~isnan(rtinds{2})
        RT2 = t(rtinds{2});
    else
        RT2 = [];
    end
    
    outdata = {RT1, RT2};
    
elseif strcmp(dataType, 'RTInd')
    
    outdata = data.subject(subjN).exp(expN).condition(condN).RTInd(trialN,:);
  
elseif strcmp(dataType, 'Direction')
    
    outdata = data.subject(subjN).exp(expN).condition(condN).(dataType)(trialN,:);
    
elseif strcmp(dataType, 'Derived_X')
    
    outdata = data.subject(subjN).exp(expN).condition(condN).('Xacc_nobaseline'){trialN};
%     outdata = outdata-mean(outdata);
    outdata = cumsum(cumsum(outdata))/1000000;
    
else
    %if none of the above conditions apply, pull the data requested.

    outdata = data.subject(subjN).exp(expN).condition(condN).(dataType){trialN};

end

varargout{1} = outdata;

end

function [outdata, outheader] = pullDataMulti(data, subjN, expN, condN, trialN, dataType, varargin)
%This function runs for loops to pull multiple chunks of data:

% try

%first do a loop to count up the total amount of data:
checkN = false;
if ischar(trialN)
    if strcmpi(trialN, 'all') || strcmpi(trialN, 'included')
        checkN = true;
    else
        error('Unknown input string for trial number. Acceptable strings are "all" and "included".')
    end
end


if checkN
    totaln = 0;
    for thisS = subjN
        for thisE = expN
            thiscondN = checkConds(thisS, thisE, condN); 
            for thisC = thiscondN
                if strcmpi(trialN, 'all')
                    theseTs = length(data.subject(thisS).exp(thisE).condition(thisC).filename);
                elseif strcmpi(trialN, 'included')
                    isincluded = cell2mat(pullData(data, thisS, thisE, thisC, 'all','IsIncluded'));
                    theseTs = length(find(isincluded == true));
                elseif ~ischar(trialN)
                    theseTs = length(trialN);
                end
                totaln = totaln+theseTs;
            end
        end
    end
else
    totaln = length(subjN)*length(expN)*length(condN)*length(trialN);
end

outdata = cell(totaln, 1);
outheader = zeros(totaln, 4);

count = 1;
for thisS = subjN
    for thisE = expN
        thiscondN = checkConds(thisS, thisE, condN); 
        for thisC = thiscondN
             
            if strcmpi(trialN, 'all')
                trialNs = 1:length(data.subject(thisS).exp(thisE).condition(thisC).filename);
            elseif strcmpi(trialN, 'included')
                isincluded = cell2mat(pullData(data, thisS, thisE, thisC, 'all','IsIncluded'));
                trialNs = (find(isincluded == true))';
            elseif ~ischar(trialN)
                trialNs = trialN;
                [m,n] = size(trialNs);
                if m > 1 && n == 1
                    trialNs = trialNs';
                end
            end
            
            for thisT = trialNs

%                 try
                    
                thisdata = pullData(data, thisS, thisE, thisC, thisT, dataType);
                if iscell(thisdata)
                    outdata(count,1:length(thisdata)) = thisdata;
                else
                    outdata{count,1} = thisdata;
                end
                outheader(count,:) = [thisS, thisE, thisC, thisT];
                count = count + 1;
                
%                 catch ME
%                    ab = 1; 
%                 end
            end
        end
    end
end

if ~isempty(varargin)
   if strcmpi(varargin{1},'clipped')
      %if this is true, clip all the data in the outdata to the same length and put it all in a cell
      
      for ii = 1:length(outheader(:,1))
         %step through each row
         %pull this data's t vector
         pullData();
          
         %find the t = 0
         
         %save the t
      end
      
      %find the lowest number of samples before and after t = 0
      
      %clip all data to that minimum n samples before and after t = 0
      
      %save the data as an array
   end
end

% catch ME
%     a = 1;
% end

    function result = checkConds(sn, en, cond_input)
        %this function checks the condition input and converts any string
        %inputs to numeric.

        if isnumeric(cond_input)
            result = cond_input;
        elseif ischar(cond_input)
            if strcmpi(cond_input, 'All')
                result = [1:length(data.subject(sn).exp(en).condition)];
            elseif strcmpi(cond_input, 'Perturb')
                result = [2:length(data.subject(sn).exp(en).condition)];
            else
                error(['Unknown input string ' cond_input, ' for condition number. Acceptable strings are "all" and "perturb".']);
            end
        else
            error('Unknown input data type for condition number. Acceptable types are string or numeric.');
        end
    end
end