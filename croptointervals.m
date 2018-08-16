function [bounds_out,classes_out,signal_names_out] = croptointervals(bounds,classes,signal_names)
%function [bounds_out,classes_out,signal_names_out] = croptointervals(bounds,classes,signal_names)
%
% Limit algorithm output to evaluation intervals used in the zerospeech
% challenge.

% Load evalution intervals

[signname_eva,onset_eva,offset_eva] = textread('english.split','%s %f %f','delimiter',' ');



all_signalnames = unique(signname_eva);
% Convert signal names to signal integers

if(iscell(signal_names(1)))
    sign_numbers = zeros(length(signal_names),1);
    for k = 1:length(signal_names)
        sign_numbers(k) = find(strcmp(all_signalnames,signal_names{k}));
    end
else
    sign_numbers = signal_names;
end

% Gather split information into arrays specific to each signal to speed up
% computations

all_bounds = [onset_eva offset_eva];

B = cell(length(unique(signname_eva)),1);
for k = 1:length(unique(signname_eva))
    B{k} = all_bounds(strcmp(signname_eva,all_signalnames{k}),:);        
end

toremove = zeros(length(bounds),1);
for k = 1:length(bounds)   % Go through each boundary
    
     d1 = bounds(k,1)-B{sign_numbers(k)}(:,1);
     
     % Get current boundaries 
     curbounds = bounds(k,:);
     % Find the first interval onset that is bigger than the current onset
     onset_loc = find(d1 < 0,1);     
     if(isempty(onset_loc))
         onset_loc = 1;
     end
    
     % Select potential intervals for the given segment
     intervals = B{sign_numbers(k)}(max(onset_loc-1,1):min(onset_loc+1,length(B{sign_numbers(k)})),:);
     
     % Find which interval is most covered with the current segment
     dur = zeros(size(intervals,1),1);
     for ii = 1:size(intervals,1)   
         len = intervals(ii,2)-intervals(ii,1); % Length of the interval
         if(curbounds(1) < intervals(ii,1) && curbounds(2) > intervals(ii,2))  % Onset and offset outside interval
            dur(ii) = len;
        elseif(curbounds(1) < intervals(ii,1) && curbounds(2) < intervals(ii,2)) % Offset in interval
            dur(ii) = curbounds(2)-intervals(ii,1);
         elseif(curbounds(1) > intervals(ii,1) && curbounds(2) > intervals(ii,2)) % Onset in interval
            dur(ii) = intervals(ii,2)-curbounds(1);
        elseif(curbounds(1) > intervals(ii,1) && curbounds(2) < intervals(ii,2))  % Onset and offset within interval
            dur(ii) = curbounds(2) - curbounds(1);
         else % Not within interval at all
             dur(ii) = 0;
        end        
     end
     
     % Pick the interval where most of the segment is
     [val,bestii] = max(dur);
     
     if(val > 0) % If has finite duration within some interval
         bounds(k,1) = max(bounds(k,1),intervals(bestii,1));
         bounds(k,2) = min(bounds(k,2),intervals(bestii,2));     
     else % Remove the segment if it's not within any evaluation interval
         toremove(k) = 1;
     end
end

toremove = boolean(toremove);

bounds_out = bounds;
classes_out = classes;
signal_names_out = signal_names;

bounds_out(toremove,:) = [];
classes_out(toremove,:) = [];
signal_names_out(toremove) = [];
