function [bounds_out,classes_out,signal_names_out] = croptointervalsTsonga(bounds,classes,signal_names)
%function [bounds_out,classes_out,signal_names_out] = croptointervalsTsonga(bounds,classes,signal_names)
%
% Limit algorithm output to evaluation intervals used in the zerospeech
% challenge.


% Load evalution intervals

[signname_eva,onset_eva,offset_eva] = textread('xitsonga.split','%s %f %f','delimiter',' ');

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
    
    
     
     % Get current boundaries 
     curbounds = bounds(k,:);
     % Find the first interval onset that is bigger than the current onset
     refbounds = B{sign_numbers(k)};
     
     
     if(curbounds(2) < refbounds(1))
         toremove(k) = 1;
     elseif(curbounds(1) > refbounds(2))
         toremove(k) = 1;
     elseif(curbounds(1) < refbounds(1))
         curbounds(1) = refbounds(1);
     elseif(curbounds(2) > refbounds(2))
         curbounds(2) = refbounds(2);
     end
     
     bounds(k,:) = curbounds;

end

toremove = boolean(toremove);

bounds_out = bounds;
classes_out = classes;
signal_names_out = signal_names;

bounds_out(toremove,:) = [];
classes_out(toremove,:) = [];
signal_names_out(toremove) = [];
