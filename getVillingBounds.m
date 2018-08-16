function [bounds,bounds_t] = getVillingBounds(x,fs)
% function [bounds,bounds_t] = getVillingBounds(x)
%
% Segments audio waveforms in cell-array x into syllables using the VSeg
% algorithm by Villing, Ward, Timoney & Costello (2004). 
if nargin <2
    fs = 16000;
end


if(fs ~= 16000)    
    for k = 1:length(x)
    x{k} = resample(x{k},16000,fs);
    end
    fs = 16000;
end

bounds = cell(length(x),1);
bounds_t = cell(length(x),1);
for signal = 1:length(x)
    
    % Scale envelopes between 0-1, limiting recording artifacts to 4 SDs
    % above the mean amplitude.
    tmp = x{signal};
    tmp(tmp > mean(tmp)+4.*std(tmp)) = mean(tmp)+4.*std(tmp);
    tmp(tmp < mean(tmp)-4.*std(tmp)) = mean(tmp)-4.*std(tmp);
    tmp = tmp./max(abs(tmp));
    
    % Run actual VSegh
    tmp = villingSeg(tmp);
    bounds{signal} = tmp.*1000;
    bounds_t{signal} = tmp;    
end