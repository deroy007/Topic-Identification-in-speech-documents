function [bounds,bounds_t] = getEnvelopeBounds(envelopes,threshold)
% function [bounds,bounds_t] = getEnvelopeBounds(envelopes,threshold)
%
% Performs minima detection from envelopes (a cell-array) using a threshold
% for the minimum-difference between the amplitude of the current minima 
% and the preceding local maxima.


if nargin <2
threshold = 0.12;
end

n_signals = length(envelopes);

bounds = cell(n_signals,1);
bounds_t = cell(n_signals,1);


for signal = 1:n_signals
    [maxtab,bounds_forward] = peakdet(envelopes{signal},threshold);
    if(~isempty(bounds_forward))
    bounds_forward = bounds_forward(:,1);
    end
    [maxtab_bw,bounds_bw] = peakdet(flipud(envelopes{signal}),threshold);
    
    if(~isempty(bounds_bw))
    bounds_bw = length(envelopes{signal})-bounds_bw(:,1)+1;    
    bounds{signal} = union(bounds_forward(:,1),bounds_bw(:,1));
    end
    
    tmp = max(find(envelopes{signal} > 0.05,1)-10,1);
    tmp2 = min(length(envelopes{signal})-find(flipud(envelopes{signal}) > 0.05,1)+10,length(envelopes{signal}));
    
    bounds{signal} = [tmp;bounds{signal};tmp2];
    
    if(~isempty(bounds{signal}))        
    bounds_t{signal} = bounds{signal}./1000;
        
    
    % Find boundaries closer than 50 ms
    a = find(diff(bounds_t{signal}) < 0.05);
    
    
    torem = boolean(zeros(length(bounds_t{signal}),1));
    for j = 1:length(a)
        t1 = bounds{signal}(a(j));
        t2 = bounds{signal}(a(j)+1);
        v1 = envelopes{signal}(t1);
        v2 = envelopes{signal}(t2);
        
        if(v1 < v2)% if a(j):th boundary is at a smaller minima
            torem(a(j)+1) = 1;
        else
            torem(a(j)) = 1;
        end             
    end    
    bounds_t{signal}(torem) = [];
    bounds{signal}(torem) = [];
    end
end