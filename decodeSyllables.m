function [classes,times,signal_names] = decodeSyllables(vq_sequences,bounds_t,talkers,orders)
%function [classes,times,signal_names] = decodeSyllables(vq_sequences,bounds_t,talkers,orders)
%
% Decodes vq_sequences into word-pattern hypotheses by starting from the
% highest n-gram order defined in "orders", and starting from the
% highest-frequency n-grams in each order.
%
%   Inputs:
%       vq_sequences        : discrete sequences of syllables for
%                             each speech signal. Nx1 cell-array for
%                             N signals.
%       bounds_t            : onset timestamps for syllables. Nx1
%                             cell-array for N signals.
%       talkers             : talker-id:s for signals (vector of Nx1)
%       orders              : n-gram orders to decode
%
%   Outputs 
%       classes             : Mx1 vector of pattern class-IDs
%       times               : Mx2 matrix of pattern onset and offset times
%       signal_names        : Mx1 cell array of signal names for each
%                             pattern

classes = zeros(10000,1);
times = zeros(10000,2);
signal_number = zeros(10000,1);
n = 1;
classnum = 1;

talker_ids = unique(talkers);

if nargin <4
    orders = 1:3;
end

% Minimum token frequency for each n-gram order.
minfreqs = [1,2,2,2]; 

orders = sort(orders,'descend');


for talker = 1:length(talker_ids)
    % Get signals for the given talker and concatenate into a long sequence
    talker_id = talker_ids(talker);        
    a = find(talkers == talker_id);
    fullseq = [];
    for k = 1:length(a)
        if(a(k) <= length(vq_sequences))
            fullseq = [fullseq;vq_sequences{a(k)}];
        end
    end   
    if(~isempty(fullseq))                
        for order_to_decode = orders
            % Find n-grams and their frequencies in the sequence
            [ngrams,freqs] = getngrams(fullseq,order_to_decode,max(fullseq));
                        
            if(max(freqs) >= minfreqs(order_to_decode))
                % Get n-gram locations in the decreasing order of frequency
                % and remove participating VQ-states from further n-grams.  
                % Assign unique class-ID for each ngram.
                for frekki = max(freqs):-1:minfreqs(order_to_decode)
                    b = find(freqs == frekki);
                    
                    for id = 1:length(b)
                        if(order_to_decode > 1)
                            grammi = b(id);
                            ngram = getngram(ngrams(grammi),max(fullseq));
                            
                        else
                            ngram = ngrams(b(id));
                        end
                        
                        % n-gram locations in each original sequence
                        locs = cell(length(a),1);
                        for k = 1:length(a)
                            if(length(vq_sequences{a(k)}) > 1)
                                if(order_to_decode > 1)
                                    locs{k} = findngram(vq_sequences{a(k)},ngram);
                                else
                                    locs{k} = find(vq_sequences{a(k)} == ngram);
                                end
                            end
                        end
                        if(sum(cellfun(@length,locs)) > 0)
                            for k = 1:length(a)
                                if(a(k) < length(vq_sequences))
                                    
                                    for j = 1:length(locs{k})
                                        if(locs{k}(j) < length(bounds_t{a(k)})) 
                                            onset = bounds_t{a(k)}(locs{k}(j));
                                            offset = bounds_t{a(k)}(locs{k}(j)+order_to_decode);
                                            times(n,:) = [onset offset];
                                            classes(n) = classnum;
                                            signal_number(n) = a(k);
                                            vq_sequences{a(k)}(locs{k}(j):locs{k}(j)+order_to_decode-1) = NaN;
                                            n = n+1;
                                        end
                                    end
                                    
                                end
                                
                            end                            
                           
                            classnum = classnum+1; % 
                        end
                    end
                end
            end
        end
    end    
end

% Cut off unused elements
classes = classes(1:n-1);
times = times(1:n-1,:);
signal_number = signal_number(1:n-1);

signal_names = convertSignNames(signal_number);
