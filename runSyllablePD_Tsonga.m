% Scripts for running syllable-based word learning for Zero Resource Speech
% Challenge @ Interspeech'2015.
%
% Please cite the following paper if you make public use of these codes:
% O. Rasanen, G. Doyle & M.C. Frank: "Unsupervised word discovery from 
% speech using automatic segmentation into syllable-like units," 
% submitted to Interspeech-2015, Dresden, Germany, 2015.
%
% See the same paper for a description of the system.
%
% (c) Okko Rasanen, 2015
%
% Please send all comments and bug-reports to okko.rasanen@aalto.fi .
%
% v0.11, 23.3.2015
% 
% % Note: you need wav-files for the challenge subset of Xitsonga corpus.

%% Basic parameters

n_signals = 138;         % 4058 corresponds to full set, 138 to first talker etc.
n_cuts = 5;            % How many sub-segments per syllable? (default: 5)
cluster_proportion = 0.30; % How many syllable clusters relative to the number of syllable tokens
ngram_orders = [1,2,3]; % N-gram orders to decode as word hypothesis

folder = '/Users/orasanen/speechdb/zerospeech/xitsonga_wavs/';       % Location of audio files

savefile = 'patterns_tsonga.txt';

% Which segmentation algorithm to use: 'oscillator', 'EnvMin', 'VSeg'

segment_algorithm = 'oscillator';

verbose = 1;

%/////////////////////////////////////////////////////////////////////////
%% Get signal filenames and talker IDs for each signal

% Get filenames for audio files
a = dir([folder '*.wav']);
signal_filenames = cell(length(a),1);
for k = 1:length(a)
    signal_filenames{k} = [folder a(k).name];
end

% Get talker IDs for each signal.
% Variable "talkers" keeps track which talker speaks which of the signals.
talker_names = cell(length(signal_filenames),1);
for k = 1:length(signal_filenames)
   talker_names{k} = signal_filenames{k}(end-12:end-9);  % Tsonga
end
unique_talkers = unique(talker_names);
talkers = zeros(length(signal_filenames),1);
for k = 1:length(signal_filenames)
   talkers(k) = find(strcmp(unique_talkers,signal_filenames{k}(end-12:end-9))); % Tsonga
end

n_talkers = length(unique_talkers);

%/////////////////////////////////////////////////////////////////////////
%% Calculate (or load) signal features at 100 Hz sampling rate
if(verbose)
    fprintf('extracting features...');
end

if(~exist('data/feats_tsonga.mat','file'))
    % Add your own feature extraction algorithm here
    F = getFeatures(signal_filenames); % MFCCs
    save('data/feats_tsonga.mat','F');
else
    load('data/feats_tsonga.mat','F');
end

if(verbose)
    fprintf(' Done!\n');
end

%/////////////////////////////////////////////////////////////////////////
% Load audio files and compute full-wave rectified amplitude envelopes at
% 1000 Hz
if(verbose)
    fprintf('calculating envelopes...');
end

x = cell(n_signals,1);
original_audio = cell(n_signals,1);
tmp = cell(n_signals,1);
for signal = 1:n_signals
    [original_audio{signal},audio_fs] = wavread(signal_filenames{signal});
    x{signal} = abs(original_audio{signal});
    x{signal} = resample(x{signal},1000,audio_fs);
    fs = 1000;
end

% Compute low-pass minimum-phase shift envelope using mutual information
% based temporal filter (Rasanen & Laine, 2012).
% Scales envelope between 0-1, limiting recording artifacts to 4 SDs
% above the mean amplitude.

load mif_filter.mat mif_filter
envelopes = cell(n_signals,1);
for signal = 1:n_signals
    envelopes{signal} = filter(mif_filter,1,abs(x{signal}));
    envelopes{signal}(envelopes{signal} > mean(envelopes{signal})+4.*std(envelopes{signal})) = mean(envelopes{signal})+4.*std(envelopes{signal});
    envelopes{signal} = envelopes{signal}./max(envelopes{signal});
end

if(verbose)
    fprintf(' Done!\n');
end

%/////////////////////////////////////////////////////////////////////////
% Perform syllable segmentation
if(verbose)
    fprintf('performing syllable segmentation with %s...\n',segment_algorithm);
end

if(strcmp(segment_algorithm,'oscillator'))
    T = 0.25; % oscillator period
    c = 0.07957; % 8 Hz
    [bounds,bounds_t,osc_env] = oscillatorSylSeg(envelopes,c,T);
elseif(strcmp(segment_algorithm,'EnvMin'))
    [bounds,bounds_t] = getEnvelopeBounds(envelopes);
elseif(strcmp(segment_algorithm,'VSeg'))
    [bounds,bounds_t] = getVillingBounds(original_audio);
else
    error('unknown segmentation algorithm');
end

% bounds: segment boundaries in 1000 Hz samples
% bounds_t: segment boundaries in seconds

% Draw segmentation example

if(verbose)
    fprintf('syllable segmentation Done!\n');
end

%/////////////////////////////////////////////////////////////////////////
% Describe syllable segments using MFCC features
%
% Runs VAD on the data to exclude non-speech segments. This is simulated here
% using evaluation intervals of the Zerospeech challenge.
% Works also without non-speech masking but increases computation time due
% to extra frames (just scales up the number of tokens and thereby the
% number of clusters, but the non-voice regions get excluded in any case
% due to the cropping to evaluation intervals at the end).

if(verbose)
    fprintf('getting syllable-specific features...');
end

vadmasks_100hz = runVAD(F,100,'tsonga');     % Get non-speech and speech regions with a VAD

S = cell(n_signals,1); % <-- this is where features of each syllable from each signal goes

feat_comps_to_take = 1:12; % Which coefficients of the original features to include

span = 0;       % How many frames to include outside segment boundaries
for signal = 1:n_signals
    
    bounds_lp = round(bounds_t{signal}.*100);
    S{signal} = cell(length(bounds_lp)-1,1);
    n = 1;
    toremove = [];
    for k = 1:length(bounds_lp)-1
        
        seglen = bounds_lp(k+1)-bounds_lp(k);
        S{signal}{n} = F{signal}(max(1,bounds_lp(k)-span):min(bounds_lp(k+1)+span,size(F{1},1)),feat_comps_to_take);
        
        if(seglen > 100) % Treat > 1s segments as zeros ("silence")
            S{signal}{n} = zeros(size(S{signal}{n}));
        end
        
        % Use VAD to determine if the current syllable contains only non-speeh frames
        if(sum(vadmasks_100hz{signal}(max(1,bounds_lp(k)-span):min(bounds_lp(k+1)+span,size(F{1},1)))) == 0)
           toremove = [toremove;n];
        end
        n = n+1;
        
    end
    S{signal} = S{signal}(1:n-1);
    
    % Remove frames that do not contain speech (based on VAD)
    S{signal}(toremove) = [];
    bounds_t{signal}(toremove) = [];
    bounds{signal}(toremove) = [];
end

%/////////////////////////////////////////////////////////////////////////
% Convert variable-length feature representations ("spectrograms") into
% fixed length feature vectors by dividing the syllable into N sub-segments
% and averaging feature vectors over those segments, finally concatenating
% all N representations into one long feature vector.

% Store talker specific syllable features

f = cell(n_signals,1); % <-- this is where syllable feature vectors go

featlen = size(S{1}{1},2);

for signal = 1:n_signals
    f{signal} = zeros(length(S{signal}),n_cuts*featlen+1);
    
    for k = 1:length(S{signal})
        
        % Divide into N sub-segments uniformly in time
        
        len = size(S{signal}{k},1);
        
        if(len > 1)
            % Get sub-segment boundaries
            local_bounds = zeros(n_cuts+1,1);
            local_bounds(1) = 1;
            local_bounds(end) = len;
            for cc = 1:n_cuts-1
                local_bounds(cc+1) = max(1,round(cc/n_cuts*len));
            end
            
            % Get mean features across sub-segments
            for t = 1:length(local_bounds)-1
                f{signal}(k,(t-1)*featlen+1:t*featlen) = mean(S{signal}{k}(local_bounds(t):local_bounds(t+1),:));
            end
            f{signal}(k,end) = log(len).*n_cuts/3;      % Add scaled log-duration
        end
    end
end

if(verbose)
    fprintf(' Done!\n');
end

%/////////////////////////////////////////////////////////////////////////
% Cluster syllable features into syllable categories.
% Clustering is done separately for each talker.
%
% Currently uses kmeans with random sampling initialization
if(verbose)
    fprintf('clustering syllables...\n');
end

codebooks = cell(n_talkers,1);   % <-- talker-specific codebooks

for talker_id = 1:n_talkers
    a = find(talkers == talker_id);
    
    % Concatenate all segment features for the given talker across signals
    datacloud = [];
    for j = 1:length(a)
        if(a(j) <= length(f))
            datacloud = [datacloud;f{a(j)}];
        end
    end
    
    % Run clustering if talker_id is included in the list of processed signals
    if(~isempty(datacloud))
        tmp = find(sum(datacloud(:,1:10),2) == 0);  % Remove zero-vectors (if any)
        datacloud(tmp,:) = [];
        
        totsyls = size(datacloud,1);   % Total number of syllable segments for the talker
        
        [~,codebooks{talker_id}] = kmeans(datacloud+randn(size(datacloud)).*0.001,round(totsyls*cluster_proportion),'start','sample','MaxIter',100,'EmptyAction','drop');
        
    end
    procbar(talker_id,n_talkers);
end

% Add cluster with zero-valued vectors (for zero feature vectors).
for k = 1:length(codebooks)
    if(~isempty(codebooks{k}))
        codebooks{k} = [codebooks{k};zeros(1,size(codebooks{k},2))];
    end
end

% Get VQ indices for syllables
vq_sequences = cell(n_signals,1);
for signal = 1:n_signals
    talker_id = talkers(signal);
    vq_sequences{signal} = getLabels(f{signal},codebooks{talker_id}',1);
end

if(verbose)
    fprintf(' Done!\n');
end

%/////////////////////////////////////////////////////////////////////////
% Decode word hypotheses from the VQ-sequences by finding n-grams
%
% Note: the n-gram algorithm in decodeSyllables() uses fast matrix-based
% n-gram indexing, but requires quite a lot of memory to initialize sparse
% matrices. Replace with your own algorithm for high-order n-grams in low
% RAM systems.

if(verbose)
    fprintf('decoding n-grams...');
end

[classes,times,signal_names] = decodeSyllables(vq_sequences,bounds_t,talkers(1:n_signals),ngram_orders);

if(verbose)
    fprintf(' Done!\n');
end

%[classes,times,signal_names] = pruneClusters(S,vq_sequences,talkers,classes,times,signal_names,unigrammap);
%/////////////////////////////////////////////////////////////////////////
% Crop full algorithm output to Zerospeech evaluation intervals
[times_eva,classes_eva,signal_names_eva] = croptointervalsTsonga(times,classes,signal_names);

%/////////////////////////////////////////////////////////////////////////
% Pruning of short segments or classes with only one token (off by default)
%[times_eva,classes_eva,signal_names_eva] = removeShorts(times_eva,classes_eva,signal_names_eva,0.08);
%[times_eva,classes_eva,signal_names_eva] = removeSingles(times_eva,classes_eva,signal_names_eva);

%/////////////////////////////////////////////////////////////////////////
% Write patterns to Zerospeech evaluation kit-compatible output file.
writeResultOutputTsonga(savefile,signal_names_eva,times_eva,classes_eva);

