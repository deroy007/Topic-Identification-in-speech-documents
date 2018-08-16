n_signals = 5;         % 68 corresponds to the full set, 5 to the first talker etc.
n_cuts = 5;            % How many sub-segments per syllable? (default: 5)
cluster_proportion = 0.30; % How many syllable clusters relative to the number of syllable tokens
ngram_orders = [1,2,3]; % N-gram orders to decode as word hypothesis

folder = '/home/aniket/Singlespeaker_data/WAV';
       % Location of audio files

savefile = 'patterns_singlespeaker.txt';

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
disp(signal_filenames)
%/////////////////////////////////////////////////////////////////////////
%% Calculate (or load) signal features at 100 Hz sampling rate
if(verbose)
    fprintf('extracting features...');
end

if(~exist('data/feats.mat','file'))
    % Add your own feature extraction algorithm here
    F = getFeatures(signal_filenames); % MFCCs
    save('data/feats.mat','F');
else
    load('data/feats.mat','F');
end

if(verbose)
    fprintf(' Done!\n');
end
