function F = getFeatures(signal_filenames) 
% function F = getFeatures(signal_filenames) 
%
% Computes MFCCs using scripts by Saikat Chatterjee (KTH)

WinShift = 0.01;
WinLen = 0.025;

MFCC_flag=1;
MMFCC_flag=0;

% Dynamic parameter. If this flag is set, then the output vector will be 39
% dimensional as E+FeatureVec+Vel+Acc. If this flag is off, then the output
% will be only 12 dimensional feature.

Dynamic_flag=1;

% Cepstrum Mean and Normalization (CMVN). If set, then output feature will be Cepstrum Mean and
% Variance Normalized (CMVN). The CMVN is carried out in
% utterance-by-utterance basis (i.e. file-by-file).

CMVN_flag=1;

% Specification about Fs, Window Length and Window Shift
SigOprFs=16000;                  % Signal Operating Sampling Frequency

% Noise parameters

NoiseType = 'white';   % (babble, pink, volvo and white
SNR = 1000;        % SNR in dB  (Choose SNR=1000 for clean testing)
NoiseDBPath = '';
train_flag = 1; % set 1 if clean training

F = cell(length(signal_filenames),1);
N = length(signal_filenames);


for k = 1:length(signal_filenames);    
    
    filename = signal_filenames{k};
    
    
    [input,fs] = wavread(filename);
    
    
    if(fs ~= 16000)
        input = resample(input,16000,fs);
        fs = 16000;
    end
    
    
    O = make_MFCC_And_MMFCC2_features(input,WinLen,WinShift,fs,SigOprFs,MFCC_flag,MMFCC_flag,Dynamic_flag,CMVN_flag,NoiseType,SNR,NoiseDBPath,train_flag);
    
    F{k} = O;
    
    
    procbar(k,N);
    
    
end
fprintf('\n');