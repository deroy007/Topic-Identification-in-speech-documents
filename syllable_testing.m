

clear all;
close all;
clc;

[y Fs] = audioread('/home/aniket/Hindi/02/wav/4_02.wav');


%figure;
%plot(y)
soundsc(y,Fs);

xx = abs(y);
x = resample(xx,1000,16000);
fs = 1000;

%figure;
%subplot(2,1,1);
%t = 1/fs:1/fs:(length(x)/fs);
%plot(t,x,'b');
%subplot(2,1,2); 
%t = 1/fs:1/fs:(length(x)/fs);
%plot(t,x,'r');


load mif_filter.mat mif_filter
load('mif_filter.mat')
envelopes = filter(mif_filter,1,abs(x));
envelopes(envelopes > mean(envelopes)+4.*std(envelopes)) = mean(envelopes)+4.*std(envelopes);
envelopes = envelopes./max(envelopes);

%figure;
%subplot(3,1,1);
%xx = resample(original_audio{1},1000,audio_fs);
%t = 1/fs:1/fs:(length(x{1})/fs);
%plot(t,x,'b');
%subplot(3,1,2); 
%t = 1/fs:1/fs:(length(x)/fs);
%plot(t,x,'r');
%subplot(3,1,3);
%plot(t,envelopes,'r');


% yyy = medfilt1(envelopes);
% plot(yyy);
% yyy = medfilt1(envelopes,10);
% plot(yyy);
T = 0.25; % oscillator period
c = 0.07957; % 8 Hz
[bounds,bounds_t,osc_env] = oscillatorSylSeg(envelopes,c,T);
%plot(osc_env{1,1})
% figure;plot(yyy)


[ta1,tfs] = audioread('/home/aniket/Hindi/02/wav/4_02.wav');
tp1=load('/home/aniket/Hindi/posterior1/MFCC/4_02.post');
        
phnFilename = '/home/aniket/Hindi/02/syl/4.syl';
LPA1 = zeros(size(tp1,1),1);
winS = 0.010*tfs;
[FNF,S_exp,S_obt,B2,B23] = obtainFNFmat(phnFilename,LPA1,winS);
Z = zeros(B2(size(B2,1),size(B2,2)),1);
B22 = [1;B2(:,2)];
Z(B22)=1;
t = 1/tfs:1/tfs:(length(ta1)/tfs);
figure,subplot(5,1,1),
bar(t(1:length(Z)),Z,30.0,'r');
set(gca,'FontSize', 14);
ht = text(B22/16000, 0.8*ones(1,length(B22)),[B23']);hold on
hold on,plot(t(1:length(Z)),ta1(1:length(Z)),'b'),grid off,hold on
subplot(5,1,2),
% xx = resample(original_audio{1},1000,audio_fs);
t = 1/fs:1/fs:(length(x)/fs);
plot(t,x,'b');
subplot(5,1,3); 
t = 1/fs:1/fs:(length(x)/fs);
plot(t,x,'r');
subplot(5,1,4); 
plot(t,envelopes,'r');
subplot(5,1,5); 
plot(t,osc_env{1,1},'r');hold on;
bb = bounds{1,1};
for k = 1:length(bb)
   line([bb(k)./1000 bb(k)./1000],[-1 1],'Color','red'); 
end
% bar(t(1:length(Z)),Z,30.0,'r');
% set(gca,'FontSize', 14);
% ht = text(B22/16000, 0.8*ones(1,length(B22)),[B23']);hold on
% h[pks,locs,widths,proms] = findpeaks(envelopes.*-1,'Annotate','Extents');old on,plot(t(1:length(Z)),ta1(1:length(Z)),'b'),grid off,hold on
[pks,locs,widths,proms] = findpeaks(envelopes.*-1,'Annotate','Extents');
figure;
findpeaks(envelopes.*-1,'Annotate','Extents');
