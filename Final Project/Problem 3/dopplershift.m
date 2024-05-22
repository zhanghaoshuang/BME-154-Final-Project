% dopplershift.m
% BME 154 Final Project
% Problem 3 Optional

% Use Doppler shift to calculate velocity of
%% load data
clear all;
close all;

file = 'A';
a = load(['mmode', lower(file), '.mat']); % load filename
b = load('mmodeb.mat');
c = load('mmodec.mat');

rows = dot(size(a.mmode), [1, 0]); % number of rows
cols = dot(size(a.mmode), [0, 1]); % number of columns
a.c = 1540;     % set the speed of sound in tissue

% fourier transform of the time series data at a specific depth
fftdataA = fftshift(fft(a.mmode'));
cf = max(abs(fftdataA)); % maximum values of each power spectrum
cfind = zeros(size(cf)); % vector of indices of center frequency

% identify the center frequency for each row
for m = 1:length(cf)
   cfind(m) = dot(find(abs(fftdataA(:,m)) == cf(m)), [0, 1]);
end

% finding the frequency range. Sampling frequency of the time
% series data is the pulse repeititon frequency
freq = linspace(-a.prf/2, a.prf/2, cols);

% % % % plot all power spectrums
% % % figure
% % % plot(freq, 20*log10(abs(fftdataA)/max(max(abs(fftdataA)))))
% % % title(['(Dataset ', file, ') Power spectrum of all time series signal'])
% % % xlabel('Frequency (Hz)')
% % % ylabel('Magnitude (dB)')
% % % 
% % % print -dpng powerspectrumtimeA

avgcf = mean(freq(cfind)); % finding the average center frequency
avgcf - freq(cfind(1));

v = a.c/(2*a.f0)*avgcf % calculate the velocities

%% Calculating direction of motion

fftdataRF = mean(fftshift(fft(a.mmode))');
fftdataRFb = mean(fftshift(fft(b.mmode))');
fftdataRFc = mean(fftshift(fft(c.mmode))');

cfRF = max(abs(fftdataRF));
cfRFb = max(abs(fftdataRFb));
cfRFc = max(abs(fftdataRFc));

cfindRF = dot(find(abs(fftdataRF) == cfRF), [0, 1]);
cfindRFb = dot(find(abs(fftdataRFb) == cfRFb), [0, 1]);
cfindRFc = dot(find(abs(fftdataRFc) == cfRFc), [0, 1]);

freqRF = linspace(-a.fs/2, a.fs/2, rows);
dir = freqRF(cfindRF) > a.f0 % test whether center frequency is greater than f0

% % % % plot all power spectrums of each RF line
% % % figure
% % % subplot(3,1,1)
% % % plot(freqRF, abs(fftdataRF), 'k-', [a.f0, a.f0], [0, 1.1*cfRF], 'r-', ... 
% % %     freqRF(cfindRF), cfRF, 'ro')
% % % title('(Dataset A) Avg Fourier transform')
% % % xlabel('Frequency (Hz)')
% % % ylabel('Magnitude')
% % % axis([0, 0.5e7, 0, 1.2*cfRF])
% % % 
% % % subplot(3,1,2)
% % % plot(freqRF, abs(fftdataRFb), 'k-', [a.f0, a.f0], [0, 1.1*cfRFb], 'r-', ... 
% % %     freqRF(cfindRFb), cfRFb, 'ro')
% % % title('(Dataset B) Avg Fourier transform')
% % % xlabel('Frequency (Hz)')
% % % ylabel('Magnitude')
% % % axis([0, 0.5e7, 0, 1.2*cfRFb])
% % % 
% % % subplot(3,1,3)
% % % plot(freqRF, abs(fftdataRFc), 'k-', [a.f0, a.f0], [0, 1.1*cfRFc], 'r-', ... 
% % %     freqRF(cfindRFc), cfRFc, 'ro')
% % % title('(Dataset C) Avg Fourier transform')
% % % xlabel('Frequency (Hz)')
% % % ylabel('Magnitude')
% % % axis([0, 0.5e7, 0, 1.2*cfRFc])
% % % 
% % % print -dpng dopplershiftdirection

