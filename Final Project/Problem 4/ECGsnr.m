%ECGsnr.m

% BME 154 final project
% Problem 4

% Finding ECG SNR

function snr = ECGsnr(signal, ind)
%% load data
n = ind;
ecg = signal(ind);

% % % fileID = fopen('ecgdata.bin');
% % % data = fread(fileID, inf, 'float32');
% % % fclose('all');
% % % 
% % % t = data(1:2:3400*2);
% % % ecg = data(2:2:3400*2);
% % % ecg = detrend(ecg, 'linear');

% % % cleanecg = load('cleanECG.mat')
% % % t = cleanecg.t(n);
% % % ecg = cleanecg.clean_ecg(n);

%% COMPUTING SNR
maxqrs = max(ecg); % find maximum qrs value

[~, qrsloc] = findpeaks(ecg, 'minpeakheight', 0.8*maxqrs, 'minpeakdistance', 100);

% parse signals and zero-pad
meshint = max(diff(qrsloc)) + 1; % finding the largest interval
ker = zeros(1, meshint); ker(end) = 1; % create a delta function to correlate
parsearray = zeros([length(qrsloc), length(xcorr(ker, ker))]); % array containing all parsed cycles

for k = 1:length(qrsloc)-1   
    % zero-padding using xcorr. parsearray rows are each qrs to qrs cycle
    parsearray(k, :) = xcorr(ecg(qrsloc(k):qrsloc(k+1)), ker);
end

% cut off extraneous zero-padding
parsearray = parsearray(:, 1:meshint);
avgsig = mean(parsearray);

% create an array where each row is the mean cycle
avgarray = ones([length(qrsloc), 1])*avgsig;

% parsed cycle adjusted by mean
noisearray = parsearray - avgarray;

% % % figure
% % % plot(mean(parsearray))
% % % title('averaged signal')
% % % 
% % % figure
% % % plot(parsearray')
% % % title('signals')
% % % 
% % % figure
% % % plot(noisearray(1:5, :)')
% % % title('noise signal')

figure
subplot(2,1,1)
plot(mean(parsearray))
title('averaged signal')
subplot(2,1,2)
plot(noisearray(1:5,:)')
title('noise signal')
axis([0 800 -.6 0.6])
print -dtiff fig70

% fair to assume noise can be measured within an area we assume to be
% centered at zero
noise_adj = detrend(noisearray(:,100:400)');

% % % figure
% % % plot(noise_adj)
% % % title('noise centered at zero')

% find the standard deviation of the adjusted cycles
noise = sqrt(mean2(noise_adj.^2)); %2*mean(std(noisearray(:,200:400)'))

% find the average maximum signal range
signal =  sqrt(sum(sum(parsearray.^2))/length(ecg)'); %max(avgsig) - min(avgsig) %sqrt(mean(ecg.^2))

snr = 20*log10(signal./noise);

end