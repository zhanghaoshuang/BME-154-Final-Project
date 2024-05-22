% BME154L_FinalProject_Question4.m

% Matlab Sample Code
% BME 154 Final Project
tic
%% Initialization

clear all;
close all;

%% Load data

fileid=fopen('BME154L_S12_PROJECT_ECG.bin'); % code for opening .bin files
data = fread(fileid,inf,'float32');
fclose('all');


%% Define variables

t=data(1:2:end); % (seconds)
ecg= data(2:2:end); %(Voltage [mV])
ecg_raw = ecg;

% find number of points for 0.4 seconds. used later for 'minpeakdistance'
num_points = round(0.4/(t(end)/length(t)));

%% Plot of Original ECG

figure(1) 

plot(t(1200:1600),ecg(1200:1600))  % one ecg cycle isolated 
xlabel('Time (s)')                  
ylabel('Voltage (mV)')
title('ECG Data')

print -dpng part4fig1

%% COMPUTING SNR?
% uses function ECGsnr.m
raw_snr = ECGsnr(ecg_raw, 1000:20000); % indices 1000:20000 chosen as a random representation of normal signal


%% Removing Power Noise

% Compute FFT        %used the fourier transform to get the power spectrum

fs=1/mean(diff(t));
ft=fft(ecg);
f = linspace(-fs/2,fs/2,length(ft));
shiftFT = fftshift(ft);

figure
plot(f, 20*log10(abs(shiftFT)./max(abs(shiftFT))))
axis([0 100 -80 0]);
xlabel('Frequency (Hz)')
ylabel('Relative Power (dB)')
title('ECG FFT Power Spectrum');

print -dpng part4fig2

% Removing Artifact with notch filter

for index = 1:length(f)
    if abs(f(index))>59.75 && abs(f(index))<60.25
        shiftFT(index) = 0;
    end
end

figure
plot(f, (20*log10(abs(shiftFT)./max(abs(shiftFT)))))
axis([0 100 -80 0])
xlabel('Frequency (Hz)')
ylabel('Relative Power (dB)')
title('ECG FFT Power Spectrum w/Noise Removed');

print -dpng part4fig3

ecg = real(ifft(fftshift(shiftFT)));

% Demonstration of Signal Improvement

figure
subplot(2,1,1)
plot(t, ecg_raw)
axis([620 621.5 -4 10])
xlabel('Time (s)'); ylabel('Voltage (mV)');
title('Raw ECG Data')

subplot(2,1,2)
plot(t, ecg)
axis([620 621.5 -4 10])
xlabel('Time (s)'); ylabel('Voltage (mV)');
title('ECG Data w/ Power Noise Artifact Removed')

print -dpng part4fig4

%% Remove DC Offset and Linear Tilt From Data           

% Removing any Offset

p = polyfit (t, ecg, 1); %fit a line to the data
tilt = p(1)*t+p(2);
clean_ecg = ecg-tilt;

figure % this graph shows the signal with tilt removed
plot(t,clean_ecg);
xlabel('Time (s)');
ylabel('Voltage (mV)');
title ('ECG Data with DC Offset and Linear Tilt Removed')
axis([620 621.5 -4 10])

print -dpng part4fig5

% Compute Clean FFT

ft_clean = fft(clean_ecg);

figure
plot(f, fftshift(20*log10(abs(ft_clean)./max(abs(ft_clean)))));
axis([0 100 -80 0]);
xlabel('Frequency (Hz)');
ylabel('Relative Power (dB)')
title('Clean ECG FFT Power Spectrum')

print -dpng part4fig6

noartifact_snr = ECGsnr(clean_ecg,1000:20000);


%% Bradycardia Identification

[~, locs] = findpeaks(clean_ecg, 'MINPEAKHEIGHT', 0.7*max(clean_ecg), 'MINPEAKDISTANCE', num_points);
beattime = zeros(1, length(locs)-1);
for index = 1:length(locs)-1
    beattime(index) = t(locs(index+1))-t(locs(index));
end

figure
plot(t(locs), clean_ecg(locs), 'rx')
hold on
plot(t, clean_ecg)
xlabel('Time (s)'); ylabel('Voltage (mV)');
title('Clean ECG with Peaks Identified')
axis([620 630 -3 10])

print -dpng part4fig7

HR = 60./beattime;

figure
plot(t(locs(1:end-1)), HR, '-')
xlabel('Time (s)'); ylabel('Instantaneous Heart Rate (bpm)')
title('Heart Rate for Individual Beats')

print -dpng part4fig8

counter = 0;
for index = 1:length(HR)-1
    if HR(index) < 50 && HR(index+1) < 50
        counter = counter+1;
        bcardia_locs(counter) = locs(index-1);
    end
end

% t(bcardia_locs) gives the times when bradycardia occurs.

%% Maximizing SNR of ECG data using Boxcar Averager

boxcar_windowsizes= [2 3 5 6 7 8 9 10 15 20]; % number of samples

for i = 1:length(boxcar_windowsizes)
    ecg_boxcar(i,:) = conv2(clean_ecg,ones(boxcar_windowsizes(i),1),'same')./boxcar_windowsizes(i);
end;
% 
% for i=1:length(boxcar_windowsizes)
%     figure(20+i)
%     plot(t,ecg_boxcar(i,:))
%     axis([620 621.5 -2 8])
%     title(sprintf('Boxcar Window = %i samples', boxcar_windowsizes(i)));
%     xlabel('Time (seconds)'); ylabel('Voltage (mV)');
% end



% % From observation it appears a boxcar window size of 8 is ideal for
% % smoothing out noise but maintaining signal power

boxcarECG = ecg_boxcar(6,:);

figure
subplot(2,1,1)
plot(t, clean_ecg)
axis([620 630 -2 8])
xlabel('Time (seconds)'); ylabel('Voltage (mV)');
title(sprintf('Clean ECG before Boxcar'));
subplot(2,1,2)
plot(t,boxcarECG)
axis([620 630 -2 8])
xlabel('Time (seconds)'); ylabel('Voltage (mV)');
title(sprintf('Boxcar Window = 8 samples'));

print -dpng part4fig9

boxcar_snr = ECGsnr(boxcarECG,1000:20000);


%% Cross Correlation Heartbeat Detection (Didn't end up using this method)

% % % % ecg_reference = boxcarECG(1200:1600); % based on observation
% % % % [xco, lags] = xcorr(boxcarECG,ecg_reference);
% % % % 
% % % % % % % % % % figure(10)
% % % % % % % % % % 
% % % % % % % % % % hold on;
% % % % % % % % % % plot(t, clean_ecg);
% % % % % % % % % % xlabel('Time (s)');
% % % % % % % % % % ylabel('Voltage (mV)')
% % % % % % % % % % title('ECG Signal with Cross Correlation Heartbeat Detection')
% % % % 
% % % % % Find Local Maxima in Cross Correlation
% % % % % get rid of noise
% % % % 
% % % % % % % w=1500;
% % % % % % % xco_bca = conv2(xco,ones(w,1),'same')./w;
% % % % % find peaks
% % % % figure(10)
% % % % plot((t)+t(201), (xco(length(t):end)))
% % % % hold on
% % % % [pks,pk_locs] = findpeaks(xco,'minpeakheight',631, 'minpeakdistance', 200); % observation threshold
% % % % [pks2,pk_locs2] = findpeaks(xco,'minpeakheight',520,'minpeakdistance', 200); % observation threshold
% % % % 
% % % % pk_locs = pk_locs+200-length(t);
% % % % pk_locs2 = pk_locs2+200-length(t);
% % % % 
% % % % PVClocs = setdiff(pk_locs2, pk_locs);
% % % % 
% % % % % heartbeats = t(pk_locs);
% % % % % plot(t(pk_locs), xco(pk_locs), 'rx')
% % % % 
% % % % % unsure how to perfectly lineup correlation, actual data
% % % % plot((t(pk_locs)), (pks), 'rx')
% % % % 
% % % % xlabel('Time (s)'); ylabel('Correlation Area');
% % % % title('ECG data Correlated w/ Reference ECG, Peaks Identified')
% % % % 
% % % % % plot(heartbeats, ones(length(heartbeats),1),'rx','MarkerSize', 4,'LineWidth',1);



%% Negative Threshold PVC Detection

% Make a fair assumptiont that no PVC will occur in the first few cycles of the data


% can be assumed min(boxcarECG) will be negative as the ECG has been centered at zero


threshold = min(boxcarECG(1:2000))-3;

[dips,dip_locs] = findpeaks(-boxcarECG, 'minpeakheight', -threshold, 'minpeakdistance',num_points);
dips = -dips;

% PVC times are given by t(dip_locs)

%% AVERAGE HEART RATE OF NSR

A = [];
B = [];

for index = 1:length(HR)
    % selecting PVC times to later be removed
    if sum(abs(t(locs(index+1))-t(dip_locs))<0.3)>0 || sum(abs(t(locs(index))-t(dip_locs))<0.3)>0
        A = [A index];
    end
end

for index = 3:length(HR)+2
    % selecting bradycardia times to later be removed
    if sum(t(locs(index-1))==t(bcardia_locs))>0||sum(t(locs(index-2))==t(bcardia_locs))>0
        B = [B index];
    end
end

NSRHR = HR;
NSRHR([A B]) = [];

tNSR = t(locs(1:end-1));
tNSR([A B]) = [];


figure
plot(NSRHR,'-*')
xlabel('Time (s)'); ylabel('Instantaneous Heart Rate (bpm)')
title('NSR Region Instantaneous Heart Rate over Time')
print -dpng part4fig10

% overall mean Heart Rate
meanHR = mean(NSRHR);



%% 60 second running average HR

% find 60 seconds (roughly)
register = ceil(60*(length(locs))/t(end));

% could have used convolution as follows:
% run_avgHR = conv(NSRHR,ones(register,1)/register,'same');
 
run_avgHR = zeros(1,length(NSRHR)-register);
run_avgHR(1) = mean(NSRHR(1:register));

% register acts
for index3 = 2:length(NSRHR)-register
    run_avgHR(index3) = run_avgHR(index3-1)+(-NSRHR(index3-1)+NSRHR(index3+register))/register;
end

figure
plot(t(locs(1:length(run_avgHR))),run_avgHR)
xlabel('Time (s)'); ylabel('Heart Rate (bpm)')
title('Sixty Second Running Average of Heart Rate over Time')
print -dpng part4fig11


%% Average Over Different Heartbeat Using Phase Alignment

% Estimated the PR interval using the method of QRS peak alignment and
% flipping

maxqrs = max(clean_ecg); % find maximum qrs value

[peaks, qrsloc] = findpeaks(clean_ecg, 'minpeakheight', 0.6*maxqrs, 'minpeakdistance', num_points);
[minpeaks, ~] = findpeaks(-clean_ecg, 'minpeakheight', 2, 'minpeakdistance', num_points);

% parse signals and zero-pad
meshint = max(diff(qrsloc)) + 1; % finding the largest interval
ker = zeros(1, meshint); ker(end) = 1; % create a delta function to correlate
parsearray = zeros([length(qrsloc), length(xcorr(ker, ker))]); % array containing all parsed cycles

for k = 1:length(qrsloc)-1   
    % zero-padding using xcorr. parsearray rows are each qrs to qrs cycle
    parsearray(k, :) = xcorr(flipud(clean_ecg(qrsloc(k):qrsloc(k+1))), ker);
end

% cut off extraneous zero-padding
parsearray = parsearray(:, 1:meshint);

% create an array where each row is the mean cycle
avgarray = ones([length(qrsloc), 1])*mean(parsearray);

% parsed cycle adjusted by mean
adjarray = parsearray - avgarray;

figure
plot(mean(diff(t))*(1:length(mean(parsearray))),fliplr(mean(parsearray)))
xlabel('Scaled Time (s)'); ylabel('Voltage (mV)');
title('QRS Peaks aligned, averaged, and flipped')
axis([1.6 2.8 -1 7])
print -dpng fig20


%% The following is based on PS 7 solutions and is used to check the mean PR interval

% cross corrlation heart beat detection

ecg_reference = clean_ecg(1200:1600); % based on observation
[xco,lags] = xcorr(clean_ecg,ecg_reference);

figure;
hold on;
plot((t)+t(201), (xco(length(t):end)))
xlabel('Time (s)');
ylabel('Voltage (mV');
title('ECG Signal with Cross Correlation Heartbeat Detection')

% find peaks
[pks,pk_locs]=findpeaks(xco,'minpeakheight',500,'minpeakdistance',num_points);

pk_locs = pk_locs+200-length(t);


plot((t(pk_locs)), pks, 'rx')
print -dpng part4fig12



% average over different heart beat using correlation-based phase alignment
ecg_width = length(ecg_reference)-1; % based on our earlier reference signal

for i=1:length(pk_locs),
    
    if((pk_locs(i)- ecg_width/2) < 1)
        zero_pad = zeros(pk_locs(i),1);
        ecg_corr_ave(i,:) = [zero_pad; clean_ecg(1:pk_locs(i)+ecg_width/2)];
    
    elseif((pk_locs(i) + ecg_width/2) > length(clean_ecg)),
        zero_pad = zeros(ecg_width/2-(length(clean_ecg)-pk_locs(i)),1);
        ecg_corr_ave(i,:) = [clean_ecg(pk_locs(i)-ecg_width/2:end); zero_pad];
        
    else
        ecg_corr_ave(i,:) = clean_ecg((pk_locs(i)-ecg_width/2):(pk_locs(i)+ecg_width/2));
        
    end

end


figure
imagesc(t(1:size(ecg_corr_ave,2)),1:length(locs), ecg_corr_ave);
xlabel('Time (s)');
ylabel('Heart Beat');
title('Cross Correlation, Phase-Aligned ECG Signals');
print -dpng part4fig13



figure
plot(t(1:size(ecg_corr_ave,2)),mean(ecg_corr_ave,1)');
xlabel('Time (s)');
ylabel('Voltage (mV)');
title('Mean ECG Signal (Xcorr, Phase-Aligned)');
print -dpng part4fig14


toc