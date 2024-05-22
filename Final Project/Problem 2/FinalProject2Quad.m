%FinalProject2Quad

%% Plotting the RF Data

clear all; close all; clf
load('rf.mat')
figure(1)
imagesc(lat,depth,rfdata)
xlabel('lateral direction (mm)')
ylabel('depth (mm)')
title('Original Radio Frequency Data')
axis image
colormap hot
colorbar

%Calculate Contrast

SignalIntensity=mean(mean(rfdata(854:975,360:521)));
BackGroundIntensity=mean(mean(rfdata(1:300,1:300)));
Contrast=(SignalIntensity-BackGroundIntensity)/(max(max(rfdata))-min(min(rfdata)));


%% Quadrature Demodulation
[row, cols]= size(rfdata);

% Find Power Spectrum to see what frequency to use
T0=(depth(length(depth))-depth(1))*2/1540/1000;
N=row;
ts=T0/N;
fs=1/ts;
freq= fs* (0 : 1/N: 1-(1/N));
DataFT=abs(fft(rfdata));
DataFTNorm=abs(DataFT)/max(max(DataFT));
PowerNorm=20*log10(DataFTNorm);
figure(2)
plot(freq,PowerNorm)
title('Power Spectrum of RF Data Columns (A-Lines)')
xlabel('frequency (Hz)')
ylabel('power (dB')

% Finding the Center Frequency

for k=1:cols
    Index(k)=find(max(PowerNorm(1:round(length(freq)/2),k))==PowerNorm(1:round(length(freq)/2),k));
    CenterFreq(k)=freq(Index(k));
end

CenterFreqAvg=mean(CenterFreq);

t=(depth*2/1000/1540)'; % Finding the right time points to use
SignalCos=zeros(row,cols);
SignalSin=zeros(row,cols);

%Multiply by Sine and Cosine
for i=1:cols
 SignalCos(:,i)=cos(2*pi*CenterFreqAvg*t).*rfdata(:,i);
 SignalSin(:,i)=sin(2*pi*CenterFreqAvg*t).*rfdata(:,i);
end

%FT
SignalCosFT=fft((SignalCos));
SignalSinFT=fft(SignalSin);


figure(4)
plot(freq,abs(SignalCosFT))
title('Amplitude Response of RF Data Multiplied by Cosine')
xlabel('frequency (Hz)')
ylabel('Amplitide')

figure(5)
plot(freq,abs(SignalSinFT))
title('Amplitude Response of RF Data Multiplied by Sine')
xlabel('frequency (Hz)')
ylabel('Amplitide')

% Multiply by Rect/Low Pass Filter
CutoffFreq=0.5e7;
SignalCosFTFiltered=zeros(row,cols);
SignalSinFTFiltered=zeros(row,cols);
for i=1:cols
    SignalCosFTFiltered(:,i)=SignalCosFT(:,i).*(abs(freq)<CutoffFreq)';
    SignalSinFTFiltered(:,i)=SignalSinFT(:,i).*(abs(freq)<CutoffFreq)';
end

figure(6)
plot(freq,abs(SignalCosFTFiltered))
title('Amplitude Response of Filtered RF Data Multiplied by Sine')
xlabel('frequency (Hz)')
ylabel('Amplitide')

figure(7)
plot(freq,abs(SignalSinFTFiltered))
title('Amplitude Response of Filtered RF Data Multiplied by Sine')
xlabel('frequency (Hz)')
ylabel('Amplitide')

% Reconstruct/InverseFT
RealSignal=real(ifft(SignalCosFTFiltered));
ImaginarySignal=real(ifft(SignalSinFTFiltered));

figure(8)
SignalRecon=2*sqrt(RealSignal.^2+ImaginarySignal.^2);
imagesc(lat,depth,SignalRecon)
axis image
colormap hot
colorbar
xlabel('lateral direction (mm)')
ylabel('depth (mm)')
title('Data After Quadrature Demodulation')
caxis([0 4e-21])

% Calculate Contrast
SignalIntensity=mean(mean(SignalRecon(854:975,360:521)));
BackGroundIntensity=mean(mean(SignalRecon(1:300,1:300)));
ContrastDemod=(SignalIntensity-BackGroundIntensity)/(max(max(SignalRecon))-min(min(SignalRecon)));

%% Hilbert Transform
Transform=hilbert(rfdata);
DemodTransform=abs(Transform);

figure(9)
imagesc(lat,depth,DemodTransform)
axis image
colormap hot
colorbar
xlabel('lateral direction (mm)')
ylabel('depth (mm)')
title('Data After Hilbert Transform')
caxis([0 4e-21])

SignalIntensity=mean(mean(DemodTransform(854:975,360:521)));
BackGroundIntensity=mean(mean(DemodTransform(1:300,1:300)));
ContrastHilbert=(SignalIntensity-BackGroundIntensity)/(max(max(DemodTransform))-min(min(DemodTransform)));

LogCompressed=log10(SignalRecon);
figure (10)
imagesc(lat,depth,LogCompressed)
axis image
colorbar
title('Log Compressed Quadrature Demodulation (Eye of Sauron)')
xlabel('lateral direction (mm)')
ylabel('depth (mm)')

SignalIntensity=mean(mean(LogCompressed(854:975,360:521)));
BackGroundIntensity=mean(mean(LogCompressed(1:300,1:300)));
ContrastLog=(SignalIntensity-BackGroundIntensity)/(max(max(LogCompressed))-min(min(LogCompressed)));

