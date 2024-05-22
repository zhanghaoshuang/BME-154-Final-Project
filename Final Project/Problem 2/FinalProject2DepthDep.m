% FinalProject2DepthDep.m

%% Plotting the RF Data

clear all; clf; close all
load('rf.mat')
figure(1)
imagesc(rfdata)
axis image
colormap gray
xlabel('lateral direction (mm)')
ylabel('depth (mm)')
title('Original Radio Frequency Data')


% The first one was used to split rfdata into smaller pieces to try to
% eliminate the lesion but this didn't really work
 
rfdata1=rfdata;
% Finding the dimensions of the matrix

[NumRows1,NumCols1]=size(rfdata1);
[NumRows,NumCols]=size(rfdata);

% Calculating time
t=2*depth/1540/1000;

%Spltting the whole matrix into columns into to find depth dependent center
%frequency

cols=3;
rows=floor(NumRows/cols);

rfdatashort=zeros(rows,3*NumCols);
tshort=zeros(rows,cols);

% Goes through and splits up the rfdata into the different number of
% columns (chunks)

for Index=0:NumCols1-1
    for i=0:cols-1
        rfdatashort(:,(i+1)+cols*Index)=rfdata1(rows*i+1:rows*(i+1),Index+1);
        tshort(:,i+1)=t(rows*i+1:rows*(i+1));
    end
end

%Fourier Transform of the split rfdata

T0=tshort(rows,1)-tshort(1,1);
N=rows;
ts=T0/N;
fs=1/ts;
freq= fs* (0 : 1/N: 1-(1/N));

rfdatashortFT=abs(fft(rfdatashort));
rfdatashortPower=20*log(rfdatashortFT);

figure(2)
plot (freq,rfdatashortFT)
title('Amplitude Response of Split RF Data')
xlabel('frequency (Hz)')
ylabel('amplitude')

%Finding the different center frequencies

for k=1:NumCols1*cols
    Index(k)=find(max(rfdatashortPower(1:round(length(freq)/2),k))==...
        rfdatashortPower(1:round(length(freq)/2),k));
    CenterFreq(k)=freq(Index(k));
end

%Averaging the center frequencies for each column but across all the rows
%in the column (chunks)

for Index2=1:cols
    AverageFreq(Index2)=mean(CenterFreq(Index2:cols:NumCols1*cols));
end

figure(3)
plot(AverageFreq)
xlabel('Number of Chunk')
ylabel('Center Frequency')
title('The Trend of Center Frequency with Respect to Depth')



%% Quadrature Demodulation
n=cols;
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

figure(4)
plot(freq,PowerNorm)
title('Power Spectrum of RF Data Columns (A-Lines)')
xlabel('frequency (Hz)')
ylabel('power (dB')

%Choose Frequency
t=(depth*2/1000/1540)';
SignalCos=zeros(NumRows,NumCols);
SignalSin=zeros(NumRows,NumCols);

%Multiply by Sine and Cosine
cols=3;
for i=0:cols-1
    for k=1:NumCols
    SignalCos(floor(i*NumRows/3)+1:floor((i+1)*NumRows/3),k)=...
        cos(2*pi*AverageFreq(i+1)*t(floor(i*NumRows/3)+1:floor((i+1)*NumRows/3)))...
        .*rfdata(floor(i*NumRows/3)+1:floor((i+1)*NumRows/3),k);
    SignalSin(floor(i*NumRows/3)+1:floor((i+1)*NumRows/3),k)=...
        sin(2*pi*AverageFreq(i+1)*t(floor(i*NumRows/3)+1:floor((i+1)*NumRows/3)))...
        .*rfdata(floor(i*NumRows/3)+1:floor((i+1)*NumRows/3),k);
    end 
end

[row, cols]= size(rfdata);    
%FT
SignalCosFT=fft(SignalCos);
SignalSinFT=fft(SignalSin);


figure(5)
plot(freq,abs(SignalCosFT))
title('Amplitude Response of RF Data Multiplied by Cosine')
xlabel('frequency (Hz)')
ylabel('Amplitide')

figure(6)
plot(freq,abs(SignalSinFT))
plot(freq,abs(SignalSinFT))
title('Amplitude Response of RF Data Multiplied by Sine')
xlabel('frequency (Hz)')
ylabel('Amplitide')

% Multiply by Rect/Filter
SignalCosFTFiltered=zeros(row,cols);
SignalSinFTFiltered=zeros(row,cols);
for i=1:cols
    SignalCosFTFiltered(:,i)=SignalCosFT(:,i).*(abs(freq)<1e7)';
    SignalSinFTFiltered(:,i)=SignalSinFT(:,i).*(abs(freq)<1e7)';
end

figure(7)
plot(freq,abs(SignalCosFTFiltered))
title('Amplitude Response of Filtered RF Data Multiplied by Sine')
xlabel('frequency (Hz)')
ylabel('Amplitide')

figure(8)
plot(freq,abs(SignalSinFTFiltered))
plot(freq,abs(SignalSinFTFiltered))
title('Amplitude Response of Filtered RF Data Multiplied by Sine')
xlabel('frequency (Hz)')
ylabel('Amplitide')

% Reconstruct/InverseFT
RealSignal=real(ifft(SignalCosFTFiltered));
ImaginarySignal=real(ifft(SignalSinFTFiltered));

figure(9)
SignalRecon=2*sqrt(RealSignal.^2+ImaginarySignal.^2);
imagesc(lat, depth, SignalRecon)
xlabel('lateral direction (mm)')
ylabel('depth (mm)')
title('Image with Quadrature Demodulation including Depth Dependent Center Frequency')
axis image
colormap hot
colorbar
caxis([0 4e-21])

%% Contrast

SignalIntensity=mean(mean(SignalRecon(854:975,360:521)))
BackGroundIntensity=mean(mean(SignalRecon(1:300,1:300)))
ContrastDepth=(SignalIntensity-BackGroundIntensity)/(max(max(SignalRecon))-min(min(SignalRecon)))