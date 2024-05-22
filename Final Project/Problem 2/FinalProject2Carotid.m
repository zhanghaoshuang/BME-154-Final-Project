%FinalProject2Carotid.m

clear;clf;close all

% Producing the original image

load('rf_carotid.mat')
figure(1)
imagesc(lat,axial,double(rfdata))
xlabel('lateral direction (m)')
ylabel('axial direction (m)')
title('Original Carotid Artery Image')
colormap hot
axis image

% Hilbert Tranform
Transform=hilbert(double(rfdata));
DemodTransform=abs(Transform);

% Plotting the Hilbert Transform of the Image

figure(2)
imagesc(lat,axial,DemodTransform)
axis image
colormap gray
colorbar
xlabel('lateral direction (m)')
ylabel('axial direction (m)')
title('Carotid Artery Image with Only Hilbert Transform')
SignalIntensity=mean(mean(DemodTransform(708:817,116)));
BackGroundIntensity=mean(mean(DemodTransform(661:682,52:70)));
ContrastHilbert=(SignalIntensity-BackGroundIntensity)/(max(max(DemodTransform))-min(min(DemodTransform)))
% Plotting the log compressed data of Hilbert Transform

figure(3)
DemodTransformLog=log10(DemodTransform);
imagesc(lat,axial,DemodTransformLog)
axis image
colormap gray
colorbar
xlabel('lateral direction (m)')
ylabel('axial direction (m)')
title('Carotid Artery Image with Hilbert and Log Compression')
SignalIntensity=mean(mean(DemodTransformLog(708:817,116)));
BackGroundIntensity=mean(mean(DemodTransformLog(661:682,52:70)));
ContrastHilbertLog=(SignalIntensity-BackGroundIntensity)/(max(max(DemodTransformLog))-min(min(DemodTransformLog)))
% Normalizing the Hilbert Transformed Data

DemodTransformNorm = DemodTransform./max(max(DemodTransform));

%Plotting the histogram of the normalized data

figure(4)
hist (DemodTransformNorm)
xlabel('Intensity')
ylabel('Frequency of Occurances')
title('Histogram for Hilbert Transformed Carotid Image')

% Plotting Image Adjusted Verstion of the Normalized Data

figure (5)
ImAdjust=imadjust(DemodTransformNorm,[0, 0.1]);
imagesc(lat,axial,ImAdjust)
axis image
colorbar
colormap gray
xlabel('lateral direction (m)')
ylabel('axial direction (m)')
title('Carotid Artery Image with Hilbert and Imadjust')
SignalIntensity=mean(mean(ImAdjust(708:817,116)));
BackGroundIntensity=mean(mean(ImAdjust(661:682,52:70)));
ContrastImadjust=(SignalIntensity-BackGroundIntensity)/(max(max(ImAdjust))-min(min(ImAdjust)))

% Taking the log compressed of above

figure (6)
ImAdjustLog=log10(ImAdjust);
imagesc(lat, axial, ImAdjustLog)
axis image
colorbar
colormap gray
xlabel('lateral direction (m)')
ylabel('axial direction (m)')
title('Carotid Artery Image with Hilbert Imadjust and Log Compression')
SignalIntensity=mean(mean(ImAdjustLog(708:817,116)));
BackGroundIntensity=mean(mean(ImAdjustLog(661:682,52:70)));
ContrastImadjustLog=(SignalIntensity-BackGroundIntensity)/(max(max(ImAdjustLog))-min(min(ImAdjustLog)))