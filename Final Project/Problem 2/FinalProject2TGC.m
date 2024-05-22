% FinalProject2TGC.m

%% Plotting the RF Data

clear;clf; close all
load ('rf_to_tgc.mat')
figure(1)
imagesc(lat, axial, RfData)
axis image
colormap gray
xlabel('lateral direction (m)')
ylabel('axial direction (m)')
title('Original Data')
[rows,cols]=size(RfData);

%% Hilbert Transform
Transform=hilbert(RfData);
DemodTransform=abs(Transform);

% Plotting the HilbertTransform

figure(2)
imagesc(lat, axial, DemodTransform)
axis image
colormap hot
colorbar
xlabel('lateral direction (m)')
ylabel('axial direction (m)')
title('Hilbert Transformed Image of Multiple Lesions')

SignalIntensity=mean(mean(DemodTransform(2250:2336,1:80)));
BackGroundIntensity=mean(mean(DemodTransform(631:693,161)));
ContrastDemod=(SignalIntensity-BackGroundIntensity)/(max(max(DemodTransform))-min(min(DemodTransform)))

% Finding the mean of every single row in the DemodTransform matrix
for k=1:rows
    Average(k)=mean(DemodTransform(k,:));
end

% Normalization Factor

NormalFactor=max(Average)./Average;

% Multiplying everything by the normalziation Factor

for i=1:rows
    DemodTransformTGC(i,:)=NormalFactor(i)*DemodTransform(i,:);
end

% Plotting the time gain compensated (normalized) data

figure(3)
imagesc(lat, axial, DemodTransformTGC)
axis image
colormap hot
colorbar
xlabel('lateral direction (m)')
ylabel('axial direction (m)')
title('Hilbert Transformed and Time Gain Compensated Image of Multiple Lesions')

SignalIntensity=mean(mean(DemodTransformTGC(637:701,161)));
BackGroundIntensity=mean(mean(DemodTransformTGC(2190:2336,1:41)));
ContrastTGC=(SignalIntensity-BackGroundIntensity)/(max(max(DemodTransformTGC))-min(min(DemodTransformTGC)))