%mmodecorrv.m

% BME 154 Final Project
% Question 3

%% Initialization
clear all;
close all;

%% Load data
file = 'C';
a = load(['mmode', lower(file), '.mat']); % load filename

% % % figure
% % % subplot(1,3,1)
% % % imagesc(a.mmode)
% % % colormap gray
% % % title('Dataset A')
% % % 
% % % subplot(1,3,2)
% % % imagesc(b.mmode)
% % % colormap gray
% % % title('Dataset B')
% % % 
% % % subplot(1,3,3)
% % % imagesc(c.mmode)
% % % colormap gray
% % % title('Dataset C')
% % % 
% % % print -dpdf 3velocityimages

% find the axial and time step lengths using linear regression
dx = dot(polyfit(1:length(a.axial), a.axial, 1), [1, 0]);
dt = dot(polyfit(1:length(a.T), a.T, 1), [1, 0]);

%% Calculating velocities
corrdata = xcorr(a.mmode(:,1), a.mmode(:,1));
Axcorrmat = zeros(size(corrdata)); % cross-correlation matrix
Acorrpeakloc = zeros(size(a.T)); % vector of peaks

for k=1:length(a.T)
    Axcorrmat(:, k) = xcorr(a.mmode(:,1), a.mmode(:, k));
    Axcorrpeakloc(k) = find(Axcorrmat(:, k) == max(Axcorrmat(:, k))) - length(a.mmode);
end

v = dot(polyfit(1:length(Axcorrpeakloc), Axcorrpeakloc, 1), [1, 0])*dx/dt

%% Produce figures
x = (1:length(a.axial))*dx;

% % comparison of two RF lines
% figure
% subplot(2,1,1)
% plot(x, a.mmode(:, 1), 'k-');
% title(['(Dataset ', file, ') Single RF line at t = 0s'])
% ylabel('Magnitude')
% 
% subplot(2,1,2)
% plot(x, a.mmode(:,32), 'k-');
% title(['(Dataset ', file, ') Single RF line at t = 0.032s'])
% xlabel('Depth (m)')
% ylabel('Magnitude')
% 
% print -dpng 3RFlinecompB
% 
% % comparison of autocorrelation and xcorrlation
% x2 = (1:length(Axcorrmat(:,1)))*dx;
% shift0 = Axcorrmat(:, 1); 
% shift0max = find(shift0 == max(shift0));
% shift32 = Axcorrmat(:, 32);
% shift32max = find(shift32 == max(shift32));
% 
% figure
% subplot(2,1,1)
% plot(x2, shift0,'k-', [1,1]*shift0max*dx, [-2e-36, 2e-36], 'r-')
% title(['(Dataset ', file, ') Autocorrelation of the RF line at t = 0s'])
% xlabel('Depth (m)')
% ylabel('Magnitude')
% 
% subplot(2,1,2)
% plot(x2, shift32, 'k-', [1,1]*shift32max*dx, [-2e-36, 2e-36], 'r-')
% title(['(Dataset ', file, ') Cross-correlation of RF lines at t = 0 and 0.032s'])
% xlabel('Depth (m)')
% ylabel('Magnitude')
% 
% print -dpng 3corrcompB
% 
% % finding velocity
% figure
% plot((1:length(Axcorrpeakloc))*dt, (Axcorrpeakloc)*dx, 'k-')
% title(['(Dataset ', file, ') Correlation peak shifts through time'])
% xlabel('Time (s)')
% ylabel('Depth(m)')
% annotation('textbox',...
%     [0.516071428571428 0.75852380952381 0.0767857142857143 0.0666666666666667],...
%     'String',{['v = ', num2str(v), ' m/s']}, 'LineStyle', 'none');
% 
% print -dpng velocityplotB
