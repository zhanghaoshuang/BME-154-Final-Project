% FinalProject1.m

clf; clear all; close all;
%% Plotting Impulse Responses

figure(1)
order = {'First'; 'Second'; 'Third'; 'Fourth'};
for index = 1:4
    F(index,:,:) = load(sprintf('TransImpResp%d.asc',index));
    time(index,:) = F(index,:,1);
    voltage(index,:) = F(index,:,2);
    subplot(2,2,index)
    plot(time(index,:),voltage(index,:))
    title([order{index} ' Impulse Response'])
    xlabel('time (s)')
    ylabel('voltage (V)')
end


%% Plotting and Calculating Power Spectrums

figure(2)
for index = 1:4
    T0(index) = time(index,end)-time(index,1);
    N = length(time(index,:));
    ts = T0(index)/N;
    fs = 1/ts;
    freq = fs*(-0.5:1/N:0.5-(1/N));
    freqplus = freq(length(freq)/2+1:end);
    FT_Voltage(index,:) = fft(voltage(index,:));
    FT_Shift_Voltage(index,:) = fftshift(FT_Voltage(index,:));
    Mag_FT_Voltage(index,:) = abs(FT_Shift_Voltage(index,:)/max(abs(FT_Shift_Voltage(index,:))));
    Power_Voltage(index,:) = 20*log10(Mag_FT_Voltage(index,:));
    Power_Voltage_Plus(index,:) = Power_Voltage(index,length(freq)/2+1:end); % Splits matrix in half
    subplot(2,2,index)
    plot(freq,Power_Voltage(index,:),'.-')
    xlabel('Frequency (Hz)'); ylabel('Power (dB)'); title(['Power Spectrum for ' order{index} ' Transducer'])
end





%% Finding -3dB

for index = 1:4
    maxindex(index)=find(max(Power_Voltage_Plus(index,:))==Power_Voltage_Plus(index,:)); % Finds where the maximum is 
    FirstTrans3dB1(index) = interp1(Power_Voltage_Plus(index,2:maxindex(index)),freqplus(2:maxindex(index)),-3,'linear'); % interpolates to find one -3dB intercept
    FirstTrans3dB2(index) = interp1(Power_Voltage_Plus(index, maxindex(index):end),freqplus(maxindex(index):end),-3,'linear'); % interpolates to find one -3dB intercept
    centerfreq(index) = geomean([FirstTrans3dB1(index),FirstTrans3dB2(index)]); % Uses the geometric mean of the -3dB points to calculate center freqquency
    FractionalBandwidth(index) = abs((FirstTrans3dB1(index)-FirstTrans3dB2(index)))/centerfreq(index);
end


%% Plotting Input
NumberofPeriods=22;

figure(3)

for index = 1:4
    timestep(index) = mean(diff(time(index,:)));
    period(index) = 1/centerfreq(index);
end

% Each time and input vector is of a different length; hence hardcoding is
% done
% This code makes sure to maintain the same time step in the transducer
t1 = linspace(0, NumberofPeriods*period(1), floor(NumberofPeriods*period(1)/timestep(1)));
t2 = linspace(0, NumberofPeriods*period(2), floor(NumberofPeriods*period(2)/timestep(2)));
t3 = linspace(0, NumberofPeriods*period(3), floor(NumberofPeriods*period(3)/timestep(3)));
t4 = linspace(0, NumberofPeriods*period(4), floor(NumberofPeriods*period(4)/timestep(4)));

Input1=2*sin(2*pi*centerfreq(1)*t1);
Input2=2*sin(2*pi*centerfreq(2)*t2);
Input3=2*sin(2*pi*centerfreq(3)*t3);
Input4=2*sin(2*pi*centerfreq(4)*t4);

for index = 1:4
    subplot(2,2,index)
    t = eval(sprintf('t%d',index));
    Input = eval(sprintf('Input%d',index));
    plot(t,Input)
    xlabel('time (s)')
    ylabel('Voltage (V)')
     title([order{index} ' Sinusoidal Input'])
end

%% Plotting the Output

figure(4)

for index = 1:4
    Input = eval(sprintf('Input%d',index));
    Output(index,:) = conv(voltage(index,:),Input,'same');
    subplot(2,2,index)
    plot(time(index,:),Output(index,:))
    xlabel('time (seconds)')
    ylabel('voltage (V)')
    title(['Output of ' order{index} ' Transducer'])
end

%% Determine When Steady State Occurs

% Hardcoding was used here because the findpeaks command would produce a
% different amount of peaks 

%This produces the maximum different between the peaks 

MaxSinePeaksDiff1=max(abs(diff(findpeaks(Output(1,:)))));
MaxSinePeaksDiff2=max(abs(diff(findpeaks(Output(2,:)))));
MaxSinePeaksDiff3=max(abs(diff(findpeaks(Output(3,:)))));
MaxSinePeaksDiff4=max(abs(diff(findpeaks(Output(4,:)))));

%This simply returns a 1 or 0 for if the output is at steady state or not
Steady1=MaxSinePeaksDiff1<0.01*max(Output(1,:))
Steady2=MaxSinePeaksDiff2<0.01*max(Output(2,:))
Steady3=MaxSinePeaksDiff3<0.01*max(Output(3,:))
Steady4=MaxSinePeaksDiff4<0.01*max(Output(4,:))

% This was hardcoded in to get a fit between cycles needed and center
% frequency

CyclesNeeded=[13 13 22 20];
figure (5)
plot (centerfreq,CyclesNeeded)
title('Plot of Cycles Needed vs, Cneter Frequency')
xlabel('Center Frequency (Hz)')
ylabel('Cycles Needed')