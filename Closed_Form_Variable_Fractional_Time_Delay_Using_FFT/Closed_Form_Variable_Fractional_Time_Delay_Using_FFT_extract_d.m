% A script that implements 
% S. -C. Pei and Y. -C. Lai, "Closed Form Variable Fractional Time Delay Using FFT," in IEEE Signal Processing Letters, vol. 19, no. 5, pp. 299-302, May 2012, doi: 10.1109/LSP.2012.2191280.
% URL: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6171834&isnumber=6167335

% Author: Germain Pham
% Date: Sept. 2023

clear all
close all

phi             = @(t,p)    1/sqrt(2^p*factorial(p)*sqrt(pi))*exp(-t.^2/2).*hermiteH(p,t);
test_function   = @(t)      phi(t,5) + phi(t,12);


N   = 256;
t   = linspace(-10,10,N);
Ts  = diff(t(1:2));
Fs  = 1/Ts;
y_0 = test_function(t);

% plot(t,y_0,'s-')

% Bandwidth of the signal : 0.8 Hz (when Fs=9Hz)

delay_vec = 0:0.5:50;
% for each delay, compute equation 4
for i = 1:length(delay_vec)
    delay             = delay_vec(i);
    y_delayed_truth   = test_function(t-delay*Ts);
    
    % get angles
    angle_truth       = angle(fft(y_delayed_truth)./fft(y_0));
    slope_truth       = diff(unwrap(angle_truth));
    
    % get the slope of the phase in the bandwith of the signal
    slope_avg(i) = mean(slope_truth(1:round(0.8/Fs*N)));
    
end

% transform the slope into a delay
delay_est = -slope_avg/(2*pi)*N;

% compute the error
error = abs(delay_est-delay_vec).^2;
plot(delay_vec,error,'s-')
xlabel('Delay (samples wrt Ts)')
ylabel('Error')
title('Error of the delay estimation')

% Create inset
axes1 = axes('Parent',gcf,...
    'Position',[0.264285714285714 0.387019230769231 0.425 0.483173076923077]);

% plot the time domain signals
plot(t,y_0,'s-')
hold on
plot(t,test_function(t-delay_vec(delay_vec==max(delay_vec))*Ts),'s-')
hold off
xlabel('Time (s)','BackgroundColor',[1 1 1])
ylabel('Amplitude','BackgroundColor',[1 1 1])
legend('Original signal','Delayed signal','Location','northoutside')

% Create arrow
annotation(gcf,'arrow',[0.707142857142857 0.860714285714286],...
    [0.696115384615385 0.884615384615385]);