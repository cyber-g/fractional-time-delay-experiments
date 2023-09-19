% A script to analyze fractional time delay effects on FFT. 

% Author: Germain Pham
% Date: Sept. 2023

clear all
close all

% Generate a simple signal : a sinc function
N         = 2^13;
SINC_BW   = 10; % in Hz
Fs        = 10000; % in Hz
t         = (-N/2:N/2-1)/Fs;
delay     = (1023/7)/Fs; % non-integer delay

% spectral_win = blackman(N)';
spectral_win = ones(1,N);


% Generate the signals
x = 2*SINC_BW*t;
y = sinc(x);

x_delayed = 2*SINC_BW*(t-delay);
y_delayed = sinc(x_delayed);

% Plot the signal
figure
plot(t,[y(:) y_delayed(:)]);
xlabel('Time (s)');
ylabel('Amplitude');
title('Sinc function');

% Compute the raw FFT
% Nfft = 2^(nextpow2(N)+1);
Nfft = N;
Y_raw = fft(y,Nfft)/sqrt(Nfft);
Y_delayed_raw = fft(y_delayed,Nfft)/sqrt(Nfft);

% Compute the windowed FFT
Y = fft(y.*spectral_win,Nfft)/sqrt(Nfft);
Y_delayed = fft(y_delayed.*spectral_win,Nfft)/sqrt(Nfft);

% Compute the frequency axis
f = (-Nfft/2:Nfft/2-1)*Fs/Nfft;

% Plot the raw spectrum
figure
plot(f,fftshift(10*log10(abs([Y_raw(:) Y_delayed_raw(:)]).^2)));
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Raw FFT of the sinc function');

figure
plot(f,fftshift(angle([Y_raw(:) Y_delayed_raw(:)])));
xlabel('Frequency (Hz)');
ylabel('Phase');
title('Phase of Raw FFT of the sinc function');


% Plot the spectrum
figure
plot(f,fftshift(10*log10(abs([Y(:) Y_delayed(:)]).^2)));
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Spectrum of the sinc function');

figure
plot(f,fftshift(angle([Y(:) Y_delayed(:)])));
xlabel('Frequency (Hz)');
ylabel('Phase');
title('Phase of the sinc function');



% Ratios of the FFTs
Z_raw = Y_delayed_raw./Y_raw;
Z = Y_delayed./Y;

% Plot the ratio
figure
plot(f,fftshift(10*log10(abs([Z_raw(:) Z(:)]).^2)));
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Ratio of the FFTs');
% ax1 = gca;
% % zoom on the useful part bandwidth
% figure
% copyobj(ax1, gcf);
% xlim([-SINC_BW SINC_BW]);

figure
plot(f,fftshift(angle([Z_raw(:) Z(:)])));
xlabel('Frequency (Hz)');
ylabel('Phase');
title('Phase of the ratio of the FFTs');
ax1 = gca;
% zoom on the useful part bandwidth
figure
copyobj(ax1, gcf);
xlim([-SINC_BW SINC_BW]);

% Bandwidth range for computation
compute_bw = 0.9*SINC_BW;

% Slope of phase:
bin_indx_useful_band = logical(abs(f) < compute_bw);
freq_useful_band = f(bin_indx_useful_band);
phase_slope = diff(angle(Z_raw(bin_indx_useful_band)));%./diff(2*pi*freq_useful_band);

plot(freq_useful_band(1:end-1),phase_slope)

avg_phase_slope = mean(phase_slope);
delay_estimate = -avg_phase_slope