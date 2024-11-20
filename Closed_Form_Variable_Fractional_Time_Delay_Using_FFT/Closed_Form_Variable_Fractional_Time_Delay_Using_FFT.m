% A script that implements 
% S. -C. Pei and Y. -C. Lai, "Closed Form Variable Fractional Time Delay Using FFT," in IEEE Signal Processing Letters, vol. 19, no. 5, pp. 299-302, May 2012, doi: 10.1109/LSP.2012.2191280.
% URL: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6171834&isnumber=6167335

% Author: Germain Pham
% Date: Sept. 2023

clear all
close all

phi             = @(t,p)    1/sqrt(2^p*factorial(p)*sqrt(pi))*exp(-t.^2/2).*hermiteH(p,t);
test_function   = @(t)      phi(t,5) + phi(t,12);


N   = 128;
t   = linspace(-7,7,N);
Ts  = diff(t(1:2));
y_0 = test_function(t);

% plot(t,y_0,'s-')


delay_vec = 0:0.5:10;
% for each delay, compute equation 4
for i = 1:length(delay_vec)
    delay       = delay_vec(i);
    H           = [ 1 ... 
                    exp(-1i*2*pi*delay/N*(1:((N/2)-1))) ...
                    cos(delay*pi) ...
                    exp(1i*2*pi*delay/N*(N-(((N/2)+1):(N-1))))]; % equation 4
    H_simpler   = ifftshift(exp(-1i*2*pi*delay/N*(-N/2:(N/2-1))));
    H_simpler(N/2+1) = H_simpler(N/2+1) - 1i*sin(delay*pi); % compensate Tarczynski error (Cain and Yardim, 1994)
    % H == H_simpler
    y_delayed_fft_based = ifft(fft(y_0).*H); % equation 5
    y_delayed_truth     = test_function(t-delay*Ts);
    RMSD(i) = sqrt(mean(abs(y_delayed_fft_based-y_delayed_truth).^2));

    % plot angles
    angle_truth       = angle(fft(y_delayed_truth)./fft(y_0));
    angle_fft_based   = angle(H);
    % figure
    % plot([angle_fft_based(:) angle_truth(:)])
    % legend('fft-based','truth')
    % plot slopes
    slope_truth       = diff(unwrap(angle_truth));
    slope_fft_based   = diff(unwrap(angle_fft_based));
    figure
    plot([slope_fft_based(:) slope_truth(:)])
    legend('fft-based','truth')
end

figure
plot(delay_vec,RMSD)

figure
plot([y_delayed_fft_based(:) y_delayed_truth(:)])
legend('fft-based','truth')

figure
plot(10*log10(abs(fft(y_0))))
hold on
plot(10*log10(abs(fft(y_0.*blackman(N).'))))
hold off
legend('Window:rectangular','Window:blackman')