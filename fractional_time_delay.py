# https://dsp.stackexchange.com/questions/60476/how-to-do-fft-fractional-time-delay-solved

import numpy as np
import matplotlib.pyplot as plt

f1 = 12.8
f2 = 22.6
samples = 1024
tDelay = .00938
tstart = 0.0
tend = 1.0

# 0. Example waveform to demonstrate the time shift
timeList = np.linspace(tstart, tend, samples)
waveform = np.sin(2 * np.pi * f1 * timeList) + 1*np.sin(2 * np.pi * f2 * timeList)

# 1. Take the FFT
fftData = np.fft.fft(waveform)

# 2. Construct the phase shift
samplePeriod = (tend - tstart) / (samples)
tDelayInSamples = tDelay / samplePeriod
N = fftData.shape[0]
k = np.linspace(0, N-1, N)
timeDelayPhaseShift = np.exp(((-2*np.pi*1j*k*tDelayInSamples)/(N)) + (tDelayInSamples*np.pi*1j))

# 3. Do the fftshift on the phase shift coefficients
timeDelayPhaseShift = np.fft.fftshift(timeDelayPhaseShift)

# 4. Multiply the fft data with the coefficients to apply the time shift
fftWithDelay = np.multiply(fftData, timeDelayPhaseShift)

# 5. Do the IFFT
shiftedWaveform = np.fft.ifft(fftWithDelay)

print("\nThe sampling period is %f seconds" % samplePeriod)
print("The time delay is %f seconds" % tDelay)
print("The time delay in samples is %f samples" % tDelayInSamples)
print("The correction phase shift is %f pi" % (tDelayInSamples))

plots = 1
plt.subplot(plots, 1, 1)
plt.plot(waveform)
plt.plot(shiftedWaveform)
plt.show()