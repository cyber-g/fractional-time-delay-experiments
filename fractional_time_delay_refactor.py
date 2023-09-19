# https://dsp.stackexchange.com/questions/60476/how-to-do-fft-fractional-time-delay-solved

import numpy as np
import matplotlib.pyplot as plt

f1 = 12.8
f2 = 22.6
samples = 2**13
fs = 10000.
tDelay = 0.1234/f1
tstart = 0.0
tend = samples/fs

# 0. Example waveform to demonstrate the time shift
timeList = np.arange(tstart, tend, 1/fs)
waveform = np.sin(2 * np.pi * f1 * timeList) + 1*np.sin(2 * np.pi * f2 * timeList)

waveform_delayed = np.sin(2 * np.pi * f1 * (timeList-tDelay)) + 1*np.sin(2 * np.pi * f2 * (timeList-tDelay))

# 1. Take the FFT
fftData = np.fft.fft(waveform)
fftData_delayed = np.fft.fft(waveform_delayed)

# 2. Construct the phase shift
samplePeriod = 1/fs
tDelayInSamples = tDelay / samplePeriod
N = fftData.shape[0]
k = np.linspace(0, N-1, N)
timeDelayPhaseShift = np.exp(((-2*np.pi*1j*k*tDelayInSamples)/(N)) + (tDelayInSamples*np.pi*1j))

# 3. Do the fftshift on the phase shift coefficients
timeDelayPhaseShift = np.fft.fftshift(timeDelayPhaseShift)

# 4. Multiply the fft data with the coefficients to apply the time shift
fftWithDelay = np.multiply(fftData, timeDelayPhaseShift)

plt.figure()
plt.plot(np.unwrap(np.angle(np.fft.fftshift(fftData))),label='original')
plt.plot(np.unwrap(np.angle(np.fft.fftshift(fftWithDelay))),label='delayed via fft')
plt.plot(np.unwrap(np.angle(np.fft.fftshift(fftData_delayed))),label='time delayed')
plt.grid()
plt.legend()
# plt.show()

plt.figure()
plt.plot(np.angle(np.fft.fftshift(fftData)),label='original')
plt.plot(np.angle(np.fft.fftshift(fftWithDelay)),label='delayed via fft')
plt.plot(np.angle(np.fft.fftshift(fftData_delayed)),label='time delayed')
plt.grid()
plt.legend()
# plt.show()

# 5. Do the IFFT
shiftedWaveform = np.fft.ifft(fftWithDelay)

print("\nThe sampling period is %f seconds" % samplePeriod)
print("The time delay is %f seconds" % tDelay)
print("The time delay in samples is %f samples" % tDelayInSamples)
print("The correction phase shift is %f pi" % (tDelayInSamples))

plt.figure()
plots = 1
plt.subplot(plots, 1, 1)
plt.plot(waveform,label='original')
plt.plot(shiftedWaveform,label='delayed via fft')
plt.plot(waveform_delayed,label='time delayed')
plt.grid()
plt.legend()
plt.show()