# https://dsp.stackexchange.com/questions/509/what-effect-does-a-delay-in-the-time-domain-have-in-the-frequency-domain

import numpy as np
import scipy.fft as sf
import matplotlib.pyplot as plt

# Generate a signal delayed or not
def signal_gen(length, delay=0, aliasing=False):
    indices = np.arange(-delay, -delay+length)
    if aliasing: indices = indices % length
    values = np.linspace(1, 0, length) ** 2 + 0.2
    xn = np.where((indices<0) | (indices>=length), 0, values[indices])
    return xn

# Plot signal
def plot(xn, Xk, title='Signal', unwrap=False):
    # Plot signal
    fig, axes = plt.subplots(ncols=3, squeeze=True, figsize=(8, 2.5))
    ax = axes[0]
    ax.set_title(title)
    ax.stem(xn)

    # Plot DFT magnitude
    ax = axes[1]
    ax.set_title('DFT magnitude')
    ax.stem(abs(Xk))

    # Plot DFT angle
    ax = axes[2]
    ax.set_title('DFT angle')
    angles = np.angle(Xk)
    if unwrap: angles = np.unwrap(angles)
    ax.stem(angles)

N = 9 # number of samples to create
nd = 4 # time-delay in samples

# Original and time-shifted signals
xn = signal_gen(N, 0)
xn_s = signal_gen(N, nd, aliasing=True)

# Plot
Xk, Xk_s = [sf.fft(signal) for signal in [xn, xn_s]]
plot(xn, Xk, title='Unshifted signal')
plot(xn_s, Xk_s, title='Time-shifted signal')

# Phase differences between DFT of shifted and unshifted signals
phase_wrapped = np.angle(Xk_s) - np.angle(Xk)
phase = np.unwrap(phase_wrapped)

# Print wrapped and unwrapped phase differences and increments
print(np.array2string(phase_wrapped, formatter={'float_kind': lambda x: "%.2f" % x}))
print(np.array2string(phase, formatter={'float_kind': lambda x: "%.2f" % x}))
print(np.array2string(np.diff(phase), formatter={'float_kind': lambda x: "%.2f" % x}))

# Compute time-delay
pd = phase[1]-phase[0]
nd2 = -N * pd / (2*np.pi)
print(f'nd = {nd2:.2f}')

plt.show()