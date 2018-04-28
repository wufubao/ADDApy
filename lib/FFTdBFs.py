# -*- coding: utf-8 -*-
#!/usr/bin/python3
# author: wufubao
import numpy as np

def FFTdBFs(wave , win = None, fs=245.76, ref = 0x8000):
    '''
    Calculate spectrum in dBFs scale
    FFT sample should be large enough.
    Args:
        x: input signal
        fs: sampling frequency
        win: vector containing window samples (same length as x).
             If not provided, then rectangular window is used by default.
        ref: reference value used for dBFS scale. 32768 for int16 and 1 for float
    Returns:
		freq: frequency vector
        yfp: spectrum in dBFs scale
    '''
    fft_size = len(wave)
    if win is None:
        win = np.ones(fft_size)
    if fft_size != len(win):
            raise ValueError('Signal and window must be of the same length')
    wave = wave * win
    yf = np.fft.fft(wave)/np.sum(win)
    yf = np.fft.fftshift(yf)
    freq = np.linspace(-fs/2,fs/2, fft_size)
    yfp = 20*np.log10(abs(yf)/ref)
    return freq, yfp