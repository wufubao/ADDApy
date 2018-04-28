# -*- coding: utf-8 -*-
#!/usr/bin/python3
# author: wufubao
from __future__ import division
from lib import FFTdBFs
import numpy as np
import matplotlib.pyplot as plt
def file_read(filname = '20MHzChA_M10dBm.txt'):
    '''
    Read the file and demap
    LMFS = 4211
    Usage:
        sig = file_read(f = filname)
        sigA = sig['ChA']
        sigB = sig['ChB']
    '''
    sampleChA = []
    sampleChB = []
    with open(filname,'r') as f:
        sig = f.read().splitlines()
    # samp_num = len(sig)
    for s in sig:
        # iTemp =  sig[x]
        Lane0_r = s[24:32]
        Lane1_r = s[16:24]
        Lane2_r = s[8:16]
        Lane3_r = s[0:8]
        sampleChA.append(Lane0_r[6:8]+Lane1_r[6:8])
        sampleChA.append(Lane0_r[4:6]+Lane1_r[4:6])
        sampleChA.append(Lane0_r[2:4]+Lane1_r[2:4])
        sampleChA.append(Lane0_r[0:2]+Lane1_r[0:2])
        sampleChB.append(Lane2_r[6:8]+Lane3_r[6:8])
        sampleChB.append(Lane2_r[4:6]+Lane3_r[4:6])
        sampleChB.append(Lane2_r[2:4]+Lane3_r[2:4])
        sampleChB.append(Lane2_r[0:2]+Lane3_r[0:2])
    return sampleChA, sampleChB

def Normalization(sig):
    '''
    Normolize the signal. full-scale is 1.7 vpp, 16 bit-res
    '''
    sample = np.array([])
    for s in sig:
        num = int(s, 16)
        if num > 0x8000:
            sample = np.append(sample, -(2**16-num))
        else:
            sample = np.append(sample, num)
    return sample

def main():
    sampleChA, sampleChB = file_read(filname = "AD_DATA/16QAM_0427.txt")
    waveA = Normalization(sampleChA)
    waveB = Normalization(sampleChB)
    needSave = True
    with open('AD_DATA/sigI.txt','w') as f:
        for x in waveA:
            f.write(str(x)+'\n')
    with open('AD_DATA/sigQ.txt','w') as f:
        for x in waveB:
            f.write(str(x)+'\n')
    w = np.hanning(len(waveB))
    fs=737.28
    f,p = FFTdBFs.FFTdBFs(waveA-waveB*1j,w, fs)
    plt.plot(f,p)
    plt.grid(True)
    plt.xlabel('Frequency [Hz] (Fs=%sMhz)'%fs)
    plt.ylabel('Amplitude [dBFS]')
    plt.xlim(-fs/2,fs/2)
    plt.show()
if __name__ == '__main__':
    main()

