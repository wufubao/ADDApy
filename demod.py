# -*- coding: utf-8 -*-
#!/usr/bin/python3
# author: wufubao
from __future__ import division
from lib import FFTdBFs,Filter,demodulation
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
    '''
    FUNC calibaration and FUNC evm is not open source
    '''
    fs=737.28
    symbolRate = 46.08
    OverSampling = 16
    fc = fs/4
    ridLeng = 3
    # TODO: You can replace you data to waveA and waveB for I and Q
    sampleChA, sampleChB = file_read(filname = "AD_data/DMA.txt")
    waveA = Normalization(sampleChA)
    waveB = Normalization(sampleChB)
    needSave = True
    with open('AD_data/sigI.txt','w') as f:
        for x in waveA:
            f.write(str(x)+'\n')
    with open('AD_data/sigQ.txt','w') as f:
        for x in waveB:
            f.write(str(x)+'\n')
    

    If = demodulation.carrier(fc, fs, len(waveB))
    I_outR=((waveA-1j*waveB)*If).real
    Q_outR=((waveA-1j*waveB)*If).imag

    Hf_rrc = Filter.rrcosFilter(len(I_outR), fs, symbolRate, rolloff=0.25)
    ssIFix_r=np.fft.ifft(np.fft.fft(I_outR) * Hf_rrc).real
    ssQFix_r=np.fft.ifft(np.fft.fft(Q_outR) * Hf_rrc).real
    offset = 1.3
    needOpt = False
    if needOpt:
        from scipy import optimize, special
        f = lambda t: demodulation.evm(demodulation.sample(ssIFix_r, ssQFix_r, offset = t)[0],\
            demodulation.sample(ssIFix_r, ssQFix_r, offset = t)[1], False)
        for x in range(0,OverSampling):
            offset = optimize.fminbound(f,x,x+1)
            if f(offset) < 0.1:
                break
    It,Qt = demodulation.sample(ssIFix_r, ssQFix_r, offset = offset)
    ssI_r = ssIFix_r[ridLeng*OverSampling:-ridLeng*OverSampling]
    ssQ_r = ssQFix_r[ridLeng*OverSampling:-ridLeng*OverSampling]
    It = It[ridLeng:-ridLeng]
    Qt = Qt[ridLeng:-ridLeng]
    needCal = False
    if needCal:
        ssI_r, ssQ_r, It,Qt = demodulation.calibaration(ssIFix_r, ssQFix_r, It,Qt)
        evm = demodulation.evm(It, Qt)*100
        print("EVM is %f%%"%evm)
        plt.figure(2)
        plt.subplot(211)
        X = np.arange(len(It))+1
        plt.bar(X,It)
        plt.subplot(212)
        plt.bar(X,Qt)

    plt.figure(1)
    plt.plot(ssI_r, ssQ_r,c='green',alpha=0.3)
    plt.scatter(It,Qt,c='red', alpha=1,s=10)
    
    # w = np.hanning(len(waveB))
    # f,p = FFTdBFs.FFTdBFs(ssIFix_r*1j-ssQFix_r,w, fs)
    # plt.figure(2)
    # plt.plot(f,p)
    # plt.grid(True)
    # plt.xlabel('Frequency [Hz] (Fs=%sMhz)'%fs)
    # plt.ylabel('Amplitude [dBFS]')
    # plt.xlim(-fs/2,fs/2)
    plt.show()
if __name__ == '__main__':
    main()

