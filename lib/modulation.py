# -*- coding: utf-8 -*-
#!/usr/bin/python3
# author: wufubao
import numpy as np

def QAMModem(size = 1024, M=16):
    """ Creates a Quadrature Amplitude Modulation (QAM) Modem object."""
    if M!= 4 and M!=16 and M!=64 and M!=256:
        raise ValueError('M is not supported!')
    else:
        M = int(np.sqrt(M))
        return np.random.randint(0, M,size)*2 - M + 1 + (np.random.randint(0, M,size)*2 - M + 1)*1j

def insertion(sB, N = 4):
    '''
    Nyquist theorem
    '''
    if N<2:
        raise ValueError('N should above or equal to 2')
    bb_s = np.array([])
    for x in range(0,len(sB)):
        for i in range(0,N):
            bb_s = np.append(bb_s, sB[x])
    return bb_s

def duplicate(sB, time = 3):
    '''
    Joint
    '''
    sig = np.array([])
    for x in range(0,time):
        sig = np.append(sig,sB)
    return sig

def getIQ(sB):
    return sB.real, sB.imag

def test():
    for x in QAMModem(M=256):
        print(x)

def mixer(IFc,sB):
    ssB = np.array([])
    for x in range(0,len(sB)):
        ssB = np.append(ssB, IFc[x]*sB[x])
    return ssB
if __name__ == '__main__':
    test()

