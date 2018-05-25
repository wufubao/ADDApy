import numpy as np
import matplotlib.pyplot as plt
def carrier(fc, fs, N):
    t = np.array(range(0,N))/fs
    return np.exp(-1j*2*np.pi*fc*t)

def sample(sB1,sB2,offset = 8.5):
    fs = 737.28
    sample_dT=1/fs
    sampleLen = len(sB1)
    if np.mod(offset,1) == 0.0:
        offset = offset + 0.00000001
    fc = 60
    ssB1 = np.array([])
    ssB2 = np.array([])
    t_arry=sample_dT*offset+np.array(range(0,int(sampleLen/16)))*sample_dT*16
    for t in t_arry:
        ht = 2*fc*sample_dT * np.sin(2 * np.pi * fc * ( t-np.array(range(0,sampleLen))*sample_dT ))/(2 * np.pi * fc * ( t-np.array(range(0,sampleLen))*sample_dT ))
        ssB1 = np.append(ssB1, np.sum(ht*sB1))
        ssB2 = np.append(ssB2, np.sum(ht*sB2))
    return ssB1,ssB2