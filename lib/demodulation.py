import numpy as np

import matplotlib.pyplot as plt
def CarryIF(N):
    fc = N/4
    CarrierIF_I =np.array([])
    CarrierIF_Q =np.array([])
    for x in range(0,N):
        CarrierIF_I = np.append(CarrierIF_I, np.sin(2*np.pi*fc*x/N))
        CarrierIF_Q = np.append(CarrierIF_Q, np.cos(2*np.pi*fc*x/N))
    return CarrierIF_I, CarrierIF_Q

def calibaration(sI,sQ):
    #coarse calibaration
    index=np.where(sI == np.max(sI))
    phase_off1=np.arctan(sQ[index]/sI[index])-np.pi/4
    sI=((sI+1j*sQ)*np.exp(-1j*phase_off1)).real;
    sQ=((sI+1j*sQ)*np.exp(-1j*phase_off1)).imag;

    Mag_max=np.sqrt((sI[index]**2+sQ[index]**2))/np.sqrt(2)
    index11=np.where((sI>(Mag_max/2))&(sQ>(Mag_max/2)))
    index22=np.where((sI<(-Mag_max/2))&(sQ<(-Mag_max/2)))
    I11mean=np.mean(sI[index11])
    Q11mean=np.mean(sQ[index11])
    I22mean=np.mean(sI[index22])
    Q22mean=np.mean(sQ[index22])
    phase_off=np.arctan((Q11mean-Q22mean)/(I11mean-I22mean))-np.pi/4
    sI=((sI+1j*sQ)*np.exp(-1j*phase_off)).real;
    sQ=((sI+1j*sQ)*np.exp(-1j*phase_off)).imag;
    return sI,sQ


if __name__ == '__main__':
    CarrierIF_I, CarrierIF_Q = CarryIF(N=1024)
    w = np.hanning(len(CarrierIF_I))
    import FFTdBFs
    f,p = FFTdBFs.FFTdBFs(CarrierIF_I, win=w, fs = 737.28, ref=100)

    # plt.figure(3)
    # plt.title('square-wave after shaping filter(Freq Domain)')
    # plt.plot(f,p)
    # plt.show()
    a = np.array([1,3,5])
    b = np.array([2,3,4])
    print(a*b)
    print(np.mod(2.0,1))

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