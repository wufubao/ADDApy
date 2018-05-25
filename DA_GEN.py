# -*- coding: utf-8 -*-
#!/usr/bin/python3
# author: wufubao
from scipy import optimize, special
import matplotlib.pyplot as plt
import numpy as np
from lib import FFTdBFs,shapingFilter,modulation,demodulation

def nomi(sig,ref = 2**14):
    max_var = max(list(map(abs,sig)))
    nomi_sig = np.array([])
    for s in sig:
        if abs(s) == max_var:
            nomi_sig = np.append(nomi_sig, int((ref-1)*s/max_var))
        else:
            nomi_sig = np.append(nomi_sig, int((ref)*s/max_var))
    return nomi_sig

def ToHex(num):
    num = int(num)
    if num>=0:
        return hex(num)[2:]
    else:
        return hex(2**16+num)[2:]

def VH_Write(sigI,sigQ):
    Len = len(sigI)
    with open('DA_data/signali.vh', 'w') as f:
        for x in range(0,Len):
            f.write("assign i[%d]= 16'h%s;\n"%(x,ToHex(sigI[x])))
    Len = len(sigQ)
    with open('DA_data/signalq.vh', 'w') as f:
        for x in range(0,Len):
            f.write("assign q[%d]= 16'h%s;\n"%(x,ToHex(sigQ[x])))

def symbol_import(f = 'DA_data_0420_1/sig.txt'):
    with open(f,'r') as f:
        sigs = f.read().splitlines()
    sym = np.array([])
    for s in sigs:
        sym = np.append(sym,complex(s))
    return sym


def main():
    # BUG: THIS SCRIPT DOSE NOT WORKING
    symbolRate = 46.08 # symbol rate is 46.08 M Baud per second
    N = 16  # Sample insertion
    symbol_num = 768
    f_sample_rate = N*symbolRate # The DAC sample data rate is 737.28MHz

    need_symbol_import = 1  # import the data from the txt file
    if need_symbol_import:
        sB = symbol_import('DA_data/sig.txt')
    else:
        sB = modulation.QAMModem(symbol_num) # Or use 16QAM modulation, generate the baseband signal randomly
    # backup the baseband signal, so we can check the data in demodulation
    with open('DA_data/sig.txt','w') as f:
        for x in sB:
            comp = str(x)
            comp = comp.strip(')')
            comp = comp.strip('(')
            f.write(comp+'\n')
    sB = modulation.insertion(sB, N)    # N must be over or equal to 2 according to the Nyquist theorem
    sBsize = len(sB)
    sB = modulation.duplicate(sB,time=3)       # joint the baseband wave and to make sure DAC plays in continus [preFix,data,postFix]
    sB_I,sB_Q = modulation.getIQ(sB)    # TODO: add support for complex shapping
    t,rc = shapingFilter.rrcosfilter(int(symbol_num)*16, 0.25, 1/symbolRate, f_sample_rate) # generate a RRC FIR which Roll off factor is 0.25
    print(len(sB_I))
    sI = np.convolve(rc, sB_I,'same')   # Get the shaped wave
    sQ = np.convolve(rc, sB_Q,'same')
    # Convert to the DAC format(Verilog HDL files)
    sI_hlf = sI[sBsize:sBsize*2]    # [preFix,data,postFix] -> [data]
    sQ_hlf = sQ[sBsize:sBsize*2]
    sI_hlf= nomi(sI_hlf)        # convert to DAC data format
    sQ_hlf= nomi(sQ_hlf)
    VH_Write(sI_hlf,sQ_hlf)     # convert to verilog type

    # demod
    needDemod =0
    if needDemod:
        sI_hlf = modulation.duplicate(sI_hlf,time=3)
        sQ_hlf = modulation.duplicate(sQ_hlf,time=3)
        t,rc = shapingFilter.rrcosfilter(int(symbol_num), 0.25, 1/symbolRate, f_sample_rate)
        sI = np.convolve(rc, sI_hlf,'same')   # Get the shaped wave
        sQ = np.convolve(rc, sQ_hlf,'same')
        sI = sI[sBsize:sBsize*2]    # [preFix,data,postFix] -> [data]
        sQ = sQ[sBsize:sBsize*2]
        f = lambda t: np.abs(demodulation.sample(sI, sQ, offset = t)[0]+demodulation.sample(sI, sQ, offset = t)[0]*1j).var()
        offset = optimize.fminbound(f,0,15.99999)

        It,Qt = demodulation.sample(sI, sQ, offset = offset)
            


    needUpConvert = 0
    if needUpConvert:
        CarrierIF_I, CarrierIF_Q = demodulation.CarryIF(len(sI_hlf))
        sI_hlf = modulation.mixer(CarrierIF_I,sI_hlf)
        sI_hlf = modulation.mixer(CarrierIF_I,sI_hlf)

        sQ_hlf = modulation.mixer(CarrierIF_Q,sQ_hlf)
        sQ_hlf = modulation.mixer(CarrierIF_Q,sQ_hlf)
        sI_hlf = np.convolve(rc, sI_hlf,'same')
        sQ_hlf = np.convolve(rc, sQ_hlf,'same')


        # sI_hlf, sQ_hlf = demodulation.calibaration(sI_hlf, sQ_hlf)
        # sI_hlf = demodulation.sample(sI_hlf)
        # print(sI_hlf)

    isPlotNeed = 1
    if isPlotNeed:
        plt.figure(1)
        w = np.hanning(len(sB))
        f,p = FFTdBFs.FFTdBFs(sB, win=w, fs = f_sample_rate, ref=1)
        plt.title('The RRC')
        plt.plot(f,p)

        # w = np.hanning(len(sI_hlf))
        # f,p = FFTdBFs.FFTdBFs(sI_hlf, win=w, fs = f_sample_rate, ref=2**14)
        # plt.figure(2)
        # plt.title('square-wave after shaping filter(Freq Domain)')
        # plt.plot(f,p)

        # plt.figure(5)
        # plt.title('I')
        # plt.plot(sI)
        # plt.figure(6)
        # plt.title('I')
        # plt.plot(sI_hlf)
        # if needDemod:
        #     plt.figure(3)
        #     plt.title('Constellation symbols')
        #     plt.plot(sI,sQ,c='green',alpha=0.3)
        #     plt.scatter(It,Qt,c='red', alpha=1,s=10)

        #     w = np.hanning(len(sI))
        #     f,p = FFTdBFs.FFTdBFs(sI, win=w, fs = f_sample_rate, ref=2**14)
        #     plt.figure(4)
        #     plt.title('After modulation')
        #     plt.plot(f,p)



        plt.show()
if __name__ == '__main__':
    main()