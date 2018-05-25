import numpy as np
import matplotlib.pyplot as plt
def rrcosFilter(sampleLen, fs, fd=1, rolloff=0.25):
	Tsys = 1/fd
	Uf = (1+rolloff)/(2*Tsys)
	Df = (1-rolloff)/(2*Tsys)
	Hf_rrc = np.array([])
	fLPE=(np.array(range(1,sampleLen+1))/sampleLen-0.5*(sampleLen+1)/sampleLen)*fs
	for f in fLPE:
		if np.abs(f) < Df:
			Hf_rrc = np.append(Hf_rrc, 1)
		elif (np.abs(f) >= Df) and (np.abs(f) <= Uf):
			Hf_rrc = np.append(Hf_rrc, np.sqrt(0.5*(1+np.cos(np.pi*(2*Tsys*np.abs(f)-1+rolloff)/2/rolloff))))
		else:
			Hf_rrc = np.append(Hf_rrc, 0)
	Hf = np.append(Hf_rrc[-int(sampleLen/2):], Hf_rrc[:int(sampleLen/2)])
	return Hf


def main():
	rrcosFilter(8192,737.28,737.28/16, 0.25)
if __name__ == '__main__':
	main()