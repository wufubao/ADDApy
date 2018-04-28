from numpy import random,log2,append,array
def bitarray2dec(in_bitarray):
    """
    Converts an input NumPy array of bits (0 and 1) to a decimal integer.
    Parameters
    ----------
    in_bitarray : 1D ndarray of ints
        Input NumPy array of bits.
    Returns
    -------
    number : int
        Integer representation of input bit array.
    """

    number = 0

    for i in range(len(in_bitarray)):
        number = number + in_bitarray[i]*pow(2, len(in_bitarray)-1-i)
    return number

class Modem:
    def __init__(self, binCode=None):
        self.__binCode = binCode
    @property
    def binCode(self):
        return self.__binCode
    def binCodeGen(self,Lenth):
        self.__binCode = random.randint(0, 2,Lenth)
    # def modulate(self, input_bits):


class QAMModem(Modem):
    __slots__ = ('__M','__symbols')
    def __init__(self, M=16, selfGen=True):
        self.__M = M
        self.__symbols = array([])
        self.num_bits_symbol = int(log2(self.__M))
        # self.symbol_mapping = arange(self.m)
        # self.num_bits_symbol
    @property
    def symbols(self):
        return self.__symbols
    def _constellation_symbol(self):
        for x in range(0,len(self.binCode), self.num_bits_symbol):
            self.__symbols = append(self.__symbols, \
                bitarray2dec(self.binCode[x:x+int(self.num_bits_symbol/2)])*2 - self.num_bits_symbol + 1 + \
                (bitarray2dec(self.binCode[x+int(self.num_bits_symbol/2):x+self.num_bits_symbol]) - self.num_bits_symbol + 1)*1j)
            # print(x)
    def getIQ(self):
        sI = np.array([])
        sQ = np.array([])
        for x in range(0, len(self.__symbols)):
            sI = np.append(sI, self.__symbols[x].real)
            sQ = np.append(sQ, self.__symbols[x].imag)
        return sI,sQ
    def insertion(self, N = 4):
    '''
    Nyquist theorem
    '''
        if N<2:
            raise ValueError('N should above or equal to 2')
        bb_s = np.array([])
        for x in range(0,len(self.__symbols)):
            for i in range(0,N):
                bb_s = np.append(bb_s, self.__symbols[x])
        self.__symbols = bb_s
    def duplicate(sB, time = 3):
        '''
        Joint
        '''
        sig = np.array([])
        for x in range(0,time):
            sig = np.append(sig,sB)
        return sig

if __name__ == '__main__':
    Q = QAMModem()
    Q.binCodeGen(32)
    Q._constellation_symbol()
    print(Q.binCode)
    print(Q.symbols)