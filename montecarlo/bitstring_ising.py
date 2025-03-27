import numpy as np
import math    
import copy as cp
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import random
random.seed(2)         

class BitString:
    """
    Simple class to implement a config of bits
    """
    def __init__(self, N):
        self.N = N
        self.config = np.zeros(N, dtype=int) 

    def __repr__(self):
        out = ""
        for i in self.config:
            out += str(i)
        return out

    def __eq__(self, other):        
        return all(self.config == other.config)
    
    def __len__(self):
        return len(self.config)

    def on(self):
        """
        Return number of bits that are on
        """
        count = 0
        for bit in self.config:
            if bit == 1:
                count += 1
        return count


    def off(self):
        """
        Return number of bits that are off
        """
        count = 0
        for bit in self.config:
            if bit == 0:
                count += 1
        return count

    def flip_site(self,i):
        """
        Flip the bit at site i
        """
        self.config[i] = 1 - self.config[i]

    
    def integer(self):
        """
        Return the decimal integer corresponding to BitString
        """
        value = 0
        i= len(self.config)
        for bit in self.config:
            if bit == 1:
                value += pow(2, i - 1)
            i -= 1
        return value
 

    def set_config(self, s:list[int]):
        """
        Set the config from a list of integers
        """
        i = 0
        for integer in s:
            self.config[i] = integer
            i += 1

    def set_integer_config(self, dec:int):
        """
        convert a decimal integer to binary
    
        Parameters
        ----------
        dec    : int
            input integer
            
        Returns
        -------
        Bitconfig
        """
        count = 0
        self.config = np.array([], dtype=int)
        while(dec > 0 or count < self.N):
            result = divmod(dec, 2)
            self.config = np.append(self.config,result[1])
            dec = result[0]
            count += 1
        self.config = self.config[::-1]
    
    def get(self, i):
        return self.config[i]
    
class IsingHamiltonian:

    
    def __init__(self, G):
        self.G = G
        self.mu = np.zeros(G.number_of_nodes())

    def set_mu(self, mu:np.array):
        self.mu = mu

    def energy(self, bs: BitString):
        energy = 0
        total_e = 0
        A = nx.adjacency_matrix(self.G).todense() #get hamiltonian

        for i in range(bs.__len__()):
            for j in range(i+1, bs.__len__()):
                if bs.get(i) == bs.get(j):
                    energy = 1
                else:
                    energy = -1

                energy *= A[i][j] #multiply the energy by the weight between each node
                total_e += energy
        #apply mu component in the sum
        for i in range(bs.__len__()):

            if bs.get(i) == 1:
                spin = 1
            else:
                spin = -1
            total_e += self.mu[i]*spin

        return total_e
    
    def get_probability(self, bs:BitString, T: float):
        beta = (1/(T)) #since k is just 1 in this approximation
        z = 0
        bs_copy = BitString(bs.__len__())

        for i in range(pow(2,bs_copy.__len__())):
            bs_copy.set_integer_config(i)
            energy_ = self.energy(bs_copy)
            z += pow(math.e, -beta*energy_)
        
        return (pow(math.e, -beta*self.energy(bs)))/z

    def compute_average_values(self, T: float):

        E  = 0.0
        M  = 0.0
        HC = 0.0
        MS = 0.0

        bs = BitString(self.G.number_of_nodes())

        # Write your function here!
        for i in range(pow(2,bs.__len__())):
            bs.set_integer_config(i)
            prob = self.get_probability(bs, T)
            E +=  prob*self.energy(bs)
            M += (bs.on() - bs.off()) * prob
            HC += prob*self.energy(bs)*self.energy(bs)
            MS += pow(bs.on() - bs.off(), 2) * prob
        HC -= (E*E)
        HC = HC/(T*T)
        MS = MS/T

        return E, M, HC, MS
    
    