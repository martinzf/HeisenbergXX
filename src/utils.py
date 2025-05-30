# utils.py
# Helper functions and classes for ease of use. Allow calling of Fortran code from Python.
from subprocess import call
import os
import numpy as np 
import pandas as pd 
from scipy import special 


CWD = os.path.join(os.getcwd(), 'src')


def write(df: pd.DataFrame, file: str):
    '''Writes data table to specified file.'''
    df.to_csv(f'{CWD}/{file}', sep=' ', header=False, index=False)

def read(file: str) -> np.ndarray:
    '''Reads data from specified file.'''
    return pd.read_csv(f'{CWD}/{file}', sep=' ', header=None, dtype=np.float64).to_numpy()


def error(a: np.ndarray, b: np.ndarray) -> tuple[float,float]:
    '''Statistics on two data sets.'''
    N = len(a)
    assert N == len(b)
    RMSE = np.sqrt(1/(N-1)*np.sum((a - b)**2))
    MEAN = np.mean(a)
    return RMSE, MEAN


def Renyi(nu: np.ndarray, a: float) -> np.ndarray:
    '''Rényi entanglement entropies for blocks of all sizes in a given spin chain.'''

    def s2(x: float, a: float) -> float:
        '''Binary (Rényi) entropy.'''
        if np.isclose(x, 0.) or np.isclose(x, 1.):
            return 0.
        if np.isclose(a, 1):
            return - x * np.log(x) - (1 - x) * np.log(1 - x)
        return np.log(x**a + (1 - x)**a) / (1 - a) 
    
    Sa = [np.sum([s2(pi, a) for pi in p]) for p in nu.T]
    return np.array(Sa)


class SpinChain:
    '''
    Class to represent an inhomogeneous Heisenberg XX spin chain with couplings J and external field B.
    
    Args:
        N (int): Number of particles
        J (np.ndarray): Length N-1 strictly positive array
        B (np.ndarray): Length N array
    '''
    def __init__(self, N: int, J: np.ndarray, B: np.ndarray):
        self.N = N
        self.J = J
        self.B = B
        # Save already calculated entropies for the chain
        self._S_cache = dict()

    def diagonalise(self):
        '''Calculate eigenvalues, eigenvectors and orthogonal polynomials'''
        df = pd.DataFrame([self.J, self.B])
        write(df, 'spin_chain.dat')
        call(['./spin_chain', 'spin_chain.dat',  f'{self.N}'], cwd=CWD)
        self.Phi = read('spin_chain_data/eigenvectors.dat')
        self.E, self.w, self.g = read('spin_chain_data/poly_data.dat')
        self.P = read('spin_chain_data/polynomials.dat')

    def filling(self, M: int):
        '''Set mode filling, calculate NxN correlation matrix'''
        self.M = M
        df = pd.DataFrame(self.Phi)
        write(df, 'spin_chain_data/eigenvectors.dat')
        call(['./correl',  f'{self.N}', f'{M}'], cwd=CWD)
        self.C = read('correl_data/correl.dat')
        self.nu = read('correl_data/eigvals.dat')

    def isothermal_length(self, x: np.ndarray, origin: str):
        '''
        Change to isothermal coordinates

        Args:
            x (np.ndarray): particle position array transformed into isothermal coordinates
            origin (str): l = [0, a, 2a, ..., l] (left), c = [-l/2, -l/2+a, ..., l/2] (centre)
        '''
        self.x = x
        self.origin = origin

        # Quantity used in entropy calculations f_N(L)
        if origin == 'l':
            self.f = 2 * x[-1] / np.pi * self.J * np.sin(np.pi * x[:-1]/x[-1])
        elif origin == 'c':
            self.f = 4 * x[-1] / np.pi * self.J * np.cos(np.pi/2 * x[:-1]/x[-1])
        else:
            raise ValueError('Origin can only be l or c.')
    
    def S(self, a: float=1) -> np.ndarray:
        '''Numerical entropy.'''
        if (a, self.M) not in self._S_cache:
            self._S_cache[(a, self.M)] = Renyi(self.nu, a)
        return self._S_cache[(a, self.M)]
    
    def S0(self, a: float=1) -> np.ndarray:
        '''Lowest order CFT entropy contribution.'''
        return 1/12*(1+1/a) * np.log(self.f)
    
    def S1(self, a: float=1) -> np.ndarray:
        '''Entropy parity oscillations.'''
        if a < 1:
            return np.zeros(self.N-1)
        if np.isclose(a, 1):
            mu = -1/4
        else:
            mu = 2**(1-2/a)/(1-a) * special.gamma(.5+.5/a) / special.gamma(.5-.5/a)
        return mu * -(-1)**np.arange(self.N-1) * self.f**(-1/a)

        
