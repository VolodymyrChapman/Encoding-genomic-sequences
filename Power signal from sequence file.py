import numpy as np
import pandas as pd
import random as rand
import math
import matplotlib.pyplot as plt
from IPython.display import Markdown
import os
from os import listdir
from os.path import isfile, join

# Extract sequence from file
def getseqsample(file):
    pure_chunk = ''
    No_of_fastas = 0
    with open (file) as seq:
        for line in seq:
            if not line.startswith('>'):
                upperLine = line.upper().strip()
                pure_upper_line = ''.join([ s for s in upperLine if s in 'ATGC'])
                pure_chunk += pure_upper_line
        return pure_chunk 


#get CGR co-ordinates and Z only
def one_cold_chaos_coord(oligomer):
    ruleX = {'A': 0, 'C' : 0, 'G': 1, 'T' : 1}
    ruleY = {'A': 0, 'C' : 1, 'G': 1, 'T' : 0}
    X = [ruleX[i] for i in oligomer ] 
    X.insert(0,0.5)
    Y = [ruleY[i] for i in oligomer ] 
    Y.insert(0,0.5)
    
    for i in range(len(X[1:])):
        h = i + 1
        X[h] = (X[i] + X[h])/2
    for i in range(len(Y[1:])):
        h = i + 1
        Y[h] = (Y[i] + Y[h])/2
    Z = []
    for x,y in zip(X[1:],Y[1:]):
         Z.append(complex(x,y))
    return X, Y, Z

    # retrieve the power signal from the Z coordinates

def powersignal(cgr_z):
    dftZ = np.fft.fft(cgr_z)
    PSZ = np.abs(dftZ)**2
    mean_signal = PSZ[0]
    signal = PSZ[1:] #skipping first signal 
    return mean_signal, signal

# even scale power signal of one length to another length - method based on Yin and Yau, 2015 (J. Theoretical Biology)
def even_scale(signal, output_scale):
    scaled_signal = [signal[0]] # first term for tm and tn the same
    initialise = range(2,output_scale + 1)
    Q = [k*(len(signal)/output_scale) for k in initialise] 
    R = [math.floor(k) for k in Q]
    for k in range(len(Q)):
        if Q[k] == R[k]:
            Tm = signal[R[k]-1]  # correction for zero-based indexing
        else:
            Tm = signal[R[k]-1] + (Q[k] - R[k])*(signal[R[k]] - signal[R[k]-1]) # correction for zero-based indexing
        scaled_signal.append(Tm)    
    return scaled_signal
    
