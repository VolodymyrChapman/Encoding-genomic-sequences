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
    #array = np.array(input[0])
    dftZ = np.fft.fft(cgr_z)
    PSZ = np.abs(dftZ)**2
    mean_signal = PSZ[0]
    signal = PSZ[1:] #skipping first signal 
    #freq = np.fft.fftfreq((psd.shape[0]) )
    return mean_signal, signal

def even_scale(power_signal, output_scale):
    
    real_scaled = (power_signal * len(power_signal))/output_scale
    floor_scaled = [math.floor(freq) for freq in real_scaled]
    scaled_signal = [] # - continue here: need to choose whether Q or R applicable and then 
    # calculate scaled signal
    for freq in real_scaled:
        scaled_signal.append 