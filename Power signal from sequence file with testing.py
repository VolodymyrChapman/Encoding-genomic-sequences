import numpy as np
import pandas as pd
import random as rand
import math
import matplotlib.pyplot as plt
from IPython.display import Markdown
import os
from os import listdir
from os.path import isfile, join
os.chdir(r'C:\Users\chapmanvl\Documents\VC 2019 projects')
print(os.getcwd())

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

# seq to txt file
def seqtotxt(seq, name):
    with open(f'{name}.txt', 'w') as output:
        output.write(seq[0:])
    return

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

# Test case: X, Y = seq_to_chaosfile('ATTTGCATGAGGGGAGAT', 'testing')

# Get CGR co-ordinates and graph
def seq_to_chaosfile(seq, name):
    # clean seq first
    sequence = seq.strip().upper()        
    # warn the user if an obnoxoious character found in the file
    for uniq_char in  set(sequence):
        if uniq_char not in 'ATCGN':  
            print( f"a bad character found: {uniq_char}")
    pure_sequence = ''.join([ s for s in sequence if s in 'ATGC'])  
    X,Y, Z = one_cold_chaos_coord(pure_sequence)

    plt.plot(X, Y, 'b.', markersize=1)
    #plt.axis('off')
    plt.xlim(0,1)
    plt.ylim(0,1)
    new_name = os.path.splitext(name)[0] # to remove current file extension
    plt.savefig(f'{new_name}.png', bboxinches = 0, set_facecolor = 'white')
    return X, Y, Z 


"""
# Get a list of the full imaginary and real parts - redundant now as seq_to_chaos() and one_cold_chaos_coord() both return Z
def getcomplex(X, Y):         # skipping the first values (initiated at (0.5,0.5))
    Z = []
    for x,y in zip(X[1:],Y[1:]):
         Z.append(complex(x,y))
    return Z    

# Test case: Z = getcomplex(X, Y)
"""

# retrieve the power signal from the Z coordinates
def powersignal(cgr_z):
    input_signal = np.array(cgr_z).real
    #print ('sum_input-signal: ', np.mean(input_signal))
    dftZ = np.fft.fft(input_signal)
    print(dftZ[0])
    PSZ = np.abs(dftZ)**2
    mean_signal = PSZ[0]
    signal = PSZ[1:] #skipping first signal 
    #freq = np.fft.fftfreq((psd.shape[0]) )
    return mean_signal, signal
    

# Test case: mean_signal, signal = powersignal(Z)

"""
# plotting
X = freq
Y = psd
plt.plot(freq, psd)
plt.show()
"""
"""
# Test cases:
#### testing for GAATTC - works, therefore issues not with plotting chaos game representation
X,Y = seq_to_chaosfile("GAATTC", 'GAATTC')

####Testing with mus musculus - 
# plot not the same as published but similar to Jeffrey Nucleic Acids Res. 18:2163–2170, 1990
# Reason: T and G reversed in CGR in  (why would they do that????)
Mus = getseqsample('mus_musculus_x_chromo.fasta')
X, Y, Z = seq_to_chaosfile_plebs(Mus,'Mus')  

#### Testing with Buba buba GU571285 - mouse chromosome X CGR only identical when
# T and G reversed
def one_cold_chaos_coord_plebs(oligomer):
    ruleX = {'A': 0, 'C' : 0, 'G': 1, 'T' : 1}
    ruleY = {'A': 0, 'C' : 1, 'G': 0, 'T' : 1}
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

def seq_to_chaosfile_plebs(seq, name):
    # clean seq first
    sequence = seq.strip().upper()        
    # warn the user if an obnoxoious character found in the file
    for uniq_char in  set(sequence):
        if uniq_char not in 'ATCGN':  
            print( f"a bad character found: {uniq_char}")
    pure_sequence = ''.join([ s for s in sequence if s in 'ATGC'])  
    X,Y, Z = one_cold_chaos_coord_plebs(pure_sequence)

    plt.plot(X, Y, 'b.', markersize=1)
    #plt.axis('off')
    plt.xlim(0,1)
    plt.ylim(0,1)
    new_name = os.path.splitext(name)[0] # to remove current file extension
    plt.savefig(f'{new_name}.png', bboxinches = 0, set_facecolor = 'white')
    return X, Y, Z 

Buba = getseqsample('Buba_buba_GU571285.fasta')
X, Y, Z = one_cold_chaos_coord_plebs(Buba)

mean_signal, signal = powersignal(Z)
X_axis = range(len(signal))
plt.plot(X_axis, signal)
#plt.xlim(0,350)      # looks very similar to Hoang, Yin and Yau 2016 in Genomics
#plt.ylim(0,8000)
plt.show()


### Scaled to 1000
Buba = getseqsample('Buba_buba_GU571285.fasta')
X, Y, Z = one_cold_chaos_coord(Buba)
#real_Z = np.array(Z).real
#subtracted_Z = real_Z - np.mean(real_Z)
print(np.sum(subtracted_Z))
mean_signal, signal = powersignal(subtracted_Z)
print(mean_signal)
scaled_to_1000 = even_scale(signal, 1000)
X_axis = range(len(scaled_to_1000))
plt.plot(X_axis, scaled_to_1000)
#plt.xlim(0,350)      # looks very similar to Hoang, Yin and Yau 2016 in Genomics
#plt.ylim(0,8000)
plt.show()


scaled_psd, scaled_freq = scaledpowersignal(Z, 1000)
plt.plot(scaled_freq, scaled_psd)
plt.xlim(0,0.4)
#plt.ylim(0,8000)
plt.show()
"""

