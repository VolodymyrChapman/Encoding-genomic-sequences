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

def one_cold_chaos_coord(oligomer):
    ruleX = {'A': 0, 'C' : 0, 'G': 1, 'T' : 1}
    ruleY = {'A': 0, 'C' : 1, 'G': 1, 'T' : 0}
    X_complex = [ruleX[i] for i in oligomer ] 
    X_complex.insert(0,0.5)
    Y_complex = [ruleY[i] for i in oligomer ] 
    Y_complex.insert(0,0.5)
    
    for i in range(len(X_complex[1:])):
        h = i + 1
        X_complex[h] = (X_complex[i] + X_complex[h])/2
    for i in range(len(Y_complex[1:])):
        h = i + 1
        Y_complex[h] = (Y_complex[i] + Y_complex[h])/2
        
    coordinates = [X_complex, Y_complex]
    return coordinates

# For exporting to a file
def seq_to_chaosfile(seq, name):
    # clean seq first
    sequence = seq.strip().upper()        
    # warn the user if an obnoxoious character found in the file
    for uniq_char in  set(sequence):
        if uniq_char not in 'ATCGN':  
            print( f"a bad character found: {uniq_char}")
    pure_sequence = ''.join([ s for s in sequence if s in 'ATGC'])  
    chaos = one_cold_chaos_coord(pure_sequence)
    X = chaos[0]
    Y = chaos[1]

    plt.plot(X, Y, 'k.', markersize=1)
    plt.axis('off')
    plt.xlim(0,1)
    plt.ylim(0,1)
    new_name = os.path.splitext(name)[0] # to remove current file extension
    plt.savefig(f'{new_name}.png', bboxinches = 0, set_facecolor = 'white')
    return X, Y 

# Test case:    seq_to_chaosfile('ATTTGCATGAGGGGAGAT', 'testing')

# Function for reading all sequence files within a folder and returning CGR plots for all to specified new folder
def readandCGR(folder, output_folder_name):
    os.mkdir(f"{output_folder_name}")
    files = [file for file in listdir(folder) if isfile(join(folder, file))]
    for file in files:
        with open(join(folder,file)) as sequence:
            seq = ''.join([bases for bases in sequence]) 
            seq_to_chaosfile(seq, f'{join(output_folder_name, file)}')

#Test case:  readandCGR(r'C:\Users\chapmanvl\Documents\VC 2019 projects\testing_cut_random', 'test_folder2')


# Writes into same folder as read
def readandCGRi(folder):
    files = [file for file in listdir(folder) if isfile(join(folder, file))]
    for file in files:
        with open(join(folder,file)) as sequence:
            seq = ''.join([bases for bases in sequence]) 
            seq_to_chaosfile(seq, f'{join(folder, file)}')

#Test case:     readandCGRi(r'C:\Users\chapmanvl\Documents\VC 2019 projects\testing_cut_random')