# Code for cutting a portion of a sequence; if FASTA format, will, erase first line 
import numpy as np
import pandas as pd
import random as rand
import math
import matplotlib.pyplot as plt
from IPython.display import Markdown
import os

def getseqsample(file, output_name, max_length):
    pure_chunk = ''
    No_of_fastas = 0
    with open (file) as seq, open(f'{output_name}.txt', 'w') as output:
        for line in seq:
            if not line.startswith('>'):
                upperLine = line.upper().strip()
                pure_upper_line = ''.join([ s for s in upperLine if s in 'ATGC'])
                pure_chunk += pure_upper_line
                                         
                if len(pure_chunk) > max_length:
                    output.write(pure_chunk[:max_length] + '\n')                                       
                    return pure_chunk[:max_length]
            else:
                No_of_fastas += 1  
                if No_of_fastas > 1:
                     print("first sequence shorter than maximum length")
                        return
        # Displays a red error message if pure_chunk shorter than required length
        display (Markdown(f'<span style="color: #ff0000">Error: Max length of unambiguous code too short:</span>{len(pure_chunk)}.'))


# test sample: getseqsample("C:\\Users\\chapmanvl\\Documents\\VC 2019 projects\\HRV_input\\HRV_A.fasta", "test1", 5000)


# cut sequence into user-defined number of lengths of fixed length - needs clean sequence sto be input (have to run sequences through getseqsample first)
def cutandsample_fixed(file, output_folder_name, desired_length, number):
    # making a new directory, opening the file,determining sequence length and checking that sequence length is longer than desired length
    os.mkdir(f"{output_folder_name}")
    with open(file) as sequence:
        seq = ''.join([bases for bases in sequence if bases in 'ATGC']) 
        seq_length = len(seq)
        end_of_seq = seq_length - desired_length
        if seq_length < desired_length:
            print("Error: Sequence length is shorter than desired length - input a longer sequence or reduce desired length")
            return

    #slicing up sequences to desired length
        for n in range(number):
            with open(f'{output_folder_name}\\{n + 1}.txt', 'w') as output_file:
                startsite = (rand.randint(0, end_of_seq)) 
                endsite = startsite + desired_length
                output_file.write(seq[startsite:endsite] + '\n')

# test sample: cutandsample_fixed("C:\\Users\\chapmanvl\\Documents\\VC 2019 projects\\test1.txt", "testing_cut", 500, 10)


# cut sequence into user-defined number of lengths of random length - needs clean sequence sto be input (have to run sequences through getseqsample first)
def cutandsample_random(file, output_folder_name, number):
    # making a new directory, opening the file,determining sequence length and checking that sequence length is longer than desired length
    os.mkdir(f"{output_folder_name}")
    with open(file) as sequence:
        seq = ''.join([bases for bases in sequence if bases in 'ATGC']) 
        seq_length = len(seq)      
        lengths = []
    #slicing up sequences to desired length
        for n in range(number):
            with open(f'{output_folder_name}\\{n + 1}.txt', 'w') as output_file:
                startsite = (rand.randint(0, seq_length))
                possible_length = seq_length - startsite
                endsite = startsite + rand.randint(0,possible_length)
                output_file.write(seq[startsite:endsite] + '\n')
                lengths.append(len(seq[startsite:endsite]))
    return lengths  
# test sample: cutandsample_random("C:\\Users\\chapmanvl\\Documents\\VC 2019 projects\\test1.txt", "testing_cut_random", 10)

# Sequence to chaos game representation co-ordinates
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
    for i in range(len(one_cold_worm_Y[1:])):
        h = i + 1
        Y_complex[h] = (Y_complex[i] + Y_complex[h])/2
        
    coordinates = [X_complex, Y_complex]
    return coordinates

""" Test case of chaos game co-ordinates 
journal_sequence = one_cold_chaos_coord("GAATTC")
X = journal_sequence[0]
Y = journal_sequence[1]

lines, points = plt.subplots()
points.plot(X, Y, 'ro-')
points.set(xlim = (0, 1), ylim = (0,1))
for i, txt in enumerate(X):
    plt.annotate(txt, (X[i], Y[i])) """ 

# Sequence to chaos game representation diagram
def seq_to_chaos(seq):
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

    plt.plot(X, Y, '.', markersize=1)
    plt.show()
    
# test case: seq_to_chaos('GAATTC')