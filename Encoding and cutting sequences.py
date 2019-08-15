import numpy as np
import pandas as pd
import random as rnd
import math
import matplotlib.pyplot as plt
from IPython.display import Markdown


#Random sequence with non-nucleotides
def GenRanSeq(length, allowed_characters = ["A", "C", "G", "T", "N", "W" ]):
    return  "". join([rnd.choice(allowed_characters)  for i in range(length)])

#Random clean sequence
def GenRanSeqClean(length, allowed_characters = ["A", "C", "G", "T"]):
    return  "". join([rnd.choice(allowed_characters)  for i in range(length)])

# Write sequence to file
def write_to_file(input_seq,filer = 'randomseq.txt'):
    """ Writing to a txt file"""
    with open(filer,'w') as out:
         out.write(input_seq)
            
#One-hot encoding of sequence
def one_hot_nucl_liner(oligomer):
    # Make np array of zeros ready for values
    array_dim = (4, len(oligomer))
    One_hot_seq = np.zeros(array_dim)
    rule = {'A': 0, 'C' : 1, 'G': 2, 'T' : 3 }
    for i in range(len(oligomer)):
         if oligomer[i] in 'ACGT':
            One_hot_seq[rule[oligomer[i]] , i] = 1     
    return One_hot_seq

#Ordinal (one-cold) rather than one-hot encoding
def one_cold_nucl_liner(oligomer):
    # Make np array of zeros ready for values
    rule = {'A': 0.25, 'C' : 0.5, 'G': 0.75, 'T' : 1, 'N': 0 }
    one_cold_seq = [ rule[i] for i in oligomer ]
    return one_cold_seq

#Complex encoding on sequence
def one_complex_nucl_liner(oligomer):
    # Make np array of zeros ready for values
    rule = {'A': complex(1,1), 'C' : complex(1,-1), 'G': complex(-1,-1), 'T' : complex(-1,1), 'N': 0 }
    one_cold_seq = [ rule[i] for i in oligomer ]
    return one_cold_seq

# Complex encoding and chaos game representation co-ordinates
def one_cold_chaos_coord(oligomer):
    ruleX = {'A': 0, 'C' : 0, 'G': 1, 'T' : 1}
    ruleY = {'A': 0, 'C' : 1, 'G': 1, 'T' : 0}
    one_cold_worm_X = [ruleX[i] for i in oligomer ] 
    one_cold_worm_X.insert(0,0.5)
    one_cold_worm_Y = [ruleY[i] for i in oligomer ] 
    one_cold_worm_Y.insert(0,0.5)
    
    for i in range(len(one_cold_worm_X[1:])):
        h = i + 1
        one_cold_worm_X[h] = (one_cold_worm_X[i] + one_cold_worm_X[h])/2
    for i in range(len(one_cold_worm_Y[1:])):
        h = i + 1
        one_cold_worm_Y[h] = (one_cold_worm_Y[i] + one_cold_worm_Y[h])/2
        
    one_cold_worm = [ one_cold_worm_X, one_cold_worm_Y]
    return one_cold_worm

# Testing chaos game representation visualisation 
journal_sequence = one_cold_chaos_coord("GAATTC")
X = journal_sequence[0]
Y = journal_sequence[1]

lines, points = plt.subplots()
points.plot(X, Y, 'ro-')
points.set(xlim = (0, 1), ylim = (0,1))
for i, txt in enumerate(X):
    plt.annotate(txt, (X[i], Y[i]))
    
# Sequence to chaos game representation
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
    

# Cutting up sequences
# 1) Load file as string and clean 
def getseqsample(file):
    c
    with open (file) as seq:
        for line in seq:
            if line.startswith('>'):
                pass
            else:
                sequence = line.upper()
                pure_sequence_line = ''.join([ s for s in sequence if s in 'ATGC']
        sequence_sample = sequence[:1000]
    

"""# for multiple sequences - does not remove first line of fasta
def getseqsample_multi(file_string_list):
    pure_sequences = []
    seq_lengths = []
    for file in file_string_list:
        with open (file) as seq:
            sequence = seq.read().strip().upper()
            pure_sequence = ''.join([ s for s in sequence if s in 'ATGC']) 
            pure_sequences.append(pure_sequence)
            seq_lengths.append(len(pure_sequence))
    return pure_sequences, seq_lengths
"""

#2) find number of windows that will fit in each file
def cutandsample(sequence, length, windowsize, skipsize, minlim):
    #limit = number of windows that can be made
    limit = (length/ (windowsize + skipsize) )
    if limit < minlim:
        print( f"limit smaller than specified: {limit}")
        minlim = limit
    #slicing up sequences then 
    slicedsequence = []
    for i in range(minlim):
        startsite = (i * (windowsize + skipsize)) - (windowsize + skipsize)
        endsite = startsite + windowsize
        slicedsequence.append(sequence[startsite:endsite])
        joined = ''.join(slicedsequence) # usage of join
    return joined


# Version of cut and sample that collapses all sequences into one at the end
def cutandsample_collapsed(sequences, lengths, windowsize, skipsize, minlim):
    #limit = number of windows that can be made
    limit = min([x / (windowsize + skipsize) for x in lengths])
    if limit < minlim:
        print( f"limit smaller than specified: {limit}")
        minlim = limit
    #slicing up sequences then 
    slicedsequences = []
    for sequence in sequences:
        slicedsequence = []
        for i in range(minlim):
            startsite = (i * (windowsize + skipsize)) - (windowsize + skipsize)
            endsite = startsite + windowsize
            slicedsequence.append(sequence[startsite:endsite])
            joined = ''.join(slicedsequence) # usage of join
        slicedsequences.append(joined)
    return slicedsequences

# Code for cutting a portion of a sequence; if FASTA format, will, erase first line 

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
            else: # so that only takes from first Fasta
                No_of_fastas += 1  
                if No_of_fastas > 1:
                     print("first sequence shorter than maximum length")
                     break
        # Displays a red error message if pure_chunk shorter than required length
        display (Markdown(f'<span style="color: #ff0000">Error: Max length of unambiguous code too short:</span>{len(pure_chunk)}.'))
    

