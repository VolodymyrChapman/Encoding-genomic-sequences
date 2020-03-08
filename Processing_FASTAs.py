# Code for cutting a portion of a sequence; if FASTA format, will, erase first line 
import random as rand
from IPython.display import Markdown
import os


# read fasta sequence from file

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
            



# extract oligomer of certain length from fasta file and produce output file

def getseqsamples(file, output_name, max_length):
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
        # Displays a red error message if pure_chunk shorter than required length
        display (Markdown(f'<span style="color: #ff0000">Error: Max length of unambiguous code too short:</span>{len(pure_chunk)}.'))




# Cut 'cleaned' sequence data into user-defined lengths  
# Requires 'clean' sequences to be input 
# i.e. sequences have to be processed through getseqsample() first

def cutandsample_fixed(file, output_folder_name, desired_length, number):
    # making a new directory, opening the file,determining sequence length and checking that sequence length is longer than desired length
    os.mkdir(f"{output_folder_name}")
    with open(file) as sequence:
        seq = ''.join([bases for bases in sequence]) 
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



# Cut 'cleaned' sequence data into random lengths (to simulate shotgun sequence data) 
# Requires 'clean' sequences to be input 
# i.e. sequences have to be processed through getseqsample() first

def cutandsample_random(file, output_folder_name, number):
    # making a new directory, opening the file,determining sequence length and checking that sequence length is longer than desired length
    os.mkdir(f"{output_folder_name}")
    with open(file) as sequence:
        seq = ''.join([bases for bases in sequence]) 
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
