# initiation and functions
import math
import numpy as np
import pandas as pd
import os
from os import listdir
from os.path import isfile, join
from sklearn.pipeline import Pipeline
import pickle
os.chdir(r'C:\Users\chapmanvl\Documents\VC 2019 projects')
print(os.getcwd())

# Split a FASTA with many sequences into separate FASTAs
def fastasplit(file, output_folder):
    No_of_fastas = 0 
    start, finish = os.path.split(file)
    output_folder = join(start,output_folder)
    os.mkdir(output_folder)
    with open (file) as seq:
        for line in seq:
            if line.startswith('>'): 
                if No_of_fastas == 0: 
                    pass 
                else:
                    output.write(pure_chunk[:])
                    output.close() 
                pure_chunk = ''
                filename = ''.join(e for e in line if e.isalnum())
                output = open(f'{output_folder}\\{filename}.txt', 'w')
                No_of_fastas += 1 
                continue
            else:
                pure_chunk += ''.join([ s for s in line if s in 'ATGC'])    
        output.write(pure_chunk[:])
        output.close()                            
    return No_of_fastas

# test case:   a = fastasplit(r'C:\Users\chapmanvl\Documents\VC 2019 projects\Influenza A segment 6 NA\sequence.fasta', 'split')   


# Extract sequence name, sequence and length from FASTA file
def getseqsample(file):
    pure_chunk = ''
    with open (file) as seq:
        for line in seq:
            if line.startswith('>'):
                name = line
            else: 
                name = file
                upperLine = line.upper().strip()
                pure_upper_line = ''.join([ s for s in upperLine if s in 'ATGC'])
                pure_chunk += pure_upper_line
    return name, pure_chunk, len(pure_chunk) 


#get CGR Z co-ordinates only
def getz(oligomer):
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
    return Z

#Power signal
def powersignal(cgr_z):
    dftZ = np.fft.fft(cgr_z)
    PSZ = np.abs(dftZ)**2
    mean_signal = PSZ[0]
    signal = PSZ[1:] #skipping first signal 
    return mean_signal, signal

# Sequence to power signal   !!! Signal always n-1 of seq length as 1st signal removed!!!
def seq_to_ps(file):
    name, oligomer, length = getseqsample(file)
    z = getz(oligomer)
    mean_signal, signal = powersignal(z)
    return mean_signal, signal, length, name

#test case:
#mean_signal, signal, length = seq_to_ps(r'C:\Users\chapmanvl\Documents\VC 2019 projects\testing_cut\1.txt')

# Read all files in a folder and find power signal, save to pickle file in a new folder

def multiple_seq_to_ps_file(folder, output_folder_name):
    os.mkdir(join(folder,output_folder_name)) # make new folder
    files = [file for file in listdir(folder) if isfile(join(folder, file))] #only use files
    details = [] # initiate details list
    for file in files:
        mean_signal, signal, length_string, name = seq_to_ps(join(folder,file))
        path, name = os.path.split(file)
        length = int(length_string)
        details.append([name, length, mean_signal])
        with open(join(folder, output_folder_name,file), 'wb') as f:
            pickle.dump(signal, f)
    labels = ['tag', 'length', 'mean signal']
    details_df = pd.DataFrame.from_records(details, columns=labels)
    return details_df

# test case: details = multiple_seq_to_ps_file(r'C:\Users\chapmanvl\Documents\VC 2019 projects\testing_cut_random', 'power_signals')


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

#even scale all files in a folder
def even_scale_multi(folder, output_scale):
    #os.mkdir(join(folder,output_folder_name)) # make new folder
    files = [file for file in listdir(folder) if isfile(join(folder, file))] #only use files
    even_scaled_signals_list = []
    tags = []
    for file in files:
        name = os.path.splitext(file)[0]
        tags.append(name)
        unpickled = pickle.load(open(join(folder,file), 'rb' ))
        unpickled_list = unpickled.tolist()
        even_scaled_signal = even_scale(unpickled_list, output_scale)
        even_scaled_signals_list.append(even_scaled_signal)
    even_scaled_signals_transposed = pd.DataFrame(even_scaled_signals_list, index = tags)
    even_scaled_signals = even_scaled_signals_transposed.T
    with open(join(folder,'even_scaled_signals.txt'), 'wb') as f:
        pickle.dump(even_scaled_signals, f)
    print('even scaling complete')
    return even_scaled_signals


Number_of_fastas = fastasplit(r'C:\Users\chapmanvl\Documents\VC 2019 projects\outgroup\sequence.fasta', 'outgroup')
details = multiple_seq_to_ps_file(r'C:\Users\chapmanvl\Documents\VC 2019 projects\outgroup\outgroup', 'power_signals')
scaling_factor = 1467
folder = r'C:\Users\chapmanvl\Documents\VC 2019 projects\outgroup\outgroup\power_signals'
even_signals = even_scale_multi(folder, scaling_factor)
outgroup =['AB2209871Humanherpesvirus1UL42geneforDNApolymeraseprocessivityfactorcompletecdsstrainYMS1732',
       'AY4505301Humanrhinovirus17polyproteinmRNApartialcds',
       'KF5416361Enterovirus122isolateMNKRCMH5polyproteingenepartialcds',
       'KX1394401EchovirusE18isolateJenaVI8743103Dpolyproteingenepartialcds',
       'KY8204191HIV1isolateHP1495KTG110854fromSouthKoreaenvelopeglycoproteinenvgenepartialcds',
       'M177111Coxackieviruspolymerasegene3end',
       'M192681Vcholeraeneuraminidasegene5end',
       'M835621VibriocholeraeneuraminidasenanHgenecompletecdsfirst1348bp',
       'MF0014851PapayaringspotvirusisolateBarabanki1truncatedcoatproteingenepartialcds',
       'MF5110611UgandancassavabrownstreakvirusisolateIturi579polyproteingenepartialcds',
       'MG9834331HepacivirusCisolate152073523206NS3genepartialcds',
       'MN1788811HIV1isolateNIRT70fromIndiareversetranscriptasepolgenepartialcds',
       'S771321BGLF2ORFEpsteinBarrvirusEBVmRNA1449nt']

col = even_signals.columns
even_signals.insert(0, 'frequencies', range(1467))
even_signals.plot(kind='line',x= 'frequencies' ,y= col, legend = None)