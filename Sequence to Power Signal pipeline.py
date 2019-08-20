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
    

""" # Influenza A Neuraminidase sequence processing
Number_of_fastas = fastasplit(r'C:\Users\chapmanvl\Documents\VC 2019 projects\Influenza A segment 6 NA\extended_list_57_tagged.fasta', 'extended_list_57_tagged')
details = multiple_seq_to_ps_file(r'C:\Users\chapmanvl\Documents\VC 2019 projects\Influenza A segment 6 NA\extended_list_57_tagged', 'power_signals')
scaling_factor = details['length'].max()
folder = r'C:\Users\chapmanvl\Documents\VC 2019 projects\Influenza A segment 6 NA\extended_list_57_tagged\power_signals'
even_signals = even_scale_multi(folder, scaling_factor)

#######plotting
col = even_signals.columns
even_signals.insert(0, 'frequencies', range(1467))
even_signals.plot(kind='line',x= 'frequencies' ,y= col, legend = None)

black = ['blackCY0147881InfluenzaAvirusAturkeyMinnesota11988H7N9segment6completesequence', 
         'blackCY1860041InfluenzaAvirusAmallardMinnesotaAI0937702009H7N9neuraminidaseNAgenecompletecds',
         'blackCY2353641InfluenzaAvirusAreassortantIDCDCRG56BHongKong1252017XPuertoRico81934H7N9neuraminidaseNAgenecompletecds']
red = ['redAF5091022InfluenzaAvirusAChickenHongKong822101H5N1neuraminidaseNAgenecompletecds',
        'redAM9140171InfluenzaAvirusAdomesticduckGermanyR177207H5N1N1geneforneuraminidasegenomicRNA',
        'redEF5414641InfluenzaAvirusAchickenKoreaes2003H5N1segment6neuraminidaseNAgenecompletecds']
blue = ['blueAB4706631InfluenzaAvirusAduckHokkaidow732007H1N1NAgeneforneuraminidasecompletecds',
       'blueAB5461591InfluenzaAvirusApintailMiyagi14722008H1N1NAgeneforneuraminidasecompletecds',
       'blueAM1573581InfluenzaAvirusAmallardFrance6912002H1N1NAgeneforneuraminidasegenomicRNA']
green = ['greenAY6460801InfluenzaAvirusAchickenBritishColumbiaGSChumanB04H7N3neuraminidaseNAgenecompletecds',
       'greenCY0393211InfluenzaAvirusAavianDelawareBay2262006H7N3segment6completesequence',
       'greenCY0762311InfluenzaAvirusAAmericangreenwingedtealCalifornia442429062007H7N3segment6completesequence']
purple = ['purpleAY2099251InfluenzaAvirusATaiwan167H2N2neuraminidaseNAgenecompletecds',
       'purpleAY2099281InfluenzaAvirusAGeorgia167H2N2neuraminidaseNAgenecompletecds',
       'purpleAY2099321InfluenzaAvirusAKorea42668H2N2neuraminidaseNAgenecompletecds']       


# Making plot of 3 meaned signals of each class
test = even_signals[:]
test.insert(0, 'mean_black', test[black].mean(axis=1))
test.insert(0, 'mean_blue', test[blue].mean(axis=1))
test.insert(0, 'mean_red', test[red].mean(axis=1))
test.insert(0, 'mean_green', test[green].mean(axis=1))
test.insert(0, 'mean_purple', test[purple].mean(axis=1))
mcdf = test.iloc[:,0:5]
col = mcdf.columns
mcdf.insert(0, 'frequencies', range(1467))
mcdf.plot(kind='line',x='frequencies',y=col ,color=['magenta', 'green', 'red', 'blue', 'black'])



#### making df of labelled data
''' # Influenza A Neuraminidase sequence processing
Number_of_fastas = fastasplit(r'C:\Users\chapmanvl\Documents\VC 2019 projects\Influenza A segment 6 NA\extended_list_57_tagged.fasta', 'extended_list_57_tagged')
details = multiple_seq_to_ps_file(r'C:\Users\chapmanvl\Documents\VC 2019 projects\Influenza A segment 6 NA\extended_list_57_tagged', 'power_signals')
scaling_factor = details['length'].max()
folder = r'C:\Users\chapmanvl\Documents\VC 2019 projects\Influenza A segment 6 NA\extended_list_57_tagged\power_signals'
even_signals = even_scale_multi(folder, scaling_factor)

black = ['blackCY0147881InfluenzaAvirusAturkeyMinnesota11988H7N9segment6completesequence',
       'blackCY1860041InfluenzaAvirusAmallardMinnesotaAI0937702009H7N9neuraminidaseNAgenecompletecds',
       'blackCY2353641InfluenzaAvirusAreassortantIDCDCRG56BHongKong1252017XPuertoRico81934H7N9neuraminidaseNAgenecompletecds',
       'blackGU0604841InfluenzaAVirusAgooseCzechRepublic1848K92009H7N9segment6neuraminidaseNAgenecompletecds',
       'blackKC6098011InfluenzaAvirusAwildduckKoreaSH19472010H7N9segment6neuraminidaseNAgenecompletecds',
       'blackKC8532311InfluenzaAvirusAShanghai4664T2013H7N9segment6neuraminidaseNAgenecompletecds',
       'blackKF2596881InfluenzaAvirusAduckJiangxi30962009H7N9segment6neuraminidaseNAgenecompletecds',
       'blackKF2597341InfluenzaAvirusAchickenRizhao7132013H7N9segment6neuraminidaseNAgenecompletecds',
       'blackKF6095291InfluenzaAvirusAtreesparrowShanghai012013H7N9segment6neuraminidaseNAgenecompletecds',
       'blackKF9389451InfluenzaAvirusAchickenJiangsu10212013H7N9segment6neuraminidaseNAgenecompletecds',
       'blackKY7511241InfluenzaAvirusAchickenGuangdongGD152016H7N9segment6neuraminidaseNAgenecompletecds']
red = ['redAB6841611InfluenzaAvirusAchickenMiyazaki102011H5N1NAgeneforneuraminidasecompletecds',
       'redAF5091022InfluenzaAvirusAChickenHongKong822101H5N1neuraminidaseNAgenecompletecds',
       'redAM9140171InfluenzaAvirusAdomesticduckGermanyR177207H5N1N1geneforneuraminidasegenomicRNA',
       'redEF5414641InfluenzaAvirusAchickenKoreaes2003H5N1segment6neuraminidaseNAgenecompletecds',
       'redEU6358751InfluenzaAvirusAchickenYunnanchuxiong012005H5N1neuraminidaseNAgenecompletecds',
       'redFM1771211InfluenzaAvirusAchickenGermanyR32342007H5N1NAgeneforneuraminidase',
       'redGU1865111InfluenzaAvirusAturkeyVA505477182007H5N1segment6neuraminidaseNAgenecompletecds',
       'redHQ1853811InfluenzaAvirusAchickenEasternChinaXH2222008H5N1neuraminidaseNAgenecompletecds',
       'redHQ1853831InfluenzaAvirusAduckEasternChinaJS0172009H5N1segment6neuraminidaseNAgenecompletecds',
       'redJF6996771InfluenzaAvirusAmandarinduckKoreaK104832010H5N1segment6neuraminidaseNAgenecompletecds',
       'redKF5724351InfluenzaAvirusAwildbirdHongKong0703512011H5N1segment6neuraminidaseNAgenecompletecds']
blue = ['blueAB4706631InfluenzaAvirusAduckHokkaidow732007H1N1NAgeneforneuraminidasecompletecds',
       'blueAB5461591InfluenzaAvirusApintailMiyagi14722008H1N1NAgeneforneuraminidasecompletecds',
       'blueAM1573581InfluenzaAvirusAmallardFrance6912002H1N1NAgeneforneuraminidasegenomicRNA',
       'blueCY1385621InfluenzaAvirusAmallardNovaScotia000882010H1N1neuraminidaseNAgenecompletecds',
       'blueCY1400471InfluenzaAvirusAmallardMinnesotaSg006202008H1N1neuraminidaseNAgenecompletecds',
       'blueCY1496301InfluenzaAvirusAthickbilledmurreCanada18712011H1N1neuraminidaseNAgenecompletecds',
       'blueEU0260462InfluenzaAvirusAmallardMaryland3522002H1N1segment6neuraminidaseNAgenecompletecds',
       'blueFJ3571141InfluenzaAvirusAmallardMD262003H1N1segment6neuraminidaseNAgenecompletecds',
       'blueGQ4118941InfluenzaAvirusAdunlinAlaska444216602008H1N1segment6neuraminidaseNAgenecompletecds',
       'blueHM3709691InfluenzaAvirusAturkeyOntarioFAV11042009H1N1segment6neuraminidaseNAgenecompletecds',
       'blueHQ8979661InfluenzaAvirusAmallardKoreaKNUYP092009H1N1segment6neuraminidaseNAgenecompletecds',
       'blueKC6081601InfluenzaAvirusAduckGuangxi030D2009H1N1segment6neuraminidaseNAgenecompletecds',
       'blueKM2440781InfluenzaAvirusAturkeyVirginia41352014H1N1segment6neuraminidaseNAgenecompletecds']
green = ['greenAY6460801InfluenzaAvirusAchickenBritishColumbiaGSChumanB04H7N3neuraminidaseNAgenecompletecds',
       'greenCY0393211InfluenzaAvirusAavianDelawareBay2262006H7N3segment6completesequence',
       'greenCY0762311InfluenzaAvirusAAmericangreenwingedtealCalifornia442429062007H7N3segment6completesequence',
       'greenCY1293361InfluenzaAvirusAAmericanblackduckNewBrunswick024902007H7N3neuraminidaseNAgenecompletecds',
       'greenEU0309701InfluenzaAvirusAshorebirdDelaware2206H7N3segment6neuraminidaseNANAgenecompletecds',
       'greenEU0309861InfluenzaAvirusAlaughinggullDelaware4206H7N3segment6neuraminidaseNANAgenecompletecds',
       'greenEU5008541InfluenzaAvirusAAmericanblackduckNB25382007H7N3segment6completesequence',
       'greenFJ6383281InfluenzaAvirusAchickenKarachiSPVC52004H7N3segment6neuraminidaseNAgenecompletecds',
       'greenHM1507381InfluenzaAvirusAchickenMurreeNARC11995H7N3segment6neuraminidaseNAgenecompletecds',
       'greenKF0420671InfluenzaAvirusAduckZhejiangDK102013H7N3segment6neuraminidaseNAgenecompletecds',
       'greenKX2896531InfluenzaAvirusAchickenPueblaCPA0330916CENASA950762016H7N3segment6neuraminidaseNAgenecompletecds',
       'greenMN2080521InfluenzaAvirusAnorthernshovelerEgyptMBD695C2016H7N3segment6neuraminidaseNAgenecompletecds']
purple = ['purpleAY2099251InfluenzaAvirusATaiwan167H2N2neuraminidaseNAgenecompletecds',
       'purpleAY2099281InfluenzaAvirusAGeorgia167H2N2neuraminidaseNAgenecompletecds',
       'purpleAY2099321InfluenzaAvirusAKorea42668H2N2neuraminidaseNAgenecompletecds',
       'purpleAY2099331InfluenzaAvirusABerkeley168H2N2neuraminidaseNAgenecompletecds',
       'purpleCY0055401InfluenzaAvirusAduckHongKong3191978H2N2segment6completesequence',
       'purpleCY1219611InfluenzaAvirusAwhitefrontedgooseNetherlands221999H2N2neuraminidaseNAgenecompletecds',
       'purpleCY1258881InfluenzaAvirusANed651963H2N2neuraminidaseNAgenecompletecds',
       'purpleDQ0174871InfluenzaAvirusAmallardPostdam178483H2N2fromGermanysegment6completesequence',
       'purpleJX0811421InfluenzaAvirusAemperorgooseAlaska442972602007H2N2segment6neuraminidaseNAgenecompletecds',
       'purpleL373292InfluenzaAvirusALeningrad13457H2N2neuraminidaseNAgenecompletecds']       

labelled = even_signals[:]
labelled = labelled.T
labelled['colour'] = 0
labelled.loc[black, 'colour'] = 0.00
labelled.loc[blue, 'colour'] = 0.25
labelled.loc[green, 'colour'] = 0.50
labelled.loc[purple, 'colour'] = 0.75
labelled.loc[red, 'colour'] = 1.00
labelled.tail()

# change index (row) names to be shorter
blacknum = []
bluenum = []
greennum = []
purplenum = []
rednum = []
for i in range(14):
    blacknum.append(join(f'black{i}'))
    bluenum.append(join(f'blue{i}')) 
    greennum.append(join(f'green{i}')) 
    purplenum.append(join(f'purple{i}')) 
    rednum.append(join(f'red{i}')) 

colnames = []
colnames.extend(blacknum[1:12])
colnames.extend(bluenum[1:14]) 
colnames.extend(greennum[1:13]) 
colnames.extend(purplenum[1:11]) 
colnames.extend(rednum[1:12])

labelled.index = colnames

# splitting into test/training sets
from sklearn.model_selection import train_test_split
X = labelled.iloc[:, 0:1467]
Y = labelled['colour']

X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.3)
print(X_train.shape, y_train.shape)
print(X_test.shape, y_test.shape)
