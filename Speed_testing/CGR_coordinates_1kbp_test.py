
#Function for getting cartesian coordinates for chaos game representations of genome sequences
def ChaosGame(raw_sequence):
 
# Ensure all bases uppercase
    upper_sequence = raw_sequence.upper()

#Filter for invlaid characters
    sequence = "".join([base for base in upper_sequence if base in "ATGC"])

# Initiate lists of coordinates for relevant corners for each letter
# in the sequence 
    Lx = []        
    Ly = []

# Create dictionaries of letter coordinates
    dictionaryX = {"G": 0, "T": 0, "A": 1, "C": 1}  
    dictionaryY = {"G": 0, "T": 1, "A": 1, "C": 0}  

# Append appropriate coordinates for each letter in the sequence
    for letter in sequence:
        Lx.append(dictionaryX[letter])
        Ly.append(dictionaryY[letter])

# Initiate final X and Y coordinate lists with starting point (0.5, 0.5):
    X = [0.5]
    Y = [0.5]

# Perform the chaos encoding algorithm:   
# Xn+1 = (Xn + Ln)/2   where Xn+1 is the next X value, Xn is current X value 
# and Ln is the coordinate of the corner for the relevant letter (A, T, C or G)
    for i in range(len(sequence)-1):
        X_value = (X[-1]+ Lx[i+1])/2
        Y_value = (Y[-1]+ Ly[i+1])/2
        X.append(X_value)
        Y.append(Y_value)
    
    return X, Y

# time process for 1000 bases:
random_sequence = """ccttcgcacagccgcagtaccctctcggtttaggtgccgttcaaagcaatgccacagggg
caaaccgttttatacgtgctaaccgcacgcgcggttgtgctttactgttggagcatcaag
attagagcccaggcacgggggttagtgggggccatcgttagtgtgtagacctgaattggc
acatatgatacgtgtcccagagtagccgaacaagttgcgacgagcaggtacggcgaacaa
gaggtagccaccccggaagacttataagggggtgttcatcttgcgtttctcttagtcgaa
aatatgacgccttccgtagggaatcgcacgagcgttgtcccctccgtagaaggcaatctt
actgactcagacacccatcctaccccttgatgcttgagggacaggtcacgccacccgtag
tccacttgtggtcaagactccttacgcccgattacgtgcataaatacattaattagtctg
tccgtggcttgcggtcaatgctgtagctaccgattatgtactccggttatgcgccaatag
agggtctcatcccaagtgcccccctgaccaggatagcgtgcaaccaagacattgtatctc
tgcatcgacaagatacatagtgcataattgcctaaaaagccgaccccagattaggaacac
taggatttcgtgctggctaggtgaggtaaaagtagcagaactcacaaaccacaggatact
ttctgatgtggagtgcatcaacaatacatactaatttcgacatattacatgatgttattt
cccgtggctggctacacataacaacatggattcatccctctttgtacctagtgttatgcg
tccctcgtactcacaggattttagacgcaagctatgactcgtcacggcatgttgcaacta
tagtccctgcttatgtacatcatggagatgggaagggaaacgttatcaggcttagcaaac
gctgggagaaaacaattacgtgacccgttagctggtagaa"""

# time process
%timeit time = ChaosGame(random_sequence)

#841 µs ± 13.6 µs per loop (mean ± std. dev. of 7 runs, 1000 loops each)