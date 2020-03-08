
# A function to obtain the final chaos game coordinate of a genetic sequence
def ChaosGameFinal(raw_sequence):
 
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

# Initiate X and Y:
    X = 0
    Y = 0

# Perform sum of Ln/2^n to retrieve final chaos encoding algorithm, where 
# Ln is associated letter coordinate and n is position within the sequence   
    for i in range(len(sequence)):
        X += Lx[i]/2**(len(sequence) - i)
        Y += Ly[i]/2**(len(sequence) - i)
    return X, Y


# Time process for 1000 bases:
one_kbp_sequence = """ccttcgcacagccgcagtaccctctcggtttaggtgccgttcaaagcaatgccacagggg
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
%timeit time = ChaosGameFinal(one_kbp_sequence)
# 2.87 ms ± 143 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)
