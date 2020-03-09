
# Function to obtain X and Y distances between coordinates on CGRs 
# in order to more clearly identify repeating patterns

def ChaosDistance(raw_sequence):
 
# Ensure all bases uppercase
    upper_sequence = raw_sequence.upper()

#Filter for invalid characters
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
        X_value = (X[-1]+ Lx[i])/2
        Y_value = (Y[-1]+ Ly[i])/2
        X.append(X_value)
        Y.append(Y_value)
    
# Calculate X and Y distance between adjacent points
    distance_X = []    # initiate distance lists
    distance_Y = []
    
    # calculate the X and Y distance between adjacent points and
    # append to distance_X/distance_Y lists, skipping the first point (0.5, 0.5) 
    for i in range(1, len(X)-1):
        dist_X = X[i+1] - X[i]
        dist_Y = Y[i+1] - Y[i]
        distance_X.append(dist_X)
        distance_Y.append(dist_Y)

    
    return distance_X, distance_Y

# Test case:
code = "ATTGACGT"
X_dist, Y_dist = ChaosDistance(code)