# seq
Tools for processing sequence data using one-hot, one-cold and chaos encoding.
Files of interest may be as follows:

- Processing FASTAs.py contains a number of functions for extracting sequence data including whole-sequence extraction, of a desired length and a random length (for model building and sequence assembly training purposes) as well as some functions for chaos encoding

- even_scaled_signals_labelled.txt  pickled (using the Python Pickle library) dataset of influenza neurmanidase sequences ( plus an outgroup) that have undergone Fourier transform analysis, conversion into power signals and been even-scaled to equal sequence length (as outlined in: Numerical encoding of DNA sequences by chaos game representation with application in similarity comparison; Genomics; 2016; Hoang, Yin, Yau.

- CNN on influenza data.py contains a simple CNN (2 conv layers, 2 max pooling layers, 1 additional hidden layer, 1 dropout layer and an output layer) to classify influenza neuraminidase sequences using the even-scaled data above
 
