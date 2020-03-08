# Finding CGR coordinates

Files enclosed are for the conversion of genetic sequences into Chaos Game Representation (CGR)
coordinates.

Processing_Fastas.py contains four functions useful for processing FASTA files:
    1) Extract sequence as a string object
    2) Extract sequence of given length and save as txt file
    3) Extract sequence, cut into given lengths and save each as a separate txt file
        in a new folder
    4) Extract sequence, cut into a given number of pieces of random 
        length and save each as a separate txt file in a new folder 

For encoding sequences as CGR coordinates, two method choices are available:
    1) CGR_coordinates.py uses the CGR algorithm to return all CGR coordinates for a     sequence 
    2) CGR_final_coordinates.py uses summation to find only the final 
    CGR coordinate for the sequence (from which all the other coordinates can be derived)


Additional:
    - Speed test folder:   Some examples of testing encoding algorithms with random sequences
    (obtained using SMS: Random DNA Sequence with default settings) to determine fastest 
    encoding method.
    Speed test results: The algorithmic method used in CGR_coordinates.py was found to encode 1kbp sequences in less than 1ms (with Intel i5-3470 processor and 8GB RAM),
    ~3 times faster than the summation method used in CGR_final_coordinates.py . 
