"""

- Fetch first 5K sequences 

"""

seq =''
with open('Triticum_aestivum.IWGSC.dna.chromosome.1A.fa') as inp, open('trit.1A-1K.txt','w') as outf:
  for line in inp:
    if line.startswith('>'):
       print line
    else:
       seq += line.strip()
       if len(seq) > 5000:
          outf.write(seq[:5000] + '\n') 
          break

print len(seq)

       
             
