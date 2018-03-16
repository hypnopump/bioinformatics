dna=input("Enter a DNA sequence, please:")
# dna = "aatgaacctttaatg"

# o1 = dna.find('atg')
# hey = dna.find('atg',o1+1)
hey = dna.find('atg',dna.find('atg')+1) 

# Use dna.rfind to find the last occurrence of a seq in a string

print(hey)
print(len(dna))