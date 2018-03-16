#Sorting dict keys

dna="aatgcaacgcattaa" # input("Enter DNA sequence:")

dna_counts={'t':dna.count('t'),'c':dna.count('c'),'g':dna.count('g'),'a':dna.count('a')}

nt=sorted(dna_counts.keys())

print(nt[-1])

# EXPLANATION
"""
Letters are sorted by its hexadecimal/decimal code.
BRAINFUCK - Try every time
"""

# Delete from dictionary
dna_counts={'g': 13, 'c': 3, 't': 1, 'a': 16} 
del dna_counts['a']

print(dna_counts)


# Count most frequent base
dna_counts = {'t':dna.count('t'), 'c':dna.count('c'), 'g':dna.count('g'), 'a':dna.count('a')}
max_freq=sorted(dna_counts.values())[-1]
print(max_freq)


# Adding & Updating Elements in a dictionary
someData = { }
someData['cheese'] = 'dairy'
someData['Cheese'] = 'dairy'
someData['Cheese'] = 'Dairy'
someData['cheese'] = 'Dairy'
print(someData)