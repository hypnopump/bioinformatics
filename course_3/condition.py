dnas =["atacagacantnnngngnacnggctctcaa", "atacagacatggacggctctcaa"]

# The MOST PYTHONIC way - enumerate gives pairs of index+object
for i, dna in enumerate(dnas): 	
	if "n" in dna and len(dna)>0:
		print("dna {0} has {1} undefined bases".format(i, dna.count("n")))
	else:
		print("dna {0} has not any undefined bases".format(i))

print("")

# Alternative with classical PYTHONIC methodology
for dna in dnas: 	
	if "n" in dna and len(dna)>0:
		print("dna {0} has {1} undefined bases and is {2} bases long".format(dnas.index(dna), dna.count("n"), len(dna)))
	else:
		print("dna {0} has not any undefined bases and is {1} bases long".format(dnas.index(dna), len(dna)))

print("")

# the C-Style - don't do it this way unless there's a good reason!
for i in range(len(dnas)): 	
	if "n" in dnas[i] and len(dnas[i])>0:
		print("dna {0} has {1} undefined bases".format(i, dnas[i].count("n")))
	else:
		print("dna {0} has not any undefined bases".format(i))