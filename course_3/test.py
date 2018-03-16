def valid_dna1(dna):
    for c in dna:
        if c in 'acgtACGT':
            return True
        else:
            return False


seq= "na"
print(valid_dna1(seq))