""" A file with the 3 main applications (just started week 3):
		- script to get the reversed complimentary seq of DNA 
		- GC% of a DNA seq
		- script to check for an in-frame stop codon
"""

seq = "AATGCGGGCCATTAGCATCAGAAATCGAAGCCTG"

def rev_compliment(dna):
	"Return a reverse compliment sequence of a given DNA sequence."
	change = {"a": "t", "g": "c", "c": "g", "t": "a", "n": "n"}
	rev = dna[::-1]
	compliment = [change[b] if b in list(change.keys()) else "n" for b in rev.lower()]
	return ''.join(compliment)

def gc_percent(dna):
	"Return a GC % (percent) of a given DNA sequence."
	dna = dna.lower()
	gc = dna.count("g")+dna.count("c")
	percent = gc*100.0/(len(dna)-dna.count("n"))
	return percent

def stop_codon(dna, start=0):
	"Return a True if a Stop codon is found in a given DNA seq, False otherwise."
	dna = dna.lower()
	stops = ["tga", "tag", "taa"]
	for i in range(start, len(dna), 3):
		if dna[i:i+3] in stops:
			return True
	return False


if __name__ == "__main__": 
	print(rev_compliment(seq))
	print(gc_percent(seq))
	print(stop_codon(seq, start=0))
	# Printing help
	print(help(stop_codon))