""" A script to get the complementary chain of DNA """

chain = "AATGCGGGCCATTAGCATCAGAAATCGAAGCCT"

def reverse(chain, rna=False):
	if not rna: change = {"A": "T", "G": "C", "C": "G", "T": "A"}
	else: 		change = {"A": "U", "G": "C", "C": "G", "T": "A"}

	rev = ""
	for b in chain: rev += change[b]
	return rev

if __name__ == "__main__": 
	print("Enter your chain")
	chain = str(input())
	print("Rna? [yes/no]")
	rna = True if "s" in str(input()) else False
	print(reverse(chain, rna=rna))