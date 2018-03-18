# Step1: Use blast to search in the database for
# The procedence of the DNA seq
from Bio.Blast import NCBIWWW
dna = "TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCACCTACGGTAGAG"
print("hello")
result_handle = NCBIWWW.qblast("blastn", "nt", dna)
print(result_handle)

# Step2: Use XML module to handle the response
from Bio.Blast import NCBIXML
blast_record = NCBIXML.read(result_handle)
print(blast_record) 

# Step3: Parse BLAST output
print(len(blast_record.alignments))

E_VALUE_TRESH = 0.01
for alignment in blast_record.alignments:
	for hsp in alignment.hsps:
		if hsp.expect < E_VALUE_TRESH:
			print("****Alignment****")
			print("sequence: ", alignment.title)
			print("length: ", alignment.length)
			print("e_value: ", hsp.expect)
			print(hsp.query)
			print(hsp.match)
			print(hsp.sbjct)