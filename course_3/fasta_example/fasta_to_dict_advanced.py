import sys
import getopt

def usage():
	print("""
This Script reads a FASTA file, extracts its entries
and save them to a dict. Two arguments are accepted:

python fasta_to_dic_advanced.py [-h] [-l <length>] <filename>

[-h]: prints this message
[-l <length>]: minimum length of the DNA sequences
			   in order to be extracted
<filename>: the file has top be in FastA format
	""")

o, a = getopt.getopt(sys.argv[1:], "l:h")

opts = {}
seqlen = 0
for k,v in o:
	opts[k] = v

if "-h" in opts.keys():
	usage(); sys.exit()
if len(a) < 1:
	usage(); sys.exit("Input FastA file is missing")
if "-l" in opts.keys():
	if int(opts["-l"])<0:
		print("Length of the sequence must be greater than 0."); sys.exit(0)
	seqlen = int(opts["-l"])

name = a[0]
f = open(name, "r")

seqs = {}
for line in f:
	line = line.rstrip()
	if line[0] == ">":
		name = line[1:].split()[0]
		seqs[name] = ""
	else:
		seqs[name] += line

# Delete short sequences:
short = []
for key in seqs.keys():
	if len(seqs[key]) < seqlen:
		short.append(key)

for key in short:
	del seqs[key]

print(seqs)
f.close()