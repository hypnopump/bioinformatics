name = "example.fa"

f = open(name, "r")

seqs = {}
for line in f:
	line = line.rstrip()
	if line[0] == ">":
		name = line[1:].split()[0]
		seqs[name] = ""
	else:
		seqs[name] += line

print(seqs)