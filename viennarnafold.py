#!/usr/bin/python

import subprocess
import sys

def getstructure(seq):
	process = subprocess.Popen(["./RNAfold", "-noPS"], shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	process.stdin.write(seq)
	return process.communicate()[0].split()[1]

fn = sys.argv[1] # "5s.fasta"
f = open(fn)
conts = f.readlines()
inbody = False
#allseqs = []
curseq = ""
header = ""
i = 0
for x in conts:
	print >> sys.stderr, i, "/", len(conts)
	i += 1
	if x[0] == ">":
		if curseq != "":
			#allseqs.append((header, curseq, getstructure(curseq)))
			print header
			print curseq
			print getstructure(curseq)
			curseq = ""
		header = x[:len(x)-1] # trim off trailing \n
	else:
		curseq += x[:len(x)-1] # trim off trailing \n

if curseq != "":
	#allseqs.append((header, curseq, getstructure(curseq)))
	print header
	print curseq
	print getstructure(curseq)

#for header,curseq,structure in allseqs:
#	print header
#	print curseq
#	print structure

