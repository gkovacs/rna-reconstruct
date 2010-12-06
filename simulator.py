#!/usr/bin/python

import mutationpath
from mutationpath import PhylogeneticTree
import transitionprobs
import math
import subprocess
import sys
import random

possiblepairs = [("A", "U"), ("U", "A"), ("C", "G"), ("G", "C"), ("G", "U"), ("U", "G")]

def groupTripletBasePairs(sequence, structure):
	assert len(sequence) == len(structure)
	leftParens = [] # (index, base)
	basepairs = [] # (leftidx, rightidx, leftbase, rightbase)
	for i in range(len(sequence)):
		if structure[i] == "(":
			leftParens.append((i, sequence[i]))
		elif structure[i] == ")":
			leftidx,leftbase = leftParens.pop()
			if (leftbase,sequence[i]) in possiblepairs:
				basepairs.append((leftidx,i,leftbase,sequence[i]))
	basepairs.sort()
	#print basepairs
	triplets = [] # (left3, paired3) ie AUGCAU -> (( (0,A), (1,U), (2,G) ), ( (5,U), (4,A), (3,C) ))
	for i in range(1,len(basepairs)-1):
		if basepairs[i][0] == basepairs[i-1][0]+1 == basepairs[i+1][0]-1:
			if basepairs[i][1] == basepairs[i-1][1]-1 == basepairs[i+1][1]+1:
				left3 = ( (basepairs[i-1][0], basepairs[i-1][2]), (basepairs[i][0], basepairs[i][2]), (basepairs[i+1][0], basepairs[i+1][2]) )
				right3 = ( (basepairs[i-1][1], basepairs[i-1][3]), (basepairs[i][1], basepairs[i][3]), (basepairs[i+1][1], basepairs[i+1][3]) )
				triplets.append((left3,right3))
	return triplets

#print groupTripletBasePairs("CAUGCAUGAU", "(((())))()")

def groupstacks(sequence, structure):
	assert len(sequence) == len(structure)
	leftParens = [] # (index, base)
	basepairs = [] # (leftidx, rightidx, leftbase, rightbase)
	for i in range(len(sequence)):
		if structure[i] == "(":
			leftParens.append((i, sequence[i]))
		elif structure[i] == ")":
			leftidx,leftbase = leftParens.pop()
			if (leftbase,sequence[i]) in possiblepairs:
				basepairs.append((leftidx,i,leftbase,sequence[i]))
	if len(basepairs) == 0:
		return []
	basepairs.sort()
	allStacks = []
	lIdx = -9
	rIdx = -9
	currentStack = []
	for x in basepairs: # (leftidx, rightidx, leftbase, rightbase)
		if x[0] == lIdx+1 and x[1] == rIdx-1:
			lIdx = x[0]
			rIdx = x[1]
			currentStack.append((x[2],x[3]))
		else:
			if len(currentStack) > 0:
				allStacks.append(currentStack)
				currentStack = []
			lIdx = x[0]
			rIdx = x[1]
			currentStack.append((x[2],x[3]))
	if len(currentStack) > 0:
		allStacks.append(currentStack)
	return allStacks

#print groupstacks("CAUGCAUGAU", "(((())))()")

'''
def groupstacks(sequence, structure):
	assert len(sequence) == len(structure)
	allStacks = []
	currentStackLeft = []
	currentStackRight = []
	for i in range(len(structure)):
		if structure[i] == '(':
			currentStackLeft.append(sequence[i])
		elif structure[i] == ')':
			assert len(currentStackLeft) > len(currentStackRight)
			currentStackRight.append(sequence[i])
			if len(currentStackLeft) == len(currentStackRight):
				currentStackRight.reverse()
				allStacks.append(zip(currentStackLeft, currentStackRight))
				currentStackLeft = []
				currentStackRight = []
	return allStacks
'''

def getstructure(seq):
	process = subprocess.Popen(["./RNAfold", "-noPS"], shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	process.stdin.write(seq)
	return process.communicate()[0].split()[1]

class Simulator:
	def __init__(self, startSeq, structure=None):
		if structure == None:
			structure = getstructure(startSeq)
		self.startSeq = startSeq
		self.structure = structure
		#self.startStacks = groupstacks(startSeq, structure)
		self.startTriplets = groupTripletBasePairs(startSeq, structure)
		self.ancestor = Mutator(startSeq, structure)
	
	def fourDescendentUltrametric(self, branchLength):
		intermediate1 = self.ancestor.getResultOfMutations(branchLength)
		intermediate2 = self.ancestor.getResultOfMutations(branchLength)
		child1a = intermediate1.getResultOfMutations(branchLength)
		child1b = intermediate1.getResultOfMutations(branchLength)
		child2a = intermediate2.getResultOfMutations(branchLength)
		child2b = intermediate2.getResultOfMutations(branchLength)
		return PhylogeneticTree(self.ancestor.seq, PhylogeneticTree(intermediate1.seq, child1a.seq, child1b.seq), PhylogeneticTree(intermediate2.seq, child2a.seq, child2b.seq))
		#return (intermediate1.seq,intermediate2.seq), (child1a.seq,child1b.seq,child2a.seq,child2b.seq)
	
	def eightDescendentUltrametric(self, branchLength):
		intermediate1 = self.ancestor.getResultOfMutations(branchLength)
		intermediate2 = self.ancestor.getResultOfMutations(branchLength)
		child1a = intermediate1.getResultOfMutations(branchLength)
		child1b = intermediate1.getResultOfMutations(branchLength)
		child2a = intermediate2.getResultOfMutations(branchLength)
		child2b = intermediate2.getResultOfMutations(branchLength)
		child1a1 = child1a.getResultOfMutations(branchLength)
		child1a2 = child1a.getResultOfMutations(branchLength)
		child1b1 = child1b.getResultOfMutations(branchLength)
		child1b2 = child1b.getResultOfMutations(branchLength)
		child2a1 = child2a.getResultOfMutations(branchLength)
		child2a2 = child2a.getResultOfMutations(branchLength)
		child2b1 = child2b.getResultOfMutations(branchLength)
		child2b2 = child2b.getResultOfMutations(branchLength)
		return PhylogeneticTree(self.ancestor.seq, PhylogeneticTree(intermediate1.seq, PhylogeneticTree(child1a.seq, child1a1.seq, child1a2.seq), PhylogeneticTree(child1b.seq, child1b1.seq, child1b2.seq)), PhylogeneticTree(intermediate2.seq, PhylogeneticTree(child2a.seq, child2a1.seq, child2a2.seq), PhylogeneticTree(child2b.seq, child2b1.seq, child2b2.seq)))
	
	def sixteenDescendentUltrametric(self, branchLength):
		intermediate1 = self.ancestor.getResultOfMutations(branchLength)
		intermediate2 = self.ancestor.getResultOfMutations(branchLength)
		child1a = intermediate1.getResultOfMutations(branchLength)
		child1b = intermediate1.getResultOfMutations(branchLength)
		child2a = intermediate2.getResultOfMutations(branchLength)
		child2b = intermediate2.getResultOfMutations(branchLength)
		child1a1 = child1a.getResultOfMutations(branchLength)
		child1a2 = child1a.getResultOfMutations(branchLength)
		child1b1 = child1b.getResultOfMutations(branchLength)
		child1b2 = child1b.getResultOfMutations(branchLength)
		child2a1 = child2a.getResultOfMutations(branchLength)
		child2a2 = child2a.getResultOfMutations(branchLength)
		child2b1 = child2b.getResultOfMutations(branchLength)
		child2b2 = child2b.getResultOfMutations(branchLength)
		
		c1 = child1a1.getResultOfMutations(branchLength)
		c2 = child1a1.getResultOfMutations(branchLength)
		c3 = child1a2.getResultOfMutations(branchLength)
		c4 = child1a2.getResultOfMutations(branchLength)
		c5 = child1b1.getResultOfMutations(branchLength)
		c6 = child1b1.getResultOfMutations(branchLength)
		c7 = child1b2.getResultOfMutations(branchLength)
		c8 = child1b2.getResultOfMutations(branchLength)
		
		c9 = child2a1.getResultOfMutations(branchLength)
		c10 = child2a1.getResultOfMutations(branchLength)
		c11 = child2a2.getResultOfMutations(branchLength)
		c12 = child2a2.getResultOfMutations(branchLength)
		c13 = child2b1.getResultOfMutations(branchLength)
		c14 = child2b1.getResultOfMutations(branchLength)
		c15 = child2b2.getResultOfMutations(branchLength)
		c16 = child2b2.getResultOfMutations(branchLength)

		return PhylogeneticTree(self.ancestor.seq, PhylogeneticTree(intermediate1.seq, PhylogeneticTree(child1a.seq, PhylogeneticTree(child1a1.seq, c1.seq, c2.seq), PhylogeneticTree(child1a2.seq, c3.seq, c4.seq)), PhylogeneticTree(child1b.seq, PhylogeneticTree(child1b1.seq, c5.seq, c6.seq), PhylogeneticTree(child1b2.seq, c7.seq, c8.seq))), PhylogeneticTree(intermediate2.seq, PhylogeneticTree(child2a.seq, PhylogeneticTree(child2a1.seq, c9.seq, c10.seq), PhylogeneticTree(child2a2.seq, c11.seq, c12.seq)), PhylogeneticTree(child2b.seq, PhylogeneticTree(child2b1.seq, c13.seq, c14.seq), PhylogeneticTree(child2b2.seq, c15.seq, c16.seq))))

class Mutator:
	def __init__(self, seq, structure=None):
		if structure == None:
			structure = getstructure(seq)
		self.seq = seq
		self.structure = structure
		#self.stacks = groupstacks(seq, structure)
		self.triplets = groupTripletBasePairs(seq, structure)
	
	def getResultOfMutations(self, numMutations):
		current = self
		for i in range(numMutations):
			current = current.getResultOfOneMutation()
		return current
		
	
	def getResultOfOneMutation(self):
		randnum = random.random() # between 0 to 1.0
		mutationDist = self.getPossibleTripletMutationDistribution()
		for possibleMutation in mutationDist:
			# probability of mutation, left base index, right base index, what left base changes to, what right base changes to
			randnum -= possibleMutation[0]
			if randnum <= 0:
				break
		idxL = possibleMutation[1]
		idxR = possibleMutation[2]
		targetL = possibleMutation[3]
		targetR = possibleMutation[4]
		newseq = self.seq[:idxL] + targetL + self.seq[idxL+1:idxR] + targetR + self.seq[idxR+1:]
		return Mutator(newseq, self.structure)
	
	'''
	def getResultOfOneMutation(self):
		randnum = random.random() # between 0 to 1.0
		mutationDist = self.getPossibleStackMutationDistribution()
		for possibleMutation in mutationDist:
			randnum -= possibleMutation[0]
			if randnum <= 0:
				break
		parenNum = 0
		for i in range(possibleMutation[1]):
			parenNum += len(self.stacks[i])
		parenNum += possibleMutations[2]
		leftBaseIdx = self.idxOfLeftParen(parenNum)
		rightBaseIdx = self.idxOfRightParen(parenNum)
		newseq = self.seq[:leftBaseIdx] + possibleMutations[3][0] + self.seq[leftBaseIdx+1:rightBaseIdx] + possibleMutations[3][1] + self.seq[rightBaseIdx+1:]
		return Mutator(newseq, self.structure)
	
	def idxOfLeftParen(nthParen):
		for i in range(len(self.seq)):
			if self.seq[i] == "(":
				nthParen -= 1
				if nthParen < 0:
					return i
	
	def idxOfRightParen(nthParen):
		for i in range(len(self.seq)):
			if self.seq[i] == ")":
				nthParen -= 1
				if nthParen < 0:
					return i
	'''
	
	def getPossibleTripletMutationDistribution(self):
		# probability of mutation, left base index, right base index, what left base changes to, what right base changes to
		possible = []
		for triplet in self.triplets: # (left3, paired3) ie AUGCAU -> (( (0,A), (1,U), (2,G) ), ( (5,U), (4,A), (3,C) ))
			leftBP = triplet[0][0][1], triplet[1][0][1]
			assert leftBP in possiblepairs
			centerBP = triplet[0][1][1], triplet[1][1][1]
			assert centerBP in possiblepairs
			rightBP = triplet[0][2][1], triplet[1][2][1]
			assert rightBP in possiblepairs
			idxL = triplet[0][1][0]
			idxR = triplet[1][1][0]
			for target in possiblepairs:
				if centerBP == target:
					continue
				# startPair -> endPair -> leftNeighbor -> rightNeighbor
				probMutation = transitionprobs.probs[centerBP][target][leftBP][rightBP]
				possible.append((probMutation, idxL, idxR, target[0], target[1]))
		# normalize to 1.0:
		totalProb = sum([x[0] for x in possible])
		return [(x[0]/totalProb,) + x[1:] for x in possible]			
		
	'''
	def getPossibleStackMutationDistribution(self):
		# probability of mutation, stack index, index within stack, what it changes to
		possible = []
		for stackNum in range(len(self.stacks)):
			for bpIdx in range(len(self.stacks[stackNum])):
				for target in possiblepairs:
					probMutation = 0
					if bpIdx == 0 and bpIdx == len(self.stacks[stackNum])-1: # no left or right neighbor
						probMutation = math.e**-noLeftOrRightNeighbor[self.stacks[stackNum][bpIdx]][target]
					elif bpIdx == 0: # no left neighbor
						probMutation = math.e**-transitionprobs.noLeftNeighbor[self.stacks[stackNum][bpIdx]][target][self.stacks[stackNum][bpIdx+1]]
					elif bpIdx == len(self.stacks[stackNum])-1: # no right neighbor
						probMutation = math.e**-transitionprobs.noRightNeighbor[self.stacks[stackNum][bpIdx]][target][self.stacks[stackNum][bpIdx-1]]
					else: # have left and right neighbors
						probMutation = math.e**-transitionprobs.neglogprobs[self.stacks[stackNum][bpIdx]][target][self.stacks[stackNum][bpIdx-1]][self.stacks[stackNum][bpIdx+1]]
					possible.append((probMutation, stackNum, bpIdx, target))
		# normalize to 1.0:
		totalProb = sum([x[0] for x in possible])
		return [(x[0]/totalProb,) + x[1:] for x in possible]
	'''

def hammingdist(seq1, seq2):
	total = 0
	for x,y in zip(seq1, seq2):
		if x != y:
			total += 1
	return total

def toStockholm(sequences, structure):
	s = "# STOCKHOLM 1.0\n\n"
	s += "#=GC SS_cons " + (structure.replace('(', '<').replace(')', '>')) + "\n\n"
	for i in range(len(sequences)):
		name = "seq"+str(i)
		s += name+" "+sequences[i]+"\n"
		s += "#=GR "+name+" SS "+structure+"\n"
	return s

'''
structure = "<<<<<<<.<<<<....>>>>........<<<<.<<<<...<<<<...>>>>..<<<<<<<<<<<.....>>>>>>>.>>>>>>>>.>>>>..>>>>>>>."
seq1 = "GCUGAUGUGUCCGAGCGGGCCAAGGCGCGGAACUUCAGGCCCAGAAGCUGGAGAUCCCGCUCGCGUUCAGUGGGCGUGGAUUGAAAUUUCACCAUCAGCG"
seq2 = "GCUGGCGUGUCCGAGCGGGCCAAGGCGCGUCACUUGAGGCCAGGAAGCCUGAGAUCCCGCUCCCGUUCAGGGGGCGUGGAUUCAAAUGACACCGCCAGCG"
seq3 = "GGGGUAGUGAGCGAGCGUUCCAAGGCGCGAGACUUAAGGCCUAGAAGCUGGAGAUUCCGCCCCCGUUCAGGGGGCGUGGAUUUGAAUCUCACCUACCCCG"
seq4 = "GGUGGAGUGAGCGAGCGCUCCAAGGCGCGGGACUUGAGGCCUAGAAGCUAGAGAACCCGUCCCCGUUCAGGGGGCGUGGUUUUAAAUCUCACCUCCGCCG"
print toStockholm([seq1, seq2, seq3, seq4], structure)
'''

def splitLeftRight(s):
	leftParenIdxs = []
	for i in range(len(s)):
		x = s[i]
		if x == "(":
			leftParenIdxs.append(i)
		elif x == ")":
			leftParenIdxs.pop()
		elif x == ",":
			if len(leftParenIdxs) == 0:
				return s[:i], s[i+1:]
	assert False

#print splitLeftRight("(seq1,seq2)node_4,(seq3,seq4)node_5")

def toTree(s, nameToSeq):
	if not "(" in s:
		return nameToSeq[s]
	leftParenIdxs = []
	for i in range(len(s)):
		x = s[i]
		if x == "(":
			leftParenIdxs.append(i)
		elif x == ")":
			startIdx = leftParenIdxs.pop()
			if len(leftParenIdxs) == 0:
				nodeName = s[i+1:]
				nodeSeq = nameToSeq[nodeName]
				interior = s[startIdx+1:i]
				leftPart,rightPart = splitLeftRight(interior)
				return PhylogeneticTree(nodeSeq, toTree(leftPart, nameToSeq), toTree(rightPart, nameToSeq))
	assert False

#print toTree("((seq1,seq2)node_4,(seq3,seq4)node_5)root")
#print toTree("((seq1,seq2)node_4,(seq3,seq4)node_5)root", {"seq1": "36", "seq2": "etye", "seq3": "utm", "seq4": "857", "node_4": "fjk", "node_5": "po", "root": "kfn"})

def fromStockholm(lines):
	seqs = {} # maps from seq name to sequence
	#structure = None
	#for x in lines:
	#	if "#=GC SS_cons" in x:
	#		splitlines = x.split()
	#		structure = splitlines[len(splitlines)-1]
	#		break
	#assert structure != None
	for x in lines:
		if "ancrec_CYK_MAP" in  x and not "ancrec_CYK_MAP_PP" in x:
			splitlines = x.split()
			seqname = splitlines[1]
			sequence = splitlines[len(splitlines)-1]
			sequence = sequence.replace("u", "U").replace("t", "U").replace("a", "A").replace("g", "G").replace("c", "C")
			seqs[seqname] = sequence
	treestructure = None
	for x in lines:
		if "#=GF NH" in x:
			splitlines = x.split()
			treestructure = splitlines[len(splitlines)-1]
			break
	assert treestructure != None
	colonsRemoved = []
	afterColon = False
	for x in treestructure:
		if afterColon:
			if x == "." or x.isdigit():
				continue
			afterColon = False
		if x == ":":
			afterColon = True
			continue
		colonsRemoved.append(x)
	treestructure = "".join(colonsRemoved)
	treestructure = treestructure.replace(" ", "").replace(";", "").strip()
	return toTree(treestructure, seqs)

#f = open("test.stk")
#lines = f.readlines()
#print fromStockholm(lines)

def xratePhylogeny(sequences, structure):
	stockholmFormatInput = toStockholm(sequences, structure)
	process = subprocess.Popen(["./dart/bin/xrate", "-ar", "-nj", "-obl", "-g", "./dart/grammars/rfam.eg"], shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	process.stdin.write(stockholmFormatInput)
	stockholmFormatReconstructed = process.communicate()[0]
	return fromStockholm(stockholmFormatReconstructed.split("\n"))

for numMutations in range(10,31):
	correctTopology = 0
	correctTopologyXrate = 0
	correctTopologyBoth = 0
	totalhammingdist = 0
	totalhammingdistXrate = 0
	totalmutationcost = 0
	totalmutationcostXrate = 0
	leftChildHammingdist = 0
	rightChildHammingdist = 0
	leftChildHammingdistXrate = 0	
	rightChildHammingdistXrate = 0
	leftChildmutationcost = 0
	rightChildmutationcost = 0
	leftChildmutationcostXrate = 0
	rightChildmutationcostXrate = 0
	for x in range(100):
		print >> sys.stderr, numMutations, x
		# >RF00005;tRNA;CP000660.1/666900-666999
		seq = "GCCGGGGUGGCCGAGCGGCCCAAGGCGCGGGACUUGAGGCCCAGAAGCUGGAGAUCCCGUCCCCGUUCAGGGGGCGUGGGUUCAAAUCCCACCCCCGGCG"
		structure =	"(((((((.((((....))))........((((.((((...((((...))))..(((((((((((.....))))))).)))))))).))))..)))))))."
		s = Simulator(seq, structure)
		#tree = s.fourDescendentUltrametric(numMutations)
		#tree = s.eightDescendentUltrametric(numMutations)
		tree = s.sixteenDescendentUltrametric(numMutations)
		#print tree
		leafs = tree.leafs()
		#print leafs
		#print len(leafs)
		tree1 = mutationpath.buildPhylogenySeqs(leafs, structure)
		if leafs == tree1.leafs():
			correctTopology += 1
		tree2 = xratePhylogeny(leafs, structure)
		#print tree2
		if leafs == tree2.leafs():
			correctTopologyXrate += 1
		#if leafs == tree1.leafs() == tree2.leafs():
		totalhammingdistXrate += hammingdist(tree2.node, seq)
		totalhammingdist += hammingdist(tree1.node, seq)
		totalmutationcost += mutationpath.getMutationCostSeq(tree1.node, seq, structure)
		totalmutationcostXrate += mutationpath.getMutationCostSeq(tree2.node, seq, structure)

		if 'PhylogeneticTree' == tree.__class__.__name__ == tree1.__class__.__name__ == tree2.__class__.__name__ == tree.leftChild.__class__.__name__ == tree1.leftChild.__class__.__name__ == tree2.leftChild.__class__.__name__ == tree.rightChild.__class__.__name__ == tree1.rightChild.__class__.__name__ == tree2.rightChild.__class__.__name__:
			leftChildHammingdistXrate += min(hammingdist(tree2.leftChild.node, tree.leftChild.node), hammingdist(tree2.rightChild.node, tree.leftChild.node))
			leftChildHammingdist += min(hammingdist(tree1.leftChild.node, tree.leftChild.node), hammingdist(tree1.rightChild.node, tree.leftChild.node))
			rightChildHammingdistXrate += min(hammingdist(tree2.rightChild.node, tree.rightChild.node), hammingdist(tree2.leftChild.node, tree.rightChild.node))
			rightChildHammingdist += min(hammingdist(tree1.rightChild.node, tree.rightChild.node), hammingdist(tree1.leftChild.node, tree.rightChild.node))

			leftChildmutationcostXrate += min(mutationpath.getMutationCostSeq(tree2.leftChild.node, tree.leftChild.node, structure), mutationpath.getMutationCostSeq(tree2.rightChild.node, tree.leftChild.node, structure))
			leftChildmutationcost += min(mutationpath.getMutationCostSeq(tree1.leftChild.node, tree.leftChild.node, structure), mutationpath.getMutationCostSeq(tree1.rightChild.node, tree.leftChild.node, structure))
			rightChildmutationcostXrate += min(mutationpath.getMutationCostSeq(tree2.rightChild.node, tree.rightChild.node, structure), mutationpath.getMutationCostSeq(tree2.leftChild.node, tree.rightChild.node, structure))
			rightChildmutationcost += min(mutationpath.getMutationCostSeq(tree1.rightChild.node, tree.rightChild.node, structure), mutationpath.getMutationCostSeq(tree1.leftChild.node, tree.rightChild.node, structure))

		#correctTopologyBoth += 1
		#print tree1.leafs()
		#print leafs
		#print leafs
		#print tree1
		#print len(tree1.leafs())
		#print tree.leafs()
		#print leafs
		#for node2 in leafs):
		#	#print node1
		#	print node2
		#	if node1 == node2:
		#		print "ping"
		#		correctTopology += 1
	print "number mutations", numMutations
	print "number with correct topology", correctTopology, correctTopologyXrate
	print "hamming dist", totalhammingdist, totalhammingdistXrate
	print "mutation cost", totalmutationcost, totalmutationcostXrate
	print "hamming dist left child", leftChildHammingdist, leftChildHammingdistXrate
	print "hamming dist right child", rightChildHammingdist, rightChildHammingdistXrate
	print "mutation cost left child", leftChildmutationcost, leftChildmutationcostXrate
	print "mutation cost right child", rightChildmutationcost, rightChildmutationcostXrate

#m = Mutator("GCCGACGGUCAUAGGACGGGGGAAACACCCGGACUCAUUCCGAACCCGGAAGUUAAGCCCCGUUCCGUCCCGCACAGUACUGUGUUCCGAGAGGGCACGGGAACUGCGGGAACCGUCGGCU", "(((((((((....((((((((..(((..((((..((.....))..))))..)))...))))))))..(((((((.....((((((.((....))))))))....)))))))))))))))).")
#m2 = m.getResultOfOneMutation()
#print m.seq
#print m2.seq
