#!/usr/bin/python

import heapq
import transitionprobs
import math
import subprocess
import sys
import random
import itertools

class GraphNode:
	def __init__(self, seq, cost, heuvpluscost, predecessor):
		self.seq = seq
		self.cost = cost
		self.heuvpluscost = heuvpluscost
		self.predecessor = predecessor
	
	def getPath(self):
		path = []
		current = self
		while current != None:
			path.append(current)
			current = current.predecessor
		return path
	
	def printtraceback(self):
		current = self
		while current != None:
			print current.seq, current.cost
			current = current.predecessor
	
	def __eq__(self, other):
		if other == None:
			return False
		return self.seq == other.seq
	
	def __ne__(self, other):
		return not self.__eq__(other)
	
	def __hash__(self):
		return hash(self.seq)
	
	def __cmp__(self, other):
		return cmp(self.heuvpluscost, other.heuvpluscost)

def getheuv(startseq, targetseq):
	#return 0.0 # disable A*
	total = 0.0
	for startPair,endPair in zip(startseq, targetseq):
		if startPair != endPair:
			total += transitionprobs.heuristic[startPair][endPair]
	return total

possiblepairs = [("A", "U"), ("U", "A"), ("C", "G"), ("G", "C"), ("G", "U"), ("U", "G")]
'''
def getEdges(node, targetseq):
	assert len(node.seq) == len(targetseq)
	edges = []
	for i in range(len(node.seq)):
		if node.seq[i] == targetseq[i]:
			continue
		newedges = possiblepairs[:]
		newedges.remove(node.seq[i])
		for x in newedges:
			newseq = tuple(list(node.seq[:i])+list([x])+list(node.seq[i+1:]))
			heuv = getheuv(newseq, targetseq)
			newcost = node.cost
			if i == 0 and i == len(node.seq)-1: # no left or right neighbor
				newcost += transitionprobs.noLeftOrRightNeighbor[node.seq[i]][targetseq[i]]
			elif i == 0: # no left neighbor
				newcost += transitionprobs.noLeftNeighbor[node.seq[i]][targetseq[i]][node.seq[i-1]]
			elif i == len(node.seq)-1: # no right neighbor
				newcost += transitionprobs.noRightNeighbor[node.seq[i]][targetseq[i]][node.seq[i+1]]
			else: # have left and right neighbor
				newcost += transitionprobs.costs[node.seq[i]][targetseq[i]][node.seq[i-1]][node.seq[i+1]]
			edges.append(GraphNode(newseq, newcost, newcost+heuv, node))
			#edges.append(node.seq[:i]+[x]+node[i+1]:, cost+1, cost+1, node) # todo change cost+1 to use transitions matrix, todo use heuv for A*
	return edges
'''

def getEdges(node, targetseq):
	assert len(node.seq) == len(targetseq)
	edges = []
	for i in range(len(node.seq)):
		if node.seq[i] == targetseq[i]:
			continue
		newedges = possiblepairs[:]
		newedges.remove(node.seq[i])
		for x in newedges:
			newseq = tuple(list(node.seq[:i])+list([x])+list(node.seq[i+1:]))
			heuv = getheuv(newseq, targetseq)
			newcost = node.cost
			if i == 0 and i == len(node.seq)-1: # no left or right neighbor
				newcost += transitionprobs.noLeftOrRightNeighbor[node.seq[i]][targetseq[i]]
			elif i == 0: # no left neighbor
				newcost += transitionprobs.noLeftNeighbor[node.seq[i]][x][node.seq[i+1]]
			elif i == len(node.seq)-1: # no right neighbor
				newcost += transitionprobs.noRightNeighbor[node.seq[i]][x][node.seq[i-1]]
			else: # have left and right neighbor
				newcost += transitionprobs.costs[node.seq[i]][x][node.seq[i-1]][node.seq[i+1]]
			edges.append(GraphNode(newseq, newcost, newcost+heuv, node))
			#edges.append(node.seq[:i]+[x]+node[i+1]:, cost+1, cost+1, node) # todo change cost+1 to use transitions matrix, todo use heuv for A*
	return edges

def getEdgesNoHeu(node):
	edges = []
	for i in range(len(node.seq)):
		newedges = possiblepairs[:]
		newedges.remove(node.seq[i])
		for x in newedges:
			newseq = tuple(list(node.seq[:i])+list([x])+list(node.seq[i+1:]))
			newcost = node.cost
			heuv = 0.0
			if i == 0 and i == len(node.seq)-1: # no left or right neighbor
				newcost += transitionprobs.noLeftOrRightNeighbor[node.seq[i]][x]
			elif i == 0: # no left neighbor
				newcost += transitionprobs.noLeftNeighbor[node.seq[i]][x][node.seq[i+1]]
			elif i == len(node.seq)-1: # no right neighbor
				newcost += transitionprobs.noRightNeighbor[node.seq[i]][x][node.seq[i-1]]
			else: # have left and right neighbor
				newcost += transitionprobs.costs[node.seq[i]][x][node.seq[i-1]][node.seq[i+1]]
			edges.append(GraphNode(newseq, newcost, newcost+heuv, node))
			#edges.append(node.seq[:i]+[x]+node[i+1]:, cost+1, cost+1, node) # todo change cost+1 to use transitions matrix, todo use heuv for A*
	return edges

def getMutationCostSeq(startSeq, targetSeq, structure):
	assert len(startSeq) == len(targetSeq)
	startStackRegions = groupstacks(startSeq, structure)
	targetStackRegions = groupstacks(targetSeq, structure)
	assert len(startStackRegions) == len(targetStackRegions)
	return sum([getMutationPath(x,y).cost for x,y in zip(startStackRegions,targetStackRegions)])

# start: [(A, U), (G, C)]...
def getMutationPath(start, targetseq):
	seqlen = len(start)
	for x in start:
		if x not in possiblepairs:
			print "non-base pair present in start"
	initial = GraphNode(start, 0, getheuv(start, targetseq), None)
	if start == targetseq:
		return initial
	visited = {start: True}
	agenda = [initial] # heuv+cost, cost, seqence, predecessor
	while True:
		current = heapq.heappop(agenda)
		newedges = getEdges(current, targetseq)
		for newnode in newedges:
			if newnode.seq in visited:
				continue
			if newnode.seq == targetseq: # reached target
				return newnode
			visited[newnode.seq] = True
			heapq.heappush(agenda, newnode)

def minidx(agendas):
	mini = 0
	for i in range(len(agendas)):
		if len(agendas[i]) > 0 and agendas[i][0] < agendas[mini][0]:
			mini = i
	return mini

def alltrue(l):
	for x in l:
		if not x:
			return False
	return True

def getCenter(stacks):
	agendas = [[GraphNode(x, 0, 0, None)] for x in stacks]
	visited = [{stacks[i]: agendas[i][0]} for i in range(len(stacks))]
	while True:
		idx = minidx(agendas)
		current = heapq.heappop(agendas[idx])
		newedges = getEdgesNoHeu(current)
		for newnode in newedges:
			if newnode.seq in visited[idx]:
				continue
			visited[idx][newnode.seq] = newnode
			heapq.heappush(agendas[idx], newnode)
			if alltrue([newnode.seq in visited[i] for i in range(len(stacks))]): # connects all
				return [visited[i][newnode.seq] for i in range(len(stacks))]

def splitinto2groupsOld(full, groupsize):
	sf = set(full)
	halfgroups = set(itertools.combinations(sf,groupsize))
	return [(list(set(x)), list(sf-set(x))) for x in halfgroups]

def splitinto2groups(full, groupsize):
	halfgroups = itertools.combinations(full,groupsize)
	rv = []
	for x in halfgroups:
		half1 = x
		half2 = full[:]
		for v in half1:
			half2.remove(v)
		rv.append((list(half1), list(half2)))
	return rv

class PhylogeneticTree:
	def __init__(self, node, child1, child2):
		self.node = node
		if child1.__class__.__name__ == 'PhylogeneticTree' and child2.__class__.__name__ == 'PhylogeneticTree':
			if child1.leftMostLeafNode < child2.leftMostLeafNode:
				self.leftChild = child1
				self.rightChild = child2
				self.leftMostLeafNode = child1.leftMostLeafNode
			else:
				self.leftChild = child2
				self.rightChild = child1
				self.leftMostLeafNode = child2.leftMostLeafNode
		elif child1.__class__.__name__ == 'PhylogeneticTree':
			if child1.leftMostLeafNode < child2:
				self.leftChild = child1
				self.rightChild = child2
				self.leftMostLeafNode = child1.leftMostLeafNode
			else:
				self.leftChild = child2
				self.rightChild = child1
				self.leftMostLeafNode = child2
		elif child2.__class__.__name__ == 'PhylogeneticTree':
			if child1 < child2.leftMostLeafNode:
				self.leftChild = child1
				self.rightChild = child2
				self.leftMostLeafNode = child1
			else:
				self.leftChild = child2
				self.rightChild = child1
				self.leftMostLeafNode = child2.leftMostLeafNode
		else:
			if child1 < child2:
				self.leftChild = child1
				self.rightChild = child2
				self.leftMostLeafNode = child1
			else:
				self.leftChild = child2
				self.rightChild = child1
				self.leftMostLeafNode = child2
	
	def leafs(self):
		l = []
		if self.leftChild.__class__.__name__ == 'PhylogeneticTree':
			l += self.leftChild.leafs()
		else:
			l += [self.leftChild]
		if self.rightChild.__class__.__name__ == 'PhylogeneticTree':
			l += self.rightChild.leafs()
		else:
			l += [self.rightChild]
		return l
	
	def __str__(self):
		return str(self.node) + '\n{\n' + str(self.leftChild) + '\n|\n' + str(self.rightChild) + '\n}'

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
				allStacks.append(tuple(currentStack))
				currentStack = []
			lIdx = x[0]
			rIdx = x[1]
			currentStack.append((x[2],x[3]))
	if len(currentStack) > 0:
		allStacks.append(tuple(currentStack))
	return allStacks

def replacestacks(sequence, structure, newstacks):
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
		return sequence
	basepairs.sort()
	allStacks = []
	lIdx = -9
	rIdx = -9
	currentStack = []
	for x in basepairs: # (leftidx, rightidx, leftbase, rightbase)
		if x[0] == lIdx+1 and x[1] == rIdx-1:
			lIdx = x[0]
			rIdx = x[1]
			currentStack.append(x)
		else:
			if len(currentStack) > 0:
				allStacks.append(tuple(currentStack))
				currentStack = []
			lIdx = x[0]
			rIdx = x[1]
			currentStack.append(x)
	if len(currentStack) > 0:
		allStacks.append(tuple(currentStack))
	assert len(allStacks) == len(newstacks)
	newseq = sequence[:]
	for j in range(len(allStacks)):
		currentStackGroup = allStacks[j]
		replaceWithStackGroup = newstacks[j]
		assert len(currentStackGroup) == len(replaceWithStackGroup)
		for i in range(len(currentStackGroup)):
			#print currentStackGroup[i]
			leftidx,rightidx,origleftbase,origrightbase = currentStackGroup[i]
			leftbase,rightbase = replaceWithStackGroup[i]
			newseq = newseq[:leftidx] + leftbase + newseq[leftidx+1:rightidx] + rightbase + newseq[rightidx+1:]
	return newseq

#q = groupstacks("AUCGGAUGUG", "(((.))).()")
#newstacks = [(('G', 'U'), ('U', 'G'), ('C', 'G')), (('U', 'A'),)]
#print replacestacks("AUCGGAUGUG", "(((.))).()", newstacks)
#AUCGGAUGUG
#AUCGGAUGUG
#AUCGGGUGUA
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

def groupSplitCostPathsSeq(half1, half2, ancestor, pairwisedists): # pairwisedists is actually structure
	assert len(pairwisedists) == len(ancestor)
	half1a = half1 + [ancestor]
	half2a = half2 + [ancestor]
	half1Stacks = [groupstacks(x,pairwisedists) for x in half1a]
	half2Stacks = [groupstacks(x,pairwisedists) for x in half2a]
	half1StacksByType = []
	for i in range(len(half1Stacks[0])):
		currentStack = [x[i] for x in half1Stacks]
		half1StacksByType.append(currentStack)
	half2StacksByType = []
	for i in range(len(half2Stacks[0])):
		currentStack = [x[i] for x in half2Stacks]
		half2StacksByType.append(currentStack)
	half1ancestorByType = [getCenter(x) for x in half1StacksByType]
	half2ancestorByType = [getCenter(x) for x in half2StacksByType]
	h1cost = sum([sum([y.cost for y in z]) for z in half1ancestorByType])
	h2cost = sum([sum([y.cost for y in z]) for z in half2ancestorByType])
	return h1cost + h2cost

def groupSplitCostPaths(half1, half2, ancestor, pairwisedists):
	half1a = half1 + [ancestor]
	half2a = half2 + [ancestor]
	half1ancestor = getCenter(half1a)
	half2ancestor = getCenter(half2a)
	h1cost = sum([x.cost for x in half1ancestor])
	h2cost = sum([x.cost for x in half2ancestor])
	return h1cost + h2cost

def groupSplitCostPairsSum(half1, half2, ancestor, pairwisedists):
	half1cost = sum([pairwisedists[x][y] for x in half1 for y in half1])
	half2cost = sum([pairwisedists[x][y] for x in half2 for y in half2])
	return half1cost + half2cost

def buildPhylogenySeqs(sequences, structure):
	stackGroups = [groupstacks(x,structure) for x in sequences]
	lengths = [len(y) for y in stackGroups[0]]
	for x in stackGroups:
		assert lengths == [len(y) for y in x]
	pairwisedists = {} # seq -> seq -> -logprob path
	for x in sequences:
		if not x in pairwisedists:
			pairwisedists[x] = {}
		for y in sequences:
			xstacks = groupstacks(x,structure)
			ystacks = groupstacks(y,structure)
			pairwisedists[x][y] = sum([getMutationPath(z,w).cost for z,w in zip(xstacks,ystacks)])
	return buildPhylogenyGroupingHeuSeqs(pairwisedists, sequences, structure, groupSplitCostPairsSum, [])

def buildPhylogenySeqsSlow(sequences, structure):
	return buildPhylogenyGroupingHeuSeqs(structure, sequences, structure, groupSplitCostPathsSeq, [])

def buildPhylogenyGroupingHeuSeqs(pairwisedists, sequences, structure, groupSplitfn, higherAncestor):
	if len(sequences) == 1:
		return sequences[0]
	if len(sequences) == 0:
		return None
	sequencesWithAdded = sequences + higherAncestor
	sequencesWithAddedStacks = [groupstacks(x,structure) for x in sequencesWithAdded]
	allStacks = []
	for i in range(len(sequencesWithAddedStacks[0])):
		currentStack = [x[i] for x in sequencesWithAddedStacks]
		allStacks.append(currentStack)
	ancestorPathsForStacks = [getCenter(stackGroup) for stackGroup in allStacks]
	ancestorStacks = [ancestorPathsForStacks[i][0].seq for i in range(len(ancestorPathsForStacks))]
	ancestor = replacestacks(sequences[0], structure, ancestorStacks) # sequence
	groups = splitinto2groups(sequences, len(sequences)/2)
	bhalf1 = None
	bhalf2 = None
	lowestcost = sys.maxint
	for half1, half2 in groups:
		assert len(half1) == len(half2) == len(sequences)/2
		cost = groupSplitfn(half1, half2, ancestor, pairwisedists)
		if cost < lowestcost:
			lowestcost = cost
			bhalf1 = half1
			bhalf2 = half2
	child1 = buildPhylogenyGroupingHeuSeqs(pairwisedists, bhalf1, structure, groupSplitfn, [ancestor])
	child2 = buildPhylogenyGroupingHeuSeqs(pairwisedists, bhalf2, structure, groupSplitfn, [ancestor])
	return PhylogeneticTree(ancestor, child1, child2)

#q = buildPhylogenySeqsSlow(["AGGUU", "GGGUU", "AAGUU", "AAGUU"], "((.))")
#print q

def buildPhylogeny(stacks):
	pairwisedists = {} # stack -> stack -> -logprob path
	for x in stacks:
		if not x in pairwisedists:
			pairwisedists[x] = {}
		for y in stacks:
			pairwisedists[x][y] = getMutationPath(x,y).cost
	return buildPhylogenyGroupingHeu(pairwisedists, stacks, groupSplitCostPairsSum, [])

def buildPhylogenySlow(stacks, higherAncestor=[]):
	return buildPhylogenyGroupingHeu(None, stacks, groupSplitCostPaths, [])

def buildPhylogenyGroupingHeu(pairwisedists, stacks, groupSplitfn, higherAncestor):
	if len(stacks) == 1:
		return stacks[0]
	if len(stacks) == 0:
		return None
	stacksWithAdded = stacks + higherAncestor
	ancestorPaths = getCenter(stacksWithAdded)
	print "getCenter"
	ancestor = ancestorPaths[0].seq
	ancestorPaths = getCenter(stacksWithAdded)
	groups = splitinto2groups(stacks, len(stacks)/2)
	bhalf1 = None
	bhalf2 = None
	lowestcost = sys.maxint
	for half1, half2 in groups:
		cost = groupSplitfn(half1, half2, ancestor, pairwisedists)
		if cost < lowestcost:
			lowestcost = cost
			bhalf1 = half1
			bhalf2 = half2
	child1 = buildPhylogenyGroupingHeu(pairwisedists, bhalf1, groupSplitfn, [ancestor])
	child2 = buildPhylogenyGroupingHeu(pairwisedists, bhalf2, groupSplitfn, [ancestor])
	return PhylogeneticTree(ancestor, child1, child2)

'''
def buildPhylogenySlow(stacks, higherAncestor=[]):
	if len(stacks) == 1:
		return stacks[0]
	if len(stacks) == 0:
		return None
	stacksWithAdded = stacks + higherAncestor
	ancestorPaths = getCenter(stacksWithAdded)
	print "getCenter"
	ancestor = ancestorPaths[0].seq
	ancestorPaths = getCenter(stacksWithAdded)
	groups = splitinto2groups(stacks, len(stacks)/2)
	bhalf1 = None
	bhalf2 = None
	lowestcost = sys.maxint
	for half1, half2 in groups:
		half1a = half1 + [ancestor]
		half2a = half2 + [ancestor]
		half1ancestor = getCenter(half1a)
		half2ancestor = getCenter(half2a)
		h1cost = sum([x.cost for x in half1ancestor])
		h2cost = sum([x.cost for x in half2ancestor])
		if h1cost + h2cost < lowestcost:
			lowestcost = h1cost + h2cost
			bhalf1 = half1
			bhalf2 = half2
	child1 = buildPhylogenySlow(bhalf1, [ancestor])
	child2 = buildPhylogenySlow(bhalf2, [ancestor])
	return PhylogeneticTree(ancestor, child1, child2)
'''
'''
	ancestor = ancestorPaths[0].seq
	agenda = [(ancestor, stacks)]
	
	while True:
		
	
	ancestorPaths = getCenter(stacks)
	ancestor = ancestorPaths[0].seq
	groups = splitinto2groups(stacks, len(stacks)/2)
	bhalf1 = None
	bhalf2 = None
	lowestcost = sys.maxint
	for half1, half2 in groups:
		half1a = half1 + [ancestor]
		half2a = half2 + [ancestor]
		half1ancestor = getCenter(half1a)
		half2ancestor = getCenter(half2a)
		h1cost = sum([x.cost for x in half1ancestor])
		h2cost = sum([x.cost for x in half2ancestor])
		if h1cost + h2cost < lowestcost:
			lowestcost = h1cost + h2cost
			bhalf1 = 
'''
		

def getstructure(seq):
	process = subprocess.Popen(["./RNAfold", "-noPS"], shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	process.stdin.write(seq)
	return process.communicate()[0].split()[1]

def getMutationPathSeqs(startSeq, targetSeq, structure=None):
	assert len(startSeq) == len(targetSeq)
	if structure == None:
		structure = getstructure(startSeq)
	startStacks = groupstacks(startSeq, structure)
	targetStacks = groupstacks(targetSeq, structure)
	assert len(startStacks) == len(targetStacks)
	print startStacks
	print targetStacks
	return [getMutationPath(x,y) for x,y in zip(startStacks,targetStacks)]

#def makePhylogeny():
#	
'''
def neighborJoin4Seqs(seq1, seq2, seq3, seq4):
	path1 = 

def neighborJoin(path1, path2):
	intermediateNodes = [[None for x in range(len(path2))] for y in range(len(path1))]
	bestx = 0
	besty = 0
	bestcost = int.maxint
	for y in range(len(path1)):
		for x in range(len(path2)):
			intermediateNodes[y][x] = getMutationPath(path1[y], path2[x])
'''

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

def prod(l):
	total = 1
	for x in l:
		total *= x
	return total

'''
seq1 = "CAAGCCGGCCAUAGCGUCAGGGUGCGACCCAAUCCCAUUCCGAACUUGGAAGUCAAACCUGAUGUCGCUGUUGUGUUACUAAGAUGCGAGAGCUCUUGGGAAGCAACAGUGCUGGCAUCA"
stu1 = "...((((((....((((((((....(((.........(((((....))))))))...)))))))).((((((((....((((((.((....))))))))...))))))))))))))...."
stacks1 = groupstacks(seq1, stu1)

seq2 = "GGCCCGGCCAUAGCGGCCGGGUAACACCCGGACUCAUUUCGAACCCGGAAGUUAAGCCGGCCGCGUUGGAGGCUCCAGUGGGGUCCGAGAGGCCCUGCAGGGGCCUCCAAGCCGGGGCCG"
stu2 = "((((((((....((((((((.((((..((((..((.....))..))))..))))..))))))))..((((((((((.((((((.((....)))))))).)))))))))).)))))).))."
stacks2 = groupstacks(seq2, stu2)

assert len(stacks1) == len(stacks2)
getMutationPath(x,y) for x,y in zip(,)
'''
'''
startseq = (("A", "U"), ("U", "A"), ("G", "C"))
endseq = (("U", "A"), ("C", "G"), ("G", "C"))
end = getMutationPath(startseq, endseq)
end.printtraceback()
'''
'''
def getEdges(node, cost, target, transitions):
	assert len(node) == len(target)
	edges = []
	for i in range(len(node)):
		if node[i][0] == target[i][0] and node[i][1] == target[i][1]:
			continue
		newedges = possiblepairs[:]
		newedges.remove(node[i])
		for x in newedges:
			edges.append(node[:i]+[x]+node[i+1]:, cost+1, cost+1, node) # todo change cost+1 to use transitions matrix, todo use heuv for A*
	return edges

# start: [(A, U), (G, C)]...
def getMutationPath(start, target, transitions):
	seqlen = len(start)
	for x in start:
		if x not in possiblepairs:
			print "non-base pair present in start"
	visited = {start: True}
	agenda = [(0, 0, start, None)] # heuv+cost, cost, seqence, predecessor
	while True:
		heuvcost,cost,node,prev = heapq.heappop(agenda)
		newedges = getEdges(node, cost, target, transitions)
		for newnode,newcost,heuvpluscost,predecessor in newedges:
			if newnode in visited:
				continue
			if newnode == target: # reached target
				return newnode,newcost,predecessor
			visited[newnode] = True
			heapq.heappush(agenda, (heuvpluscost,newcost,newnode,predecessor))
'''
'''
seq1 = "GCCGAUGGUCAUAGGACGGGGGAAACACCCGGACUCAUUCCGAACCCGGAAGUUAAGCCCCGUCCCGUCCUGCACAGUACUGUGUUCCGAGAGGGCACGGGAACUGCAGGAACCAUCGGCU"
seq2 = "GCAGACGGUCAUAGGACGGUGGAAACACCCGGACUCAUUCCGAACCCGGAAGUUAAGCACCGUUCCGUCCCACACAGUACCGUGUUCCGAGAGGGCACGGGAACUGUGGGAACCGUCUGCU"
structure = "(((((((((....((((((((..(((..((((..((.....))..))))..)))...))))))))..(((((((.....((((((.((....))))))))....))))))))))))))))."
for x in getMutationPathSeqs(seq1, seq2):
	x.printtraceback()
'''
'''
#s1 = (('G', 'C'), ('C', 'G'), ('C', 'G'), ('G', 'C'), ('A', 'U'), ('U', 'A'), ('G', 'C'), ('G', 'C'), ('U', 'A'))
#s2 = (('G', 'C'), ('C', 'G'), ('A', 'U'), ('G', 'C'), ('A', 'U'), ('C', 'G'), ('G', 'C'), ('G', 'C'), ('U', 'A'))

#s1 = (('G', 'C'), ('G', 'C'), ('A', 'U'), ('C', 'G'), ('G', 'C'), ('G', 'C'), ('G', 'C'), ('G', 'C'))
#s2 = (('G', 'C'), ('G', 'U'), ('A', 'U'), ('C', 'G'), ('G', 'C'), ('G', 'C'), ('U', 'A'), ('G', 'C'))

#s1 = (('A', 'U'), ('A', 'U'), ('C', 'G'))
#s2 = (('A', 'U'), ('A', 'U'), ('C', 'G'))

#path = getMutationPath(s1, s2)
#path.printtraceback()
s1 = (('G', 'C'), ('G', 'C'), ('A', 'U'), ('C', 'G'), ('G', 'C'), ('G', 'C'), ('G', 'C'), ('G', 'C'))
s2 = (('G', 'C'), ('G', 'C'), ('A', 'U'), ('A', 'U'), ('G', 'C'), ('G', 'C'), ('G', 'C'), ('G', 'C'))
s3 = (('G', 'C'), ('G', 'C'), ('A', 'U'), ('G', 'U'), ('G', 'C'), ('G', 'C'), ('G', 'C'), ('G', 'C'))
s4 = (('G', 'U'), ('G', 'C'), ('A', 'U'), ('G', 'U'), ('G', 'C'), ('G', 'C'), ('G', 'C'), ('G', 'C'))
s5 = (('G', 'U'), ('G', 'C'), ('A', 'U'), ('G', 'U'), ('G', 'U'), ('G', 'C'), ('G', 'C'), ('G', 'C'))
s6 = (('G', 'U'), ('G', 'C'), ('A', 'U'), ('G', 'U'), ('G', 'U'), ('G', 'C'), ('G', 'U'), ('G', 'C'))
s7 = (('G', 'U'), ('G', 'C'), ('A', 'U'), ('G', 'U'), ('G', 'U'), ('G', 'C'), ('G', 'U'), ('G', 'U'))
s8 = (('G', 'U'), ('G', 'U'), ('A', 'U'), ('G', 'U'), ('G', 'U'), ('G', 'C'), ('G', 'U'), ('G', 'U'))
#s1 = (('A', 'U'), ('C', 'G'), ('G', 'C'))
#s2 = (('A', 'U'), ('A', 'U'), ('G', 'C'))

q = buildPhylogeny([s1,s2,s3,s4,s5,s6,s7,s8])
print q
print q.leafs()

#q = getCenter([s1, s2, s3])
#for x in q:
#	x.printtraceback()
'''
