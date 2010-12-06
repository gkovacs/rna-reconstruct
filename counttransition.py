#!/usr/bin/python

import subprocess
import sys
import pickle

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

def allfreqs3():
	bases = ["A", "C", "G", "U"]
	retv = []
	for x in bases:
		for y in bases:
			for z in bases:
				for i in bases:
					for j in bases:
						for k in bases:
							retv.append(x+y+z+"|"+i+j+k)
	return retv

def allfreqs2():
	bases = ["A", "C", "G", "U"]
	retv = []
	for x in bases:
		for y in bases:
			for i in bases:
				for j in bases:
					retv.append(x+y+"|"+i+j)
	return retv

def allfreqs1():
	bases = ["A", "C", "G", "U"]
	retv = []
	for x in bases:
		for y in bases:
			retv.append(x+"|"+y)
	return retv

def countstackfreqs3(stacks):
	freqs = {}
	for x in allfreqs3():
		freqs[x] = 0
	for x in stacks:
		for i in range(1,len(x)-1):
			side1 = x[i-1][0] + x[i][0] + x[i+1][0]
			side2 = x[i-1][1] + x[i][1] + x[i+1][1]
			freqs[side1+"|"+side2] += 1
	return freqs

def countstackfreqs2(stacks):
	freqs = {}
	for x in allfreqs2():
		freqs[x] = 0
	for x in stacks:
		for i in range(1,len(x)):
			side1 = x[i-1][0] + x[i][0]
			side2 = x[i-1][1] + x[i][1]
			freqs[side1+"|"+side2] += 1
	return freqs

def countstackfreqs1(stacks):
	freqs = {}
	for x in allfreqs1():
		freqs[x] = 0
	for x in stacks:
		for i in range(len(x)):
			side1 = x[i][0]
			side2 = x[i][1]
			freqs[side1+"|"+side2] += 1
	return freqs

def hammingdist(seq1, seq2):
	assert len(seq1) == len(seq2)
	total = 0
	for x,y in zip(seq1,seq2):
		if x != y:
			total += 1
	return total

def minidx(l):
	minidx = 0
	for i in range(len(l)):
		if l[i] < l[minidx]:
			minidx = i
	return minidx

def minidxnonzero(l):
	minidx = 0
	for i in range(len(l)):
		if l[minidx] > 0 and l[i] < l[minidx]:
			minidx = i
	return minidx

possiblepairs = [("A", "U"), ("U", "A"), ("C", "G"), ("G", "C"), ("G", "U"), ("U", "G")]

def emptyTransitionFreqMatrix():
	m = {}
	for startPair in possiblepairs:
		if not startPair in m:
			m[startPair] = {}
		for endPair in possiblepairs:
			if not endPair in m[startPair]:
				m[startPair][endPair] = {}
			for leftNeighbor in possiblepairs:
				if not leftNeighbor in m[startPair][endPair]:
					m[startPair][endPair][leftNeighbor] = {}
				for rightNeighbor in possiblepairs:
					m[startPair][endPair][leftNeighbor][rightNeighbor] = 0
	return m

totalcounted = 0

def countTransitionFreqs(stackRegions1, stackRegions2):
	assert len(stackRegions1) == len(stackRegions2)
	transitions = emptyTransitionFreqMatrix()
	for stacks1,stacks2 in zip(stackRegions1,stackRegions2):
		assert len(stacks1) == len(stacks2)
		for i in range(1,len(stacks1)-1):
			leftStack1 = stacks1[i-1]
			leftStack2 = stacks2[i-1]
			rightStack1 = stacks1[i+1]
			rightStack2 = stacks2[i+1]
			if leftStack1 != leftStack2 or rightStack1 != rightStack2:
				continue
			before = stacks1[i]
			after = stacks2[i]
			if before == after:
				continue
			if not before in possiblepairs:
				continue
			if not after in possiblepairs:
				continue
			if not leftStack1 in possiblepairs:
				continue
			if not rightStack1 in possiblepairs:
				continue
			transitions[before][after][leftStack1][rightStack1] += 1

			global totalcounted
			totalcounted += 1
			#print "counted", totalcounted
	return transitions

def countTransitionFreqsSingleStack(stacks1,stacks2):
	transitions = emptyTransitionFreqMatrix()
	assert len(stacks1) == len(stacks2)
	for i in range(1,len(stacks1)-1):
		leftStack1 = stacks1[i-1]
		leftStack2 = stacks2[i-1]
		rightStack1 = stacks1[i+1]
		rightStack2 = stacks2[i+1]
		if leftStack1 != leftStack2 or rightStack1 != rightStack2:
			continue
		before = stacks1[i]
		after = stacks2[i]
		#if before == after:
		#	continue
		if not before in possiblepairs:
			continue
		if not after in possiblepairs:
			continue
		if not leftStack1 in possiblepairs:
			continue
		if not rightStack1 in possiblepairs:
			continue
		transitions[before][after][leftStack1][rightStack1] += 1
		global totalcounted
		totalcounted += 1
		#print "counted", totalcounted
	return transitions

def sumTransitionFreqs(ltransitions):
	total = emptyTransitionFreqMatrix()
	for x in ltransitions:	
		for startPair in x:
			for endPair in x[startPair]:
				for leftNeighbor in x[startPair][endPair]:
					for rightNeighbor in x[startPair][endPair][leftNeighbor]:
						total[startPair][endPair][leftNeighbor][rightNeighbor] += x[startPair][endPair][leftNeighbor][rightNeighbor]
	return total

def addTransitionFreqsToExisting(existing, x):
	for startPair in x:
		for endPair in x[startPair]:
			for leftNeighbor in x[startPair][endPair]:
				for rightNeighbor in x[startPair][endPair][leftNeighbor]:
					existing[startPair][endPair][leftNeighbor][rightNeighbor] += x[startPair][endPair][leftNeighbor][rightNeighbor]

fn = sys.argv[1] # "singlefromfamily-structure.fasta"
f = open(fn)
conts = f.readlines()
families = [{} for i in range(100000)] # families -> stacklength -> stack,sequence
for x in range(1,len(conts),3):
	name = conts[x-1]
	rfnum = int(name[3:8])
	sequence = conts[x]
	structure = conts[x+1]
	stacks = groupstacks(sequence, structure)
	for stack in stacks:
		if not len(stack) in families[rfnum]:
			families[rfnum][len(stack)] = []
		families[rfnum][len(stack)].append((stack,sequence))

numseqscounted = 0
#transitions = []
sumtransitions = emptyTransitionFreqMatrix()
for x in range(len(families)):
	stacknum = 0
	print "family ", x
	print sumtransitions
	for stackLength in families[x]:
		stacksAndSeqs = families[x][stackLength]
		print "stack ", stacknum, "of", len(families[x]), "with", len(stacksAndSeqs)
		stacknum += 1
		distances = [[None for q in stacksAndSeqs] for z in stacksAndSeqs]
		for i in range(len(stacksAndSeqs)):
			for j in range(len(stacksAndSeqs)):
				if j == i:
					distances[i][j] = sys.maxint
				elif len(stacksAndSeqs[i][1]) != len(stacksAndSeqs[j][1]): # different sequence lengths
					distances[i][j] = sys.maxint
				else:
					distances[i][j] = hammingdist(stacksAndSeqs[i][1], stacksAndSeqs[j][1])
					#distances[i][j] = hammingdist(stacksAndSeqs[i][0], stacksAndSeqs[j][0]) # comparing stacks
					#distances[i][j] = hammingdist(stacksAndSeqs[i][1], stacksAndSeqs[j][1]) # comparing sequences
		for i in range(len(stacksAndSeqs)):
			
			j = minidx(distances[i]) # find the one with the least distance to seq[i]
			val = distances[i][j] # hamming distance to the closest match
			if val < 0.5*len(stacksAndSeqs[i][1]): # more than 90% sequence identity, count its transitions
				numseqscounted += 1
				addTransitionFreqsToExisting(sumtransitions, countTransitionFreqsSingleStack(stacksAndSeqs[i][0], stacksAndSeqs[j][0]))

print "number of seqs counted", numseqscounted
print "total number of transitions counted", totalcounted
#sumtransitions = sumTransitionFreqs(transitions)
output = open(fn+'.pkl', 'wb')
pickle.dump(sumtransitions, output)

