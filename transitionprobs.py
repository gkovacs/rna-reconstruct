#!/usr/bin/python

from transitionfreqs import freqs
import math

possiblepairs = [("A", "U"), ("U", "A"), ("C", "G"), ("G", "C"), ("G", "U"), ("U", "G")]

# startPair -> endPair -> leftNeighbor -> rightNeighbor
# make transitions symmetric
sfreqs = {}
for startPair in possiblepairs:
	if not startPair in sfreqs:
		sfreqs[startPair] = {}
	for endPair in possiblepairs:
		if not endPair in sfreqs[startPair]:
			sfreqs[startPair][endPair] = {}
		for leftNeighbor in possiblepairs:
			if not leftNeighbor in sfreqs[startPair][endPair]:
				sfreqs[startPair][endPair][leftNeighbor] = {}
			for rightNeighbor in possiblepairs:
				sfreqs[startPair][endPair][leftNeighbor][rightNeighbor] = freqs[startPair][endPair][leftNeighbor][rightNeighbor] + freqs[endPair][startPair][leftNeighbor][rightNeighbor]

# startPair -> endPair -> leftNeighbor -> rightNeighbor
# add 1 to all freqs to ensure that there will be no prob 0 transition
nfreqs = {}
for startPair in possiblepairs:
	if not startPair in nfreqs:
		nfreqs[startPair] = {}
	for endPair in possiblepairs:
		if not endPair in nfreqs[startPair]:
			nfreqs[startPair][endPair] = {}
		for leftNeighbor in possiblepairs:
			if not leftNeighbor in nfreqs[startPair][endPair]:
				nfreqs[startPair][endPair][leftNeighbor] = {}
			for rightNeighbor in possiblepairs:
				nfreqs[startPair][endPair][leftNeighbor][rightNeighbor] = sfreqs[startPair][endPair][leftNeighbor][rightNeighbor] + 1
'''
# startPair -> leftNeighbor -> rightNeighbor -> endPair (frequency of triplets)
tfreqs = {}
for startPair in possiblepairs:
	if not startPair in tfreqs:
		tfreqs[startPair] = {}
	for leftNeighbor in possiblepairs:
		if not leftNeighbor in tfreqs[startPair]:
			tfreqs[startPair][leftNeighbor] = {}
		for rightNeighbor in possiblepairs:
			if not rightNeighbor in tfreqs[startPair][leftNeighbor]:
				tfreqs[startPair][leftNeighbor][rightNeighbor] = {}
			for endPair in possiblepairs:
				tfreqs[startPair][leftNeighbor][rightNeighbor][endPair] = nfreqs[startPair][endPair][leftNeighbor][rightNeighbor]
'''
totals = {}
for startPair in possiblepairs:
	if not startPair in totals:
		totals[startPair] = {}
	for leftNeighbor in possiblepairs:
		if not leftNeighbor in totals[startPair]:
			totals[startPair][leftNeighbor] = {}
		for rightNeighbor in possiblepairs:
			totals[startPair][leftNeighbor][rightNeighbor] = sum([nfreqs[startPair][endPair][leftNeighbor][rightNeighbor] for endPair in possiblepairs])

'''
for startPair in possiblepairs:
	if not startPair in totals:
		totals[startPair] = {}
	for endPair in possiblepairs:
		totals[startPair][endPair] = sum(nfreqs[startPair][endPair][leftNeighbor][rightNeighbor] for leftNeighbor in nfreqs[startPair][endPair] for rightNeighbor in nfreqs[startPair][endPair][leftNeighbor])
'''

# startPair -> endPair -> leftNeighbor -> rightNeighbor
probs = {}
for startPair in possiblepairs:
	if not startPair in probs:
		probs[startPair] = {}
	for endPair in possiblepairs:
		if not endPair in probs[startPair]:
			probs[startPair][endPair] = {}
		for leftNeighbor in possiblepairs:
			if not leftNeighbor in probs[startPair][endPair]:
				probs[startPair][endPair][leftNeighbor] = {}
			for rightNeighbor in possiblepairs:
				probs[startPair][endPair][leftNeighbor][rightNeighbor] = float(nfreqs[startPair][endPair][leftNeighbor][rightNeighbor])/totals[startPair][leftNeighbor][rightNeighbor]

# startPair -> endPair -> leftNeighbor -> rightNeighbor
neglogprobs = {}
for startPair in possiblepairs:
	if not startPair in neglogprobs:
		neglogprobs[startPair] = {}
	for endPair in possiblepairs:
		if not endPair in neglogprobs[startPair]:
			neglogprobs[startPair][endPair] = {}
		for leftNeighbor in possiblepairs:
			if not leftNeighbor in neglogprobs[startPair][endPair]:
				neglogprobs[startPair][endPair][leftNeighbor] = {}
			for rightNeighbor in possiblepairs:
				neglogprobs[startPair][endPair][leftNeighbor][rightNeighbor] = -math.log(probs[startPair][endPair][leftNeighbor][rightNeighbor])

# startPair -> endPair -> leftNeighbor
noRightNeighbor = {}
for startPair in possiblepairs:
	if not startPair in noRightNeighbor:
		noRightNeighbor[startPair] = {}
	for endPair in possiblepairs:
		if not endPair in noRightNeighbor[startPair]:
			noRightNeighbor[startPair][endPair] = {}
		for leftNeighbor in possiblepairs:
			noRightNeighbor[startPair][endPair][leftNeighbor] = sum([neglogprobs[startPair][endPair][leftNeighbor][rightNeighbor] for rightNeighbor in possiblepairs])/len(possiblepairs)

# startPair -> endPair -> rightNeighbor
noLeftNeighbor = {}
for startPair in possiblepairs:
	if not startPair in noLeftNeighbor:
		noLeftNeighbor[startPair] = {}
	for endPair in possiblepairs:
		if not endPair in noLeftNeighbor[startPair]:
			noLeftNeighbor[startPair][endPair] = {}
		for rightNeighbor in possiblepairs:
			noLeftNeighbor[startPair][endPair][rightNeighbor] = sum([neglogprobs[startPair][endPair][leftNeighbor][rightNeighbor] for leftNeighbor in possiblepairs])/len(possiblepairs)

noLeftOrRightNeighbor = {}
for startPair in possiblepairs:
	if not startPair in noLeftOrRightNeighbor:
		noLeftOrRightNeighbor[startPair] = {}
	for endPair in possiblepairs:
		noLeftOrRightNeighbor[startPair][endPair] = sum([neglogprobs[startPair][endPair][leftNeighbor][rightNeighbor] for leftNeighbor in possiblepairs for rightNeighbor in possiblepairs])/(len(possiblepairs)*len(possiblepairs))

heuristic = {}
for startPair in possiblepairs:
	if not startPair in heuristic:
		heuristic[startPair] = {}
	for endPair in possiblepairs:
		heuristic[startPair][endPair] = min([neglogprobs[startPair][endPair][leftNeighbor][rightNeighbor] for leftNeighbor in possiblepairs for rightNeighbor in possiblepairs])

#print heuristic

costs = neglogprobs

