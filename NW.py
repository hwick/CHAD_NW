#def tyler_NW(seq, primer):
#for p in primers:
#print "\tLooking at primer "+primer
seq = 'GGTTGATGACAGGGCAGCTACACAGCCCATTTCGAGAGAGAGAGAGATACCTTTTATAT'
seq = 'GGTTGATGACAAAGCAGCTACACAGCCCATTTCGAGAGAGAGAGAGATACCTTTTATAT'
#seq = 'GGTTGATGACAGGGCAGCTACACAGAGTGTGAGCCCATTTCGAGAGAGAGAGAGATACCTTTTATAT' 
primer = "GGACACCTACACAGCCCAT"	
#primer = "GGACACCTACCCCACAGCCCAT"
seq = "CCTCTTGCATGCATGTTCTCCTTGGTCCTG" #actual sequence for below primer from data
primer = "TCGTCCTGATCCCTCTTCTC" #HLAB-F 3'-5' reverse directeion
seq = "CCCAGGTCTCGGTCAGGCCAGGCCTCCCG" #actual sequence for below primer from data
seq = "CCCAGGTCTCGGTCAGACCAGGCCTCCCG" # sequence with incorrect 3'-most base
primer = "GCCCAGGTCTCGGTCAGG" #HLAB-F 5'-3' forward direction

#seq = "CCCAGGTCTCGGTCCTGATCCCTCTTCTC" #actual sequence for below primer from data
#primer = 'GCCCAGGTCTCGGTCAGG' #HLAB-F 5'-3' forward direction
print seq
print primer
#primer = p
matchscore = 1
mismatch = -1
gap = -2
S = [] # sequence "matrix" (list of lists)
T = [] # traceback "matrix" (list of lists)
#build score matrix and traceback matrix (actually list (length of sam sequence) of lists (length of primer))
for x in range(len(seq)+1):
	S.append([0]*(len(primer)+1))
	T.append([0]*(len(primer)+1))
#Walk through the matrices, starting at position 1,1
#and proceding row by row, left to right/top to bottom, assigning a score.
for y in range(1,len(primer)+1):
	for x in range(1,len(seq)+1):
		if seq[x-1] == primer[y-1]:#See if it's a match
		#Set this cell equal to the best possible score out of three options:
		#Either it's a gap coming from the cell above, a gap from the left,
		#or a match from the upper-left.
			S[x][y] = max(S[x-1][y]+gap,S[x][y-1]+gap,S[x-1][y-1] + matchscore)
			if S[x][y] == S[x-1][y-1]+matchscore:#if the match is the best score, mark it
				T[x][y] = 1
			elif S[x][y] == S[x-1][y]+gap:
				T[x][y] = 2
			else:
				T[x][y] = 3
		else: #if it's not a match, check whether the best score would be from mismatch or the two gap options
			S[x][y] = max(S[x-1][y]+gap,S[x][y-1]+gap,S[x-1][y-1]+mismatch)
			if S[x][y] == S[x-1][y-1]+mismatch:
				T[x][y] = 1
			elif S[x][y] == S[x-1][y]+gap:
				T[x][y] = 2
			else:
				T[x][y] = 3
#Find the maximum score in the last row and last column and set the rest of the remaining matrix in the last row/column to that number, this lets you have trailing gaps
#(instead of just starting in the bottom right hand corner like in a "usual" NW aligner)
lastColumn = []
for x in range(len(seq)+1):
	lastColumn.append(S[x][-1])
for x in range(1,len(seq)):
	if S[x][-1] == max(lastColumn):
		S[x+1][-1] = S[x][-1]
		T[x+1][-1] = 2
for y in range(1,len(primer)):
	if S[-1][y] == max(S[-1][:]):
		S[-1][y+1] = S[-1][y]
		T[-1][y+1] = 3
score = S[-1][-1]
max_score = len(primer)
distance_from_max = max_score-score
print str(S[-1][-1]) + " out of "+str(len(primer))+" possible.\t\t\t\t\t\t\t\t"+str(len(primer)-S[-1][-1])+"<-distance from maximum score"
str(S[-1][-1]/len(primer)*100)+"% (this is not a real percent because the score could be negative if it's really bad)"
	
# This will print out the two "matrices" #was commented out in actual program
for x in range(len(seq)+1): #was commented out in actual program
    print S[x] #was commented out in actual program
print "Traceback:" #was commented out in actual program
for x in range(len(seq)+1): #was commented out in actual program
	    print T[x] #was commented out in actual program

# Build the alginment to print it out
S = [] #aligned sequence
P = [] #aligned primer
x = len(seq)
y = len(primer)
while x > 0 and y > 0:
	if T[x][y] == 1:
		S.append(seq[x-1])
		P.append(primer[y-1])
		x -= 1
		y -= 1
	elif T[x][y] == 2:
		S.append(seq[x-1])
		P.append("-")
		x -= 1
	elif T[x][y] == 3:
		S.append("-")
		P.append(primer[y-1])
		y -= 1
while x > 0:
	S.append(seq[x-1])
	P.append("-")
	x -= 1
while y > 0:
	S.append("-")
	P.append(primer[y-1])
	y -= 1
aligned_s = ""
for i in range(len(S)): #was commented out in actual program
	aligned_s += S[-i-1] #was commented out in actual program
print aligned_s+" <- aligned_s from commented out section"
#print "\t\t"+sS #was commented out in actual program
aligned_p = ""
for i in range(len(P)): #was commented out in actual program
	aligned_p += P[-i-1] #was commented out in actual program
print aligned_p+" <- aligned_p from commented out section"
stringS = "".join(S)
aligned_s = stringS[::-1]
stringP = "".join(P)
aligned_p = stringP[::-1]
print aligned_s+" <- aligned_s that would get returned"
print aligned_p+" <- aligned_s that would get returned"
print stringS+" <- stringS doesn't get returned what is this"
print stringP+" <- stringP doesn't get returned what is this"
#print "\t\t"+sP #was commented out in actual program
print "\t\t\t\t\t\t\t\t"+str(score)+"\t"+str(max_score)+"\t"+str(distance_from_max)+"<- distance from maximum score\t" #was commented out in actual program
print len(primer)-S[-1][-1] #was commented out in actual program
print len(primer) #was commented out in actual program
print S[-1][-1] #was commented out in actual program
print
print
#	return(aligned_s, aligned_p, score)
