#takes as input the original bam file and a file of mutations in the cosmic format
#outputs a bam file that is the same as the original but with the mutations in the cosmic file

def parse_mutations(cosmic):
	import os,sys,re
	mutations = {} #dictionary {chr: {pos:mut}}
	
	f = open(cosmic)
	for line in f:
		li = line.strip('\n').split('\t')
		chr,start,end,type,strand=li[0],int(li[1]),int(li[2]),li[3],li[4]

		# convert all to forward strand
		mut=''
		if strand=='-':
			#convert each base to its complement
			for i,char in enumerate(type):
				if char=='C':
					mut+='G'
				elif char=='G':
					mut+='C'
				elif char=='T':
					mut+='A'
				elif char=='A':
					mut+='T'
				else:
					mut+=type[i]
			if re.match('ins\w+',type):
                        	mut=mut[0:3]+mut[len(mut):2:-1] #reverse the order of insertion sequence       
		else:	
			mut=type	
	
		if chr in mutations:
			mutations[chr][start]=mut
		else:
			mutations[chr]={start:mut}
	
	return mutations
	
def create_mut(mutations,ar,percent):
	import random,pysam,sys,re
	tempq=ar.qual
	nr=ar
	seq = [] #list of tuples (pos, nuc)
	cigar=[] #0=M,1=I,2=D,3=N,4=S,5=H,6=P,7==,8=X
	mutseq = '' #mutated sequence
	hardclips=[]

	if ar.cigar is None:
		return ar
	
	#expand cigar string out into a list (e.g. if alignedread.cigar()=[(1,1),(2,0)] then cigar=[1,0,0])
	for tup in ar.cigar:
		cigar=cigar+[tup[0]]*tup[1]
	
	#remove hardclips at the start of a cigar string so that the sequence(which doesn't contain the hardclipped region) and cigar string match up
	while cigar[0]==5:
		cigar.remove(5)
		hardclips.append(5)

	refcount=0
	pos=ar.pos+1
	for i,nuc in enumerate(ar.seq):
		if cigar[i]!=1:
			seq.append((pos+refcount,nuc))
			refcount+=1
		else:
			seq.append((-100,nuc)) #if mutation is insertion, position is irrelevant b/c pos refers to reference

	loopseq = enumerate(seq)
	for i,nuc in loopseq:
		if nuc[0] in mutations:
			if random.random()<=percent: #mutate with a given probability
				mut=mutations[nuc[0]]
				
				if re.match('\w+>\w+',mut): #if it's a SNP
					mutseq+=mut[2] #mutate the base
				elif re.match('ins\w+',mut): #if it's an insertion
					insertion=re.sub('ins','',mut)
					mutseq+=nuc[1]+insertion #insert after current position
					for j in range(0,len(insertion)):
						cigar.insert(i+j+1,1) #insert '1' for insertion into cigar string
						tempq=tempq[0:i+1]+tempq[i]+tempq[i+1:len(tempq)] #make insertion's quality same as quality of base before insertion
				elif re.match('del\w+',mut): #if it's a deletion
					deletion=re.sub('del','',mut)
					for j in range(0,len(deletion)):
						if i+j<len(cigar):
							cigar[i+j]=2 #change the cigar list to contain deletion
						else:
							break
					tempq=tempq[0:i]+tempq[i+len(deletion):len(tempq)]
					#skip over the next few that are also deletions
					for j in range(0,len(deletion)-1):	
						#skip over next positions in sequence
						try:
							loopseq.next()
							continue 
						except StopIteration:
							break
				else:
					mutseq+=nuc[1]
			else:
				mutseq+=nuc[1]
		else:
			mutseq+=nuc[1]

	cigar=hardclips+cigar #add hardclips back into cigar list
	
	#convert cigar list back into pysam cigar format
	cig=[]
	count=1
	for i in range(0, len(cigar)-1):
		if cigar[i]==cigar[i+1]:
			count+=1
		else:
			cig.append((cigar[i],count))
			count=1
	cig.append((cigar[len(cigar)-1],count))
	cig=tuple(cig)
	
	nr.cigar=cig
	nr.seq=mutseq
	nr.qual=tempq
	return nr

#origbam: original file to mutate from
#cosmic: mutations
#frac: percent to mutate
#outbam: path to name new bam file
def mutate_sample(origbam,cosmic,frac,outbam):
	import pysam, os, sys, re

	percent = float(frac)/2
	mutations = parse_mutations(cosmic) #dictionary {chr#: {..., pos:'B>B',...}, ...}
	orig = pysam.Samfile(origbam,'rb')
	tumor = pysam.Samfile(outbam+'.sim'+frac+'.bam','wb',template=orig)

	# set of alignedReads that need to be checked to mutate
	toCheck = set([])
	for chr in mutations.iterkeys():
		for pos in mutations[chr].iterkeys():
			for ar in orig.fetch(chr,pos-1,pos):
				toCheck.add(ar.qname)

	#fetch all the reads
	for ar in orig.fetch():
		if ar.qname in toCheck: #if read spans a mutation position, mutate it with input probability
			ref=pysam.Samfile.getrname(orig,ar.tid)
			nr=create_mut(mutations[ref],ar,percent)
			tumor.write(nr)
		else:
			tumor.write(ar)
	
	tumor.close()
	orig.close()

if __name__ == "__main__":
	import sys
	mutate_sample(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
