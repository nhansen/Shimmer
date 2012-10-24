# calculates sensitivity and specificity from counts of TPs, FPs, FNs

class mutation:
        fields = []
        def __init__(self):
                self.properties={}
                for field in mutation.fields:
                        self.properties[field]=[]

def getFields(vsfile):
        import os,sys,re

	file_by_chr = {} #dictionary of chromosomes

	v = open(vsfile)
	fields = v.readline().rstrip('\n').split("\t")
	mutation.fields=fields
	chr_ind = fields.index("Chr")
	var_ind=fields.index('var_allele')
	ref_ind=fields.index('ref_allele')

	for line in v:
		li = line.split("\t")
		if li[var_ind].upper()!=li[ref_ind].upper():
			new_mut = mutation()

			#get the field properties of mutation from the vsfile fields
			for i,field in enumerate(li):
				new_mut.properties[fields[i]] = field

			#mutation is added to the chromosome to which it belongs
			if li[chr_ind] in file_by_chr:
				file_by_chr[li[chr_ind]].append(new_mut)
			else:
				file_by_chr[li[chr_ind]]=[new_mut]

        return file_by_chr

# parse out the mutations -- only care about SNPs, not indels
def parse_mutations(cosmic):
        import os,sys,re
        mutations = {} #dictionary {chr: {pos:mut}}

        f = open(cosmic)
        for line in f:
                li = line.strip('\n').split('\t')
                chr,start,end,type,strand=li[0],int(li[1]),int(li[2]),li[3],li[4]

                # convert all to forward strand
                mut=''
		if re.match('\w>\w',type): #if SNP
			if strand=='-':
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
                	else:
                        	mut=type

			# add to dict of mutations
                	if chr in mutations:
                        	mutations[chr][start]=mut
                	else:
                        	mutations[chr]={start:mut}

        return mutations

def simulation_performance(cosmic,vsfile):
	import os,re,sys
	truemuts = parse_mutations(cosmic)
	vsmuts = getFields(vsfile)
	
	totalmuts=0
	for chr in truemuts.iterkeys():
		totalmuts+=len(truemuts[chr])

	fpcount=0
	tpcount=0
	fncount=0
	
	for chr in vsmuts.iterkeys():
		for mut in vsmuts[chr]:
			pos = int(mut.properties['LeftFlank'])+1
			if (chr in truemuts) and (pos in truemuts[chr]):
				assert mut.properties['var_allele']==truemuts[chr][pos][2],'detected variant allele incorrect'
				tpcount+=1
			else:
				fpcount+=1
	fncount=totalmuts-tpcount

	print tpcount
	print fpcount
	print fncount
	
	sensitivity = float(tpcount)/float(tpcount+fncount)
	specificity = float(tpcount)/float(tpcount+fpcount)
	
	print 'sensitivity = '+`sensitivity`
	print 'specificity = '+`specificity`

if __name__ == "__main__":
        import sys
        simulation_performance(sys.argv[1],sys.argv[2])
