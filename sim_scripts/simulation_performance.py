# calculates sensitivity and specificity from counts of TPs, FPs, FNs

class mutation:
        fields = []
        def __init__(self):
                self.properties={}
                for field in mutation.fields:
                        self.properties[field]=[]

def getFields(vsfile, prog, cutoff):
        import os,sys,re

	file_by_chr = {} #dictionary of chromosomes

	v = open(vsfile)
	fields = v.readline().rstrip('\n').split("\t")
	mutation.fields=fields
	chr_ind = fields.index("Chr")
	var_ind=fields.index('var_allele')
	ref_ind=fields.index('ref_allele')
	pos_ind=fields.index('LeftFlank')

	last_chr = ""
	last_pos = 0

	for line in v:
		li = line.split("\t")
		if li[var_ind].upper()!=li[ref_ind].upper():
			if prog == "varscan":
				p_index=fields.index('p_value')
				pvalue = float(li[p_index])
			if prog == "shimmer":
				q_index=fields.index('q_value')
				qvalue = float(li[q_index])
			if prog == "somaticsniper":
				ss_index=fields.index('Somatic_Score')
				ssvalue = float(li[ss_index])
			if prog == "jsnvmix":
				ps_index=fields.index('p_somatic_genotype')
				psvalue = float(li[ps_index])
			cutoff = float(cutoff)
			if prog == "varscan" and pvalue > cutoff:
				continue

			if prog == "shimmer" and qvalue > cutoff:
				continue

			if prog == "somaticsniper" and ssvalue < cutoff:
				continue

			if prog == "jsnvmix" and psvalue < cutoff:
				continue

			new_mut = mutation()

			#get the field properties of mutation from the vsfile fields
			for i,field in enumerate(li):
				new_mut.properties[fields[i]] = field

			#mutation is added to the chromosome to which it belongs
			# check to be sure this isn't a dup:
			if li[chr_ind] == last_chr and li[pos_ind]==last_pos:
				continue

			if li[chr_ind] in file_by_chr:
				file_by_chr[li[chr_ind]].append(new_mut)
			else:
				file_by_chr[li[chr_ind]]=[new_mut]

			last_chr = li[chr_ind]
			last_pos = li[pos_ind]
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

def simulation_performance(cosmic,subdir,prog,cutoff,fh):
	import os,re,sys,glob,subprocess
	truemuts = parse_mutations(cosmic)

	vsfile = '/cluster/ifs/projects/seqSim/shimmerpaper/newsims/'+prog+'/'+subdir+'/'
	if prog == "shimmer":
		vsfile = vsfile + 'somatic_diffs.ANN.vs'
	if prog == "varscan":
		possfiles = glob.glob(vsfile+'*'+'.ann.vs')
		if len(possfiles) > 1:
			print 'More than one file in varscan directory!'
		vsfile = possfiles[0]
	if prog == "somaticsniper":
		possfiles = glob.glob(vsfile+'*'+'.finalfilter.dbSNP.filtered')
		if len(possfiles) > 1:
			print 'More than one file in somaticsniper directory!'
		vsfile = possfiles[0]
	if prog == "deepsnv":
		vsfile = vsfile + 'deepSNV.final.out'
		execstring = '/cluster/ifs/projects/seqSim/shimmerpaper/newsimplots/compare_deepSNV_to_truth.pl '+vsfile+' '+cosmic+' '+str(cutoff)+'\n'
		output = os.popen(execstring).read()
		#print output
		fh.write(output)
	if prog == "jsnvmix":
		possfiles = glob.glob(vsfile+'*'+'.dbSNP.filtered')
		if len(possfiles) > 1:
			print 'More than one file in jsnvmix directory!'
		vsfile = possfiles[0]

	if prog != "deepsnv":
		vsmuts = getFields(vsfile, prog, cutoff)
	else:
		vsmuts = {}
	
	totalmuts=0
	for chr in truemuts.iterkeys():
		totalmuts+=len(truemuts[chr])

	#print 'Total muts '+str(totalmuts)+' detected\n'
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
				print prog + '\t' + chr + '\t' + str(pos) + '\t' + mut.properties['var_allele']

	fncount=totalmuts-tpcount

	if prog != "deepsnv":
		fh.write(str(cutoff) + ' ' + str(tpcount) + ' ' + str(fpcount) + ' ' + str(fncount) + '\n')

	sensitivity = 0
	specificity = 0
	if tpcount+fncount:	
		sensitivity = float(tpcount)/float(tpcount+fncount)
	if tpcount+fpcount:	
		specificity = float(tpcount)/float(tpcount+fpcount)

	print prog + ':' + str(cutoff) + ':sens='+`sensitivity`+',spec='+`specificity`

if __name__ == "__main__":
        import os,sys

	outdir = sys.argv[3]

	# shimmer:
	if not os.path.exists(outdir):
		os.makedirs(outdir)	

	fh = open(outdir + '/shimmer_sens_spec.txt', 'w')
	shim_cuts = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5]
	for cutoff in shim_cuts:
        	simulation_performance(sys.argv[1],sys.argv[2],'shimmer',cutoff,fh)
	fh.closed
	# varscan:
	fh = open(outdir + '/varscan_sens_spec.txt', 'w')
	vs_cuts = [0.01, 0.005, 0.002, 0.0005, 0.0002, 0.00005, 0.00002]
	for cutoff in vs_cuts:
        	simulation_performance(sys.argv[1],sys.argv[2],'varscan',cutoff,fh)
	fh.closed
	# deepsnv:
	fh = open(outdir + '/deepsnv_sens_spec.txt', 'w')
	deepcuts = [0.1, 0.2, 0.9999]
	for cutoff in deepcuts:
		simulation_performance(sys.argv[1],sys.argv[2],'deepsnv',cutoff,fh)
	fh.closed
	# sniper:
	fh = open(outdir + '/sniper_sens_spec.txt', 'w')
	snipercuts = [40, 45, 50, 60, 70, 90, 110, 130]
	for cutoff in snipercuts:
        	simulation_performance(sys.argv[1],sys.argv[2],'somaticsniper',cutoff,fh)
	fh.closed
	# jsnvmix:
	fh = open(outdir + '/jsnvmix_sens_spec.txt', 'w')
	jsnvmixcuts = [1.0]
	for cutoff in jsnvmixcuts:
        	simulation_performance(sys.argv[1],sys.argv[2],'jsnvmix',cutoff,fh)
	fh.closed
