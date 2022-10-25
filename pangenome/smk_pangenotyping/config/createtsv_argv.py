import os 
import sys
import glob
import pandas as pd 
import itertools

def main(p,bam_cram, tab_f):
	# define the list
	sample=[]
	bam_cram=bam_cram
	for f in os.listdir(p):
		if os.path.isfile(os.path.join(p, f)):
			if f.endswith('.'+bam_cram):
				sample.append((os.path.splitext(f)[0], p))
	
	df = pd.DataFrame(sample, columns= ['sample_id','sample_path'])
	df.to_csv(tab_f, sep="\t", index= False)

if __name__ == '__main__':

	p=sys.argv[1]
	bam_cram=sys.argv[2]
	tab_f=sys.argv[3]
	# os.makedirs(outf,exist_ok = True)
	main(p,bam_cram, tab_f)