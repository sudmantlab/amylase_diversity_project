import pandas as pd
import os
import sys
from glob import glob
from datetime import datetime

def organizer_samples():

	'''
	Organize input-output for better usage of Snakemake rules - Illumina
	'''

	#samples

	if not os.path.exists(config['samples']):

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Error] Input sample table does not exist')
		sys.exit(1)
	
	base=os.path.dirname(os.path.abspath(os.path.dirname(config['samples'])))

	try:

		df=(pd.read_table(config['samples'], dtype={"sample_id": str, "sample_path":str})
		.set_index("sample_id", drop=False)
		.sort_index()
		)

	except:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Error] Incorrect format of the input sample table')
		sys.exit(1)


	for index, row in df.iterrows():

		x=row['sample_path']
		
		if not os.path.exists(os.path.abspath(x)):

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + ']' + '[Error] Directory ' + os.path.abspath(x) + '  does not exist')
			sys.exit(1)

		y=glob(os.path.abspath(x + "/" + index) + "*am*") #retrive file with id
		resources=os.path.abspath(base + '/resources/aligned/' + index)

		if not os.path.exists(resources):

			os.makedirs(resources)

		for fq in y: #create symbolik links

			if not os.path.exists(os.path.abspath(resources + '/' + os.path.basename(fq))):

				target=os.path.abspath(resources + '/' + os.path.basename(fq))
				print(target)
				os.symlink(os.path.abspath(fq), target)

	return df

def organizer_reference():

	'''
	Create symlink to reference file - 
	'''

	if not os.path.exists(config['reference']['localpath']):

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Error] Input reference directory does not exist')
		sys.exit(1)

	base=os.path.dirname(os.path.abspath(os.path.dirname(config['samples'])))
	resources=os.path.abspath(base + '/resources/reference')
	os.makedirs(resources,exist_ok=True)
	target=os.path.basename(config['reference']['localpath'])

	if not os.path.exists(os.path.abspath(resources + '/' + target)):
	
		os.symlink(config['reference']['localpath'], os.path.abspath(resources + '/' + target))



def organizer_odgi_input():

	'''
	Create symlink to odgi file with region - 
	'''

	if not os.path.exists(config['odgi']['chop_input_local_path']):

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + ']' + '[Error] Input reference directory does not exist')
		sys.exit(1)

	base=os.path.dirname(os.path.abspath(os.path.dirname(config['samples'])))
	resources=os.path.abspath(base + '/resources/odgi')
	os.makedirs(resources,exist_ok=True)
	target=os.path.basename(config['odgi']['chop_input_local_path'])

	if not os.path.exists(os.path.abspath(resources + '/' + target)):
	
		os.symlink(config['odgi']['chop_input_local_path'], os.path.abspath(resources + '/' + target))


