# TODO : effet du biais

# 2 modes of speciation:
#	D1: within a clade, creating a new lineage
#	D2: origination from not sampled specimen/species of a new clade

# general model:
# 3 modes of extinction:
#	E1: extinction impacting only a single lineage, immediatly killing this lineage
#	E2: extinction impacting only a genus, immediatly killing a given percentage of lineages of this clade (pE2 = proportion of lineages extincted during an E2 event)
#	E3: extinction impacting only a geographic area, immediatly killin a given percentage of lineage within a concerned geography

# colonization:
#	C1: a given lineage has a given probability of moving into a different geographic area (but still the same lineage)

from numpy.random import poisson 
from numpy.random import binomial
from numpy.random import choice
from numpy.random import randint 

from numpy import nan
from numpy import nansum
from numpy import isnan

def get_geography(genus, lineage, nLocalities):
	currently_occuring = []
	for i in range(len(genus['lineage'])):
		if genus['lineage'][i]==lineage:
			currently_occuring.append(genus['geography'][i])
	currently_occuring = set(currently_occuring)
	left_localities = [ i for i in range(nLocalities) if i not in currently_occuring ]
	if len(left_localities) >0:
		res = choice(a=left_localities, size=1, replace=False)[0]
	else:
		res = -1
	return(res)

#simulation_name = 'S1'
# E1 probability for a single lineage/deme to be extinct
# E2 probability for a genus to experiment a mass extinction
# pE2 proportion of demes mass extincted by an E2 event
# E3 probability for a geography to experiment a mass extinction
# pE3 proportion of demes within the locality for being mass extincted by an E3 event

# D1
# D2

# C1

def simulation(rate, proportion, simulation_name, result_repertory, verbose, nb_localities=10, nb_generations = 600, nb_max_genus = 2000, nb_max_lineages_in_genus = 50):
	#### parameters used for trials
	nLocalities = nb_localities # number of geographic localities where species from different genus can be found
	nGenerations = nb_generations  # number of generations to simulate

	maxNGenus = nb_max_genus # maximum number of genus
	nMaxLineages_in_genus = nb_max_lineages_in_genus # a given genus can have a maximum number of species/lineages
	nMaxDemes_in_genus = nLocalities*nMaxLineages_in_genus

	E1, E2, E3, D1, D2, C1 = rate
	pE2, pE3 = proportion

	with open(f'{result_repertory}/{simulation_name}.par', 'w') as fileIn:
		fileIn.write(f'{simulation_name}\t{E1}\t{E2}\t{E3}\t{pE2}\t{pE3}\t{D1}\t{D2}\t{C1}\t{nLocalities}\t{nGenerations}\t{maxNGenus}\t{nMaxLineages_in_genus}\n')

	living_genus = [] # list of genus which 1) were born and 2) not dead yet

	record = {}
	for genus in range(maxNGenus):
		record[genus] = {} # record['genus']
		record[genus]['nLineages'] = 0
		record[genus]['lineage'] = [nan] * nMaxDemes_in_genus
		record[genus]['oldest'] = [nan] * nMaxDemes_in_genus
		record[genus]['youngest'] = [nan] * nMaxDemes_in_genus
		record[genus]['geography'] = [nan] * nMaxDemes_in_genus
		record[genus]['extinction_statu'] = [nan] * nMaxDemes_in_genus

	living_genus.append(0)

	record[0]['nLineages'] +=1
	record[0]['lineage'][0] = 0
	record[0]['oldest'][0] = nGenerations
	record[0]['youngest'][0] = nGenerations
	record[0]['geography'][0] = 0
	record[0]['extinction_statu'][0] = 0

	record_extincted = {}
	record_extincted['genus'] = []
	record_extincted['lineage'] = []
	record_extincted['geography'] = []
	record_extincted['youngest'] = []	
	record_extincted['oldest'] = []
	record_extincted['extinction_statu'] = []

	record_extincted['genus'] = [nan]*4000000
	record_extincted['lineage'] =[nan]*4000000 
	record_extincted['geography'] = [nan]*4000000
	record_extincted['youngest'] = 	[nan]*4000000
	record_extincted['oldest'] = [nan]*4000000
	record_extincted['extinction_statu'] = [nan]*4000000
	nExtincted_lineages = 0

	mass_extinctions_genus = []
	mass_extinctions_times = []
	mass_extinctions_kind = []
	mass_extinctions_effect = []
	mass_extinctions_deme = []

	last_created_genus = 0
	for time in range(nGenerations)[::-1]:
		# if time==250:
		# 	E1 = 0.5 # probability for a single lineage/deme to be extinct
		# 	E2 = 0.5 # probability for a genus to experiment a mass extinction
		# 	pE2 = 0.5 # proportion of demes mass extincted by an E2 event
		# 	E3 = 0.5 # probability for a geography to experiment a mass extinction
		# 	pE3 = 0.5 # proportion of demes within the locality for being mass extincted by an E3 event
		# if time==150:
		# 	E1 = 0.1 # probability for a single lineage/deme to be extinct
		# 	E2 = 0.1 # probability for a genus to experiment a mass extinction
		# 	pE2 = 0.1 # proportion of demes mass extincted by an E2 event
		# 	E3 = 0 # probability for a geography to experiment a mass extinction
		# 	pE3 = 0.1 # proportion of demes within the locality for being mass extincted by an E3 event
		if verbose:
			print('Generation {0}: {1} extincted lineages'.format(time, nExtincted_lineages))
		# origination of a new genus
		test_new_genus = binomial(n=1, p=D2, size=1)[0]
		if test_new_genus:
			last_created_genus += 1
			new_genus = last_created_genus
			living_genus.append(new_genus)
	#		record[new_genus] = {}
			record[new_genus]['nLineages'] = 1
	#		record[new_genus]['lineage'] = [nan] * nMaxDemes_in_genus
	#		record[new_genus]['oldest'] = [nan] * nMaxDemes_in_genus
	#		record[new_genus]['youngest'] = [nan] * nMaxDemes_in_genus
	#		record[new_genus]['geography'] = [nan] * nMaxDemes_in_genus
	#		record[new_genus]['extinction_statu'] = [nan] * nMaxDemes_in_genus
			
			record[new_genus]['lineage'][0] = 0
			record[new_genus]['oldest'][0] = time
			record[new_genus]['youngest'][0] = time
			new_geography = randint(low=0, high=nLocalities, size=1)[0]
			record[new_genus]['geography'][0] = new_geography
			record[new_genus]['extinction_statu'][0] = 0
			
		# loop over different genus for the origination of new lineage
	#	for genus in record.keys():
		for genus in living_genus:
			# only treats genus that are not fully extincted
			if record[genus]['lineage'].count(nan)!=len(record[genus]['lineage']):
				# origination of a new lineage
				n_non_extincted_demes = len(record[genus]['extinction_statu']) - nansum(record[genus]['extinction_statu']) - record[genus]['extinction_statu'].count(nan)
				n_new_lineages = sum(binomial(n=1, p=D1, size=int(n_non_extincted_demes)))
				if n_new_lineages>0:
					for lineage in range(n_new_lineages):
						# if the genus is not "full"
						if record[genus]['lineage'].count(nan)>0:
							newPos = record[genus]['lineage'].index(nan)
							record[genus]['nLineages'] += 1
							
							record[genus]['lineage'][newPos] = record[genus]['nLineages']-1
							record[genus]['oldest'][newPos] = time
							record[genus]['youngest'][newPos] = time
							
							new_geography = randint(low=0, high=nLocalities, size=1)[0]
							record[genus]['geography'][newPos] = new_geography # the new lineage appears in a new geography
							record[genus]['extinction_statu'][newPos] = 0

				# colonization of a new geographic area
				n_non_extincted_demes = len(record[genus]['extinction_statu']) - nansum(record[genus]['extinction_statu']) - record[genus]['extinction_statu'].count(nan)
				n_colonizing_lineages = sum(binomial(n=1, p=C1, size=int(n_non_extincted_demes)))
				
				if n_colonizing_lineages>0:
					copied_deme = [ i for i in range(len(record[genus]['lineage'])) if isnan(record[genus]['lineage'][i])==False ]
					copied_deme = choice(a=copied_deme, size=n_colonizing_lineages, replace=True) # a deme can colonize multiple new localities at the same time
					for colonizer in copied_deme:
						if record[genus]['lineage'].count(nan)>0:
							new_geography = get_geography(record[genus], record[genus]['lineage'][colonizer], nLocalities)
							if new_geography != -1:
								newPos = record[genus]['lineage'].index(nan)
								record[genus]['lineage'][newPos] = record[genus]['lineage'][colonizer]
								record[genus]['oldest'][newPos] = time
								record[genus]['youngest'][newPos] = time
								record[genus]['geography'][newPos] = new_geography
								record[genus]['extinction_statu'][newPos] = 0
				
				# extinction of a given deme (lineage in a locality)
				current_demes = [ i for i in range(len(record[genus]['lineage'])) if isnan(record[genus]['lineage'][i])==False ]
				n_non_extincted_demes = len(current_demes)
				n_lineages_to_kill = sum(binomial(n=1, p=E1, size=int(n_non_extincted_demes)))
				
				if n_lineages_to_kill>n_non_extincted_demes:
					n_lineages_to_kill=n_non_extincted_demes
				
				if n_lineages_to_kill>0:
					killed_deme = choice(a=current_demes, size=n_lineages_to_kill, replace=False) # a deme cannot be extincted more that one time
					for death in killed_deme:
							record[genus]['extinction_statu'][death] = 1
					
	#						record_extincted['genus'].append(genus)
	#						record_extincted['lineage'].append( record[genus]['lineage'][death] )
	#						record_extincted['geography'].append( record[genus]['geography'][death] )
	#						record_extincted['youngest'].append( time )
	#						record_extincted['oldest'].append( record[genus]['oldest'][death] )
	#						record_extincted['extinction_statu'].append( record[genus]['extinction_statu'][death] )
							
							record_extincted['genus'][nExtincted_lineages] = genus
							record_extincted['lineage'][nExtincted_lineages] = record[genus]['lineage'][death] 
							record_extincted['geography'][nExtincted_lineages] = record[genus]['geography'][death] 
							record_extincted['youngest'][nExtincted_lineages] = time 
							record_extincted['oldest'][nExtincted_lineages] = record[genus]['oldest'][death] 
							record_extincted['extinction_statu'][nExtincted_lineages] = record[genus]['extinction_statu'][death] 
							nExtincted_lineages += 1
							
							record[genus]['lineage'][death] = nan
							record[genus]['oldest'][death] = nan 
							record[genus]['youngest'][death] = nan 
							record[genus]['geography'][death] = nan 
							record[genus]['extinction_statu'][death] = nan
		
		# loop over different genus for the extinction of a genus
		for genus in living_genus:
			if record[genus]['lineage'].count(nan)!=len(record[genus]['lineage']):
				# number of non extincted lineages
				non_extincted_demes = [ i for i in range(len(record[genus]['lineage'])) if isnan(record[genus]['lineage'][i])==False ]
				n_non_extincted_demes = len(non_extincted_demes)
				
				# test extinction
				test_extinction_genus = binomial(n=1, p=E2, size=1) ## test extinction at the start to limitate calcul
				if test_extinction_genus:
					nDemes_to_extinct = poisson(lam=pE2*n_non_extincted_demes, size=1)[0]
					if nDemes_to_extinct>n_non_extincted_demes:
						nDemes_to_extinct=n_non_extincted_demes

					if nDemes_to_extinct>0:
						mass_extinctions_genus.append(genus)
						mass_extinctions_deme.append(-1)
						mass_extinctions_times.append(time)
						mass_extinctions_kind.append('E2')
						mass_extinctions_effect.append(nDemes_to_extinct)
						# list of demes to extinct
						non_extincted_demes = [ i for i in range(len(record[genus]['extinction_statu'])) if record[genus]['extinction_statu'][i]==0 ]
						demes_to_extinct = choice( a=non_extincted_demes, size=nDemes_to_extinct, replace=False)
						
						for i in demes_to_extinct:
							record[genus]['extinction_statu'][i] = 1
							record[genus]['youngest'][i] = time

		#					record_extincted['genus'].append(genus)
		#					record_extincted['lineage'].append( record[genus]['lineage'][i] )
		#					record_extincted['geography'].append( record[genus]['geography'][i] )
		#					record_extincted['youngest'].append( record[genus]['youngest'][i] )
		#					record_extincted['oldest'].append( record[genus]['oldest'][i] )
		#					record_extincted['extinction_statu'].append( record[genus]['extinction_statu'][i] )

							record_extincted['genus'][nExtincted_lineages] = genus
							record_extincted['lineage'][nExtincted_lineages] = record[genus]['lineage'][i] 
							record_extincted['geography'][nExtincted_lineages] = record[genus]['geography'][i] 
							record_extincted['youngest'][nExtincted_lineages] = record[genus]['youngest'][i]
							record_extincted['oldest'][nExtincted_lineages] = record[genus]['oldest'][i] 
							record_extincted['extinction_statu'][nExtincted_lineages] = record[genus]['extinction_statu'][i] 
							nExtincted_lineages += 1

							record[genus]['lineage'][i] = nan
							record[genus]['oldest'][i] = nan
							record[genus]['youngest'][i] = nan
							record[genus]['geography'][i] = nan
							record[genus]['extinction_statu'][i] = nan
		
		test_extinction_deme = binomial(n=1, p=E3, size=1)
		if test_extinction_deme:
			localities_record = {}
			
			# regroup lineages by localities
			for genus in record.keys():
				for i in range(len(record[genus]['geography'])):
					current_localities = record[genus]['geography'][i]
					if isnan(current_localities) == False:
						if current_localities not in localities_record.keys():
							localities_record[current_localities] = {}
							localities_record[current_localities]['lineage'] = [(genus, record[genus]['lineage'][i])] + ([None] * (maxNGenus*nMaxDemes_in_genus))
							localities_record[current_localities]['index'] = [i] + ([None] * (maxNGenus*nMaxDemes_in_genus))
						else:
							new_index = localities_record[current_localities]['lineage'].index(None)
							localities_record[current_localities]['lineage'][new_index] = (genus, record[genus]['lineage'][i])
							localities_record[current_localities]['index'][new_index] = i

			# define a deme where extinction will take place
			if len(localities_record.keys()) != 0:
				demes_to_extinct = choice(a=list(localities_record.keys()), size=1, replace=False)[0]
				non_extincted_lineages = [ i for i in range(len(localities_record[demes_to_extinct]['lineage'])) if localities_record[demes_to_extinct]['lineage'][i] is not None ]
				n_non_extincted_lineages = len(non_extincted_lineages)
				nLineages_to_extinct = poisson(lam=pE3*n_non_extincted_lineages, size=1)[0]
				if nLineages_to_extinct > n_non_extincted_lineages:
					nLineages_to_extinct = n_non_extincted_lineages
				if nLineages_to_extinct > 0:
					mass_extinctions_genus.append(-1) # not adaptate to E3
					mass_extinctions_deme.append(demes_to_extinct)
					mass_extinctions_times.append(time)
					mass_extinctions_kind.append('E3')
					mass_extinctions_effect.append(nLineages_to_extinct)
					# list of lineages to extinct
					lineages_to_extinct = choice(a=non_extincted_lineages, size=nLineages_to_extinct, replace=False)
					for i in lineages_to_extinct:
						index = localities_record[demes_to_extinct]['index'][i]
						genus = localities_record[demes_to_extinct]['lineage'][i][0]
						record_extincted['genus'][nExtincted_lineages] = genus
						record_extincted['lineage'][nExtincted_lineages] = localities_record[demes_to_extinct]['lineage'][i][1]
						record_extincted['geography'][nExtincted_lineages] = demes_to_extinct
						record_extincted['youngest'][nExtincted_lineages] = time
						record_extincted['oldest'][nExtincted_lineages] = record[genus]['oldest'][index]
						record_extincted['extinction_statu'][nExtincted_lineages] = 1
						nExtincted_lineages += 1

						record[genus]['lineage'][index] = nan
						record[genus]['oldest'][index] = nan
						record[genus]['youngest'][index] = nan
						record[genus]['geography'][index] = nan
						record[genus]['extinction_statu'][index] = nan
						
			

	# output file
	## fossils
	outfile = open(f'{result_repertory}/simulated_fossils_{simulation_name}.txt', 'w')
	outfile.write('simulation\tspecimen\tgenus\tspecies\tgeography\tmax_ma\tmin_ma\textinction_statu\n')

	for i in range(nExtincted_lineages):
		outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t1\n'.format(simulation_name, i, record_extincted['genus'][i], record_extincted['lineage'][i], record_extincted['geography'][i], record_extincted['oldest'][i],  record_extincted['youngest'][i]))

	cnt = i
	## current species
	for genus in record:
		if len(record[genus]['lineage'])!=record[genus]['lineage'].count(nan):
			for i in range(len(record[genus]['lineage'])):
				if isnan(record[genus]['lineage'][i])==False:
					cnt += 1
					outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t0\n'.format(simulation_name, cnt, genus, record[genus]['lineage'][i], record[genus]['geography'][i], record[genus]['oldest'][i], time))
				

	outfile.close()

	## mass extinction events
	with open(f'{result_repertory}/simulated_mass_extinctions_{simulation_name}.txt', 'w') as outfile:
		outfile.write('simulation\ttime\taffected genus\taffected deme\tkind of extinction\tnumber of extincted demes\n')
		for i in range(len(mass_extinctions_genus)):
			outfile.write(f'{simulation_name}\t{mass_extinctions_times[i]}\t{mass_extinctions_genus[i]}\t{mass_extinctions_deme[i]}\t{mass_extinctions_kind[i]}\t{mass_extinctions_effect[i]}\n')

