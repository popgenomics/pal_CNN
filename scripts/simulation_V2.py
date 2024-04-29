# TODO : effet du biais

# 3 modes of speciation:
#	D1: within a clade, creating a new lineage
#	D2: origination from not sampled specimen/species of a new clade
#	D3: within a clade, creating a new lineage but replace the original lineage i the deme

# general model:
# 3 modes of extinction:
#	E1: extinction impacting only a single lineage, immediatly killing this lineage
#	E2: extinction impacting only a genus, immediatly killing a given percentage of lineages of this clade (pE2 = proportion of lineages extincted during an E2 event)
#	E3: extinction impacting only a geographic area, immediatly killin a given percentage of lineage within a concerned geography

# colonization:
#	C1: a given lineage has a given probability of moving into a different geographic area (but still the same lineage)

import numpy as np
from math import exp
from random import choices

import pandas as pd
from numpy.random import poisson
from numpy.random import binomial
from numpy.random import choice
from numpy.random import randint 

from numpy import nan
from numpy import nansum
from numpy import isnan

def get_geography_weighted_colonization(genus, lineage, nLocalities, s_deme, sigma_k, lineage_index):
	"""
	Return one geography not currently occupied by the lineage

	Parameters
	----------
	genus : dict
		The dict use in the simulation to represent the genus
	lineage : str
		Current lineage which we want estimate his new localities
	nLocalities : int
		Number of localities in the simulation
	s_deme : dict
		All the values of s_deme
	sigma_k : list of tuple
		The list where each tuple represent the sigma and the k of the trait
	lineage_index : int
		the index of parental lineage

	Returns
	-------
	int :
		The localities selected, if all localities are occupied return -1
	"""
	# define parameters
	sigma = []
	k = []
	for i in sigma_k:
		sigma.append(i[0])
		k.append(i[1])
	s_opt = genus['s_opt'][lineage_index]
	# acces to left localities
	currently_occuring = []
	for i in range(len(genus['lineage'])):
		if genus['lineage'][i] == lineage:
			currently_occuring.append(genus['geography'][i])
	currently_occuring = set(currently_occuring)
	left_localities = [i for i in range(nLocalities) if i not in currently_occuring]
	# define weight
	w = []

	for locality in left_localities:

		wi = np.zeros(len(sigma_k))
		for i in range(len(sigma_k)):
			wi[i] = exp(-sigma[i] * (np.absolute(s_opt[i] - s_deme[locality][i]) ** k[i]))
		w.append(sum(sigma * wi) / sum(sigma))
	# pick a localities
	if len(left_localities) >0:
		res = choices(left_localities, weights=w)[0]
	else:
		res = -1
	return(res)

def get_geography_weighted(s_deme, sigma_k, genus, nLocalities, lineage_index):
	"""
	Return one geography where the new lineage will be present

	Parameters
	----------
	s_deme : dict
		All the values of s_deme
	sigma_k : list of tuple
		The list where each tuple represent the sigma and the k of the trait
	genus : dict
		The dict use in the simulation to represent the genus
	nLocalities : int
		Number of localities in the simulation
	lineage_index : int
		the index of parental lineage

	Returns
	-------
	int:
		A localities
	"""
	sigma = []
	k = []
	for i in sigma_k:
		sigma.append(i[0])
		k.append(i[1])
	s_opt = genus['s_opt'][lineage_index]

	w = []

	for deme in s_deme.keys():
		wi = np.zeros(len(sigma_k))
		for i in range(len(sigma_k)):
			wi[i] = exp(-sigma_k[i][0] * (np.absolute(s_opt[i] - s_deme[deme][i])**sigma_k[i][1]))
		w.append(sum(sigma * wi) / sum(sigma))

	demes = [i for i in range(nLocalities)]
	return (choices(demes, weights=w)[0])
def get_index_living_lineages(n_non_extincted, genus):
	index_of_0 = [None] * int(n_non_extincted)
	posi = 0
	count = 0
	while index_of_0.count(None) != 0:
		temp = genus['extinction_statu'].index(0, posi)
		index_of_0[count] = temp
		count += 1
		posi = temp + 1
	return index_of_0

def simulation(rate, proportion, simulation_name, result_repertory, verbose, sigma, k, s_deme_csv=None, nb_localities=10,
			   nb_generations = 600, nb_max_genus = 2000, nb_max_lineages_in_genus = 50):
	"""
	Launch a simulation of the evolution of diversity of life, using speciation and extinction.

	Parameters
	----------
	rate : list
		The list with the rate E1, E2, E3, D1, D2, D3 and C1 (in this order).
	proportion :
		The list with the proportion of lineage to extinct with E2 (pE2) and E3 (pE3) (in this order).
	simulation_name : str
		The name of the simulation. Use to name the output file.
	result_repertory : str
		The path of the directory where the results will be store.
	verbose : bool
		Active the verbose mode if True.
	sigma : list
		The list with the value sigma of each trait. Take care to use the same order in k and in the s_deme_csv.
	k : list
		The list with the value k of each trait. Take care to use the same order in sigma and in the s_deme_csv.
	s_deme_csv : str or None
		The path of a csv file with each line is a locality and each column the s_deme of a trait in this locality.
		If None the s_deme are generated randomly and can take the value 0 or 1.
		None by default.
	nb_localities : int
		Number of localities. 10 by default
	nb_generations : int
		Max number of generation in the simulation. 600 by default.
	nb_max_genus : int
		Max number of genus in the simulation. 2,000 by default.
	nb_max_lineages_in_genus : int
		Max number of lineage for each genus. 50 by default.

	Returns
	-------
	None
	"""
	#### parameters used for trials
	nLocalities = nb_localities # number of geographic localities where species from different genus can be found
	nGenerations = nb_generations  # number of generations to simulate

	maxNGenus = nb_max_genus # maximum number of genus
	nMaxLineages_in_genus = nb_max_lineages_in_genus # a given genus can have a maximum number of species/lineages
	nMaxDemes_in_genus = nLocalities*nMaxLineages_in_genus

	E1, E2, E3, D1, D2, D3, C1 = rate
	pE2, pE3 = proportion

	assert len(sigma) == len(k), "The list sigma and the list k didn't have the same length"
	nTrait = len(sigma)
	sigma_k = []
	for i in range(len(sigma)):
		sigma_k.append((sigma[i], k[i]))

	scenar = simulation_name.split('_')[0]
	with open(f'{result_repertory}/{simulation_name}.par', 'w') as fileIn:
		fileIn.write(f'{scenar}\t{simulation_name}\t{E1}\t{E2}\t{E3}\t{pE2}\t{pE3}\t{D1}\t{D2}\t{C1}\t{nLocalities}\t{nGenerations}\t{maxNGenus}\t{nMaxLineages_in_genus}\t{nTrait}\n')

	living_genus = [] # list of genus which 1) were born and 2) not dead yet

	s_deme = {}
	if s_deme_csv is None:
		for deme in range(nLocalities):
			s_deme[deme] = randint(0, 2, size=nTrait)  # randint for the moment
		type_of_trait = ['d'] * nTrait
	else:
		type_of_trait = []
		df = pd.read_csv(s_deme_csv, sep=',', header=None)
		assert df.shape[1] == nTrait, "Not enough trait's values in the s_deme csv"
		for i in range(df.shape[0]):
			s_deme[i] = df.iloc[i,:].to_list()
		for i in range(df.shape[1]):
			if df.iloc[:,i].isin([0, 1]).all():
				type_of_trait.append('d')
			else:
				type_of_trait.append('c')

	record = {}
	for genus in range(maxNGenus):
		record[genus] = {} # record['genus']
		record[genus]['nLineages'] = 0
		record[genus]['lineage'] = [nan] * nMaxDemes_in_genus
		record[genus]['oldest'] = [nan] * nMaxDemes_in_genus
		record[genus]['youngest'] = [nan] * nMaxDemes_in_genus
		record[genus]['geography'] = [nan] * nMaxDemes_in_genus
		record[genus]['extinction_statu'] = [nan] * nMaxDemes_in_genus
		record[genus]['s_opt'] = [nan] * nMaxDemes_in_genus

	living_genus.append(0)

	record[0]['nLineages'] +=1
	record[0]['lineage'][0] = 0
	record[0]['oldest'][0] = nGenerations
	record[0]['youngest'][0] = nGenerations
	record[0]['geography'][0] = 0
	record[0]['extinction_statu'][0] = 0
	record[0]['s_opt'][0] = randint(0, 2, size=nTrait) # randint for the moment

	record_extincted = {}
	record_extincted['genus'] = []
	record_extincted['lineage'] = []
	record_extincted['geography'] = []
	record_extincted['youngest'] = []	
	record_extincted['oldest'] = []
	record_extincted['extinction_statu'] = []
	record_extincted['s_opt'] = []

	record_extincted['genus'] = [nan]*4000000
	record_extincted['lineage'] =[nan]*4000000 
	record_extincted['geography'] = [nan]*4000000
	record_extincted['youngest'] = 	[nan]*4000000
	record_extincted['oldest'] = [nan]*4000000
	record_extincted['extinction_statu'] = [nan]*4000000
	record_extincted['s_opt'] = [nan]*4000000
	nExtincted_lineages = 0

	mass_extinctions_genus = []
	mass_extinctions_times = []
	mass_extinctions_kind = []
	mass_extinctions_effect = []
	mass_extinctions_deme = []

	last_created_genus = 0
	for time in range(nGenerations)[::-1]:
		if verbose:
			print('Generation {0}: {1} extincted lineages'.format(time, nExtincted_lineages))
		# origination of a new genus
		test_new_genus = binomial(n=1, p=D2, size=1)[0]
		if test_new_genus:
			last_created_genus += 1
			new_genus = last_created_genus
			living_genus.append(new_genus)
			record[new_genus]['nLineages'] = 1

			record[new_genus]['lineage'][0] = 0
			record[new_genus]['oldest'][0] = time
			record[new_genus]['youngest'][0] = time
			new_geography = randint(low=0, high=nLocalities, size=1)[0]
			record[new_genus]['geography'][0] = new_geography
			record[new_genus]['extinction_statu'][0] = 0
			record[new_genus]['s_opt'][0] = s_deme[new_geography]
			
		# loop over different genus for the origination of new lineage
	#	for genus in record.keys():
		for genus in living_genus:
			print(f'genus = {genus}/{len(living_genus)-1}')
			# only treats genus that are not fully extincted
			if record[genus]['lineage'].count(nan)!=len(record[genus]['lineage']):
				# origination of a new lineage
				n_non_extincted_demes = len(record[genus]['extinction_statu']) - nansum(record[genus]['extinction_statu']) - record[genus]['extinction_statu'].count(nan)
				n_new_lineages = sum(binomial(n=1, p=D1, size=int(n_non_extincted_demes)))
				if n_new_lineages>0:
					index_of_living = get_index_living_lineages(n_non_extincted_demes, record[genus])
					parental_index = choice(a=index_of_living, size=n_new_lineages, replace=False)
					for lineage in parental_index:
						# if the genus is not "full"
						if record[genus]['lineage'].count(nan)>0:
							newPos = record[genus]['lineage'].index(nan)
							record[genus]['nLineages'] += 1
							record[genus]['lineage'][newPos] = record[genus]['nLineages']-1
							record[genus]['oldest'][newPos] = time
							record[genus]['youngest'][newPos] = time
							new_geography = get_geography_weighted(s_deme, sigma_k, record[genus], nLocalities, lineage)
							record[genus]['geography'][newPos] = new_geography # the new lineage appears in a new geography
							record[genus]['extinction_statu'][newPos] = 0
							record[genus]['s_opt'][newPos] = s_deme[new_geography]


				n_non_extincted_demes = int(len(record[genus]['extinction_statu']) - nansum(record[genus]['extinction_statu']) - record[genus]['extinction_statu'].count(nan))
				n_new_lineages = sum(binomial(n=1, p=D3, size=n_non_extincted_demes))
				if n_new_lineages > 0:
					index_of_living = get_index_living_lineages(n_non_extincted_demes, record[genus])
					index_to_change = choice(a=index_of_living, size=n_new_lineages, replace=False)

					for lineage in index_to_change:
						record[genus]['extinction_statu'][lineage] = 1
						record_extincted['genus'][nExtincted_lineages] = genus
						record_extincted['lineage'][nExtincted_lineages] = record[genus]['lineage'][lineage]
						record_extincted['geography'][nExtincted_lineages] = record[genus]['geography'][lineage]
						record_extincted['youngest'][nExtincted_lineages] = time
						record_extincted['oldest'][nExtincted_lineages] = record[genus]['oldest'][lineage]
						record_extincted['extinction_statu'][nExtincted_lineages] = record[genus]['extinction_statu'][lineage]
						record_extincted['s_opt'][nExtincted_lineages] = record[genus]['s_opt'][lineage]
						nExtincted_lineages += 1

						record[genus]['lineage'][lineage] = nan
						record[genus]['oldest'][lineage] = nan
						record[genus]['youngest'][lineage] = nan
						# record[genus]['geography'][lineage] = nan
						record[genus]['extinction_statu'][lineage] = nan

						record[genus]['nLineages'] += 1

						record[genus]['lineage'][lineage] = record[genus]['nLineages'] - 1
						record[genus]['oldest'][lineage] = time
						record[genus]['youngest'][lineage] = time
						# record[genus]['geography'][lineage] = record[genus]['geography'][lineage]
						record[genus]['extinction_statu'][lineage] = 0

				# colonization of a new geographic area
				n_non_extincted_demes = len(record[genus]['extinction_statu']) - nansum(record[genus]['extinction_statu']) - record[genus]['extinction_statu'].count(nan)
				n_colonizing_lineages = sum(binomial(n=1, p=C1, size=int(n_non_extincted_demes)))
				
				if n_colonizing_lineages>0:
					copied_deme = [ i for i in range(len(record[genus]['lineage'])) if isnan(record[genus]['lineage'][i])==False ]
					copied_deme = choice(a=copied_deme, size=n_colonizing_lineages, replace=True) # a deme can colonize multiple new localities at the same time
					for colonizer in copied_deme:
						if record[genus]['lineage'].count(nan)>0:
							new_geography = get_geography_weighted_colonization(record[genus], record[genus]['lineage'][colonizer], nLocalities, s_deme, sigma_k, colonizer)
							if new_geography != -1:
								print('~~~~~~~~~~~~')
								print(colonizer)
								print(record[genus]['lineage'])
								print(record[genus]['geography'])
								print(record[genus]['s_opt'])
								print(s_deme)
								print('------------')
								newPos = record[genus]['lineage'].index(nan)
								record[genus]['lineage'][newPos] = record[genus]['lineage'][colonizer]
								record[genus]['oldest'][newPos] = time
								record[genus]['youngest'][newPos] = time
								record[genus]['geography'][newPos] = new_geography
								record[genus]['extinction_statu'][newPos] = 0
								record[genus]['s_opt'][newPos] = s_deme[new_geography]
								print(record[genus]['lineage'])
								print(record[genus]['geography'])
								print(record[genus]['s_opt'])
				
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
							
							record_extincted['genus'][nExtincted_lineages] = genus
							record_extincted['lineage'][nExtincted_lineages] = record[genus]['lineage'][death] 
							record_extincted['geography'][nExtincted_lineages] = record[genus]['geography'][death] 
							record_extincted['youngest'][nExtincted_lineages] = time 
							record_extincted['oldest'][nExtincted_lineages] = record[genus]['oldest'][death] 
							record_extincted['extinction_statu'][nExtincted_lineages] = record[genus]['extinction_statu'][death]
							record_extincted['s_opt'][nExtincted_lineages] = record[genus]['s_opt'][death]
							nExtincted_lineages += 1
							
							record[genus]['lineage'][death] = nan
							record[genus]['oldest'][death] = nan 
							record[genus]['youngest'][death] = nan 
							record[genus]['geography'][death] = nan 
							record[genus]['extinction_statu'][death] = nan
							record[genus]['s_opt'][death] = nan

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

							record_extincted['genus'][nExtincted_lineages] = genus
							record_extincted['lineage'][nExtincted_lineages] = record[genus]['lineage'][i] 
							record_extincted['geography'][nExtincted_lineages] = record[genus]['geography'][i] 
							record_extincted['youngest'][nExtincted_lineages] = record[genus]['youngest'][i]
							record_extincted['oldest'][nExtincted_lineages] = record[genus]['oldest'][i] 
							record_extincted['extinction_statu'][nExtincted_lineages] = record[genus]['extinction_statu'][i] 
							record_extincted['s_opt'][nExtincted_lineages] = record[genus]['s_opt'][i]
							nExtincted_lineages += 1

							record[genus]['lineage'][i] = nan
							record[genus]['oldest'][i] = nan
							record[genus]['youngest'][i] = nan
							record[genus]['geography'][i] = nan
							record[genus]['extinction_statu'][i] = nan
							record[genus]['s_opt'][i] = nan
		
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
						record_extincted['s_opt'][nExtincted_lineages] = record[genus]['s_opt'][index]
						nExtincted_lineages += 1

						record[genus]['lineage'][index] = nan
						record[genus]['oldest'][index] = nan
						record[genus]['youngest'][index] = nan
						record[genus]['geography'][index] = nan
						record[genus]['extinction_statu'][index] = nan
						record[genus]['s_opt'][index] = nan

			

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
