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
import argparse
from math import exp
import random
import math

def poisson_pypy(lam, size):
	"""
	Genere des valeurs aléatoires suivant une distribution de Poisson.

	Args:
		lam (float): Parametre lambda de la distribution de Poisson, le taux moyen d'occurrences.
		size (int): Nombre de valeurs aleatoires à générer.

	Returns:
		list: Liste contenant les valeurs générées suivant la distribution de Poisson.
	"""
	res = [0] * size
	L = math.exp(-lam)

	for i in range(size):
		k = 0
		p = 1

		# Genere le nombre de succes jusqu'a ce que le produit des aleas soit inferieur à L
		while p > L:
			p *= random.random()
			k += 1

		res[i] = k - 1
	return res


def binomial_pypy(n, p, size):
	""" Genere des valeurs aleatoires suivant une distribution binomiale.

	Args:
		n (int): Le nombre d'essais.
		p (float): La probabilité de succes de chaque essai.
		size (int): Nombre de fois que la distribution doit etre generee.

	Returns:
		list: Liste contenant les resultats des distributions binomiales.
	"""
	res = []
	for _ in range(size):
		res.append(sum(1 for _ in range(n) if random.random() < p))
	return res


def choice_pypy(a, size=None, replace=True, p=None):
	"""
	Selectionne des elements aléatoires d'une liste. Il y a des options pour des probabilites non uniformes (p)
	mais egalement avec ou sans remplacement des elements echantillonnes (replace).

	Args:
		a (list): La liste source a partir de laquelle les elements seront tires.
		size (int, optional): Le nombre d'elements à choisir. Si None, un seul élément est tire.
		replace (bool, optional): Si True, le tirage se fait avec remplacement. Si False, sans remplacement.
		p (list of float, optional): Une liste de probabilites correspondant à chaque élément dans `a`. Les probabilites doivent sommer à 1.

	Returns:
		list: une liste d'elements tires aleatoirement.
	"""
	if not replace and size is not None and size > len(a):
		raise ValueError("Cannot take a larger sample than population when 'replace=False'")

	if p is not None:
		if len(p) != len(a):
			raise ValueError("Probabilities p and array a must have the same size")
		if not math.isclose(sum(p), 1):
			raise ValueError("Probabilities p must sum to 1")
	
	# Choix unique ou multiple
	if size is None:
		if p is None:
			return random.choice(a)
		else:
			return random.choices(a, weights=p, k=1)[0]
	else:
		if p is None:
			if replace:
				return [random.choice(a) for _ in range(size)]
			else:
				return random.sample(a, size)
		else:
			return random.choices(a, weights=p, k=size)


def randint_pypy(low, high=None, size=None):
	"""
	Retourne des entiers aleatoires dans l'intervalle [low, high) ou [0, low).

	Args:
		low (int): limite inferieure (inclusive). Si high est None, low est utilise comme limite superieure.
		high (int, optional): Borne superieure (exclusive).
		size (int or tuple of ints, optional): La taille de la sortie. Si None, retourne un seul int.

	Returns:
		int or list: Un entier ou une liste d'entiers aléatoires.
	"""
	# Gerer le cas ou "high" n'est pas fourni
	if high is None:
		high = low
		low = 0

	# Generer un seul nombre aleatoire si "size" n'est pas specifie
	if size is None:
		return random.randint(low, high-1)

	# Generer une liste ou un tableau de nombres aleatoires
	if isinstance(size, int):  # Gerer la taille unidimensionnelle
		return [random.randint(low, high-1) for _ in range(size)]
	elif isinstance(size, tuple):  # Gerer la taille multidimensionnelle
		def generate_multidimensional_random(sizes, dim=0):
			if dim < len(sizes) - 1:
				return [generate_multidimensional_random(sizes, dim + 1) for _ in range(sizes[dim])]
			else:
				return [random.randint(low, high-1) for _ in range(sizes[dim])]

		return generate_multidimensional_random(size)
	else:
		raise TypeError("size should be an int or a tuple of ints")


def nonesum(iterable):
	"""
	Calcule la somme des elements d'une liste, en ignorant les 'None'.

	Args:
		iterable (iterable): Une liste ou un autre iterable contenant des nombres et potentiellement des 'None'.

	Returns:
		float: La somme des elements qui ne sont pas 'None'.
	"""
	total = 0
	for element in iterable:
		if element is not None:
			total += element
	return total


def isnone(value):
	"""
	Teste si la valeur est None.

	Args:
		value (any): La valeur a tester.

	Returns:
		bool: True si la valeur est None, sinon False.
	"""
	return value is None

def recreate_record(record, current_genus, max_living_genus, nMaxDemes_in_genus):
    still_alive = []
    for genus in current_genus:
        if record[genus]['lineage'].count(None)!=len(record[genus]['lineage']):
            still_alive.append(genus)

    new_record = {}
    for genus in still_alive:
        new_record[genus]=record[genus]
    
    last = current_genus[-1]
    need_to_add = last + (max_living_genus - len(new_record))
    for genus in range(last+1, need_to_add+1):
        new_record[genus] = {}
        new_record[genus]['nLineages'] = 0
        new_record[genus]['lineage'] = [None] * nMaxDemes_in_genus
        new_record[genus]['oldest'] = [None] * nMaxDemes_in_genus
        new_record[genus]['youngest'] = [None] * nMaxDemes_in_genus
        new_record[genus]['geography'] = [None] * nMaxDemes_in_genus
        new_record[genus]['extinction_statu'] = [None] * nMaxDemes_in_genus
    return new_record, still_alive, len(still_alive)

def div_competition(diversification, k, current_nb_species, carrying_capacity):
	density = current_nb_species / carrying_capacity
	alpha = exp(-k*(density**2))
	return diversification*alpha

def extinction_competition(extinction, k, current_nb_species, carrying_capacity):
	density = current_nb_species / carrying_capacity
	alpha = 1-(exp(-k*(density**2)))
	return extinction*alpha

def get_geography(genus, lineage, nLocalities, species_in_deme):
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

	Returns
	-------
	int :
		The localities selected, if all localities are occupied return -1
	"""
	currently_occuring = []
	for i in range(len(genus['lineage'])):
		if genus['lineage'][i]==lineage:
			currently_occuring.append(genus['geography'][i])
	currently_occuring = set(currently_occuring)
	left_localities = [ i for i in range(nLocalities) if i not in currently_occuring ]
	deme_full = {}
	for i in left_localities:
		if species_in_deme[i].count(None) == 0:
			deme_full[i]=True
		else:
			deme_full[i]=False	
	left_localities = [i for i in left_localities if not deme_full[i]]
	if len(left_localities) >0:
		res = choice_pypy(a=left_localities, size=1, replace=False)[0]
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

def extinction(new_geography, species_in_deme, carrying_capacity, k, extinction, competition, genera, nb_test):
	count_species = len(species_in_deme[new_geography]) - species_in_deme[new_geography].count(None)
	if competition == 'universal':
		new_E = extinction_competition(extinction=extinction, k=k, current_nb_species=count_species, carrying_capacity=carrying_capacity)
		test_extinction = binomial_pypy(n=1, p=new_E, size=nb_test)
		return test_extinction
	elif competition == 'in_genera':
		species_in_genera = len(species_in_deme) 
		new_E = extinction_competition(extinction=extinction, k=k, current_nb_species=species_in_genera, carrying_capacity=carrying_capacity)
		test_extinction = binomial_pypy(n=1, p=new_E, size=nb_test)
		return test_extinction
	elif competition is None:
		test_extinction = binomial_pypy(n=1, p=extinction, size=nb_test)
		return test_extinction

def speciation(new_geography, species_in_deme, carrying_capacity, k, diversification, competition):
	count_species = len(species_in_deme[new_geography]) - species_in_deme[new_geography].count(None)
	if count_species != carrying_capacity:
		if competition == 'universal':
			new_D = div_competition(diversification=diversification, k=k, current_nb_species=count_species, carrying_capacity=carrying_capacity)
			test_new_genus = binomial_pypy(n=1, p=new_D, size=1)[0]
			return test_new_genus
		elif competition == 'in_genera':
			species_in_genera = len(species_in_deme)
			new_D = div_competition(diversification=diversification, k=k, current_nb_species=species_in_genera, carrying_capacity=carrying_capacity)
			test_new_genus = binomial_pypy(n=1, p=new_D, size=1)[0]
			return test_new_genus
		elif competition is None:
			test_new_genus = binomial_pypy(n=1, p=diversification, size=1)[0]
			return test_new_genus
	else:
		return 0
		
	

def simulation(rate, proportion, simulation_name, result_repertory, verbose, carrying_capacity, k, speciation_competition, extinction_competition,
			   nb_localities=10, nb_generations = 600, nb_max_genus = 2000, nb_max_lineages_in_genus = 50):
	"""
	Launch a simulation of the evolution of diversity of life, using speciation and extinction.

	Parameters
	----------
	rate : list
		The list with the rate E1, E2, E3, D1, D2 and C1 (in this order).
	proportion :
		The list with the proportion of lineage to extinct with E2 (pE2) and E3 (pE3) (in this order).
	simulation_name : str
		The name of the simulation. Use to name the output file.
	result_repertory : str
		The path of the directory where the results will be store.
	verbose : bool
		Active the verbose mode if True.
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

	E1, E2, E3, D1, D2, C1 = rate
	pE2, pE3 = proportion

	if speciation_competition == 0:
		D1_compet = None
		D2_compet = None
	elif speciation_competition == 1:
		D1_compet = None
		D2_compet = 'in_genera'
	elif speciation_competition == 2:
		D1_compet = 'universal'
		D2_compet = 'universal'
	
	if extinction_competition == 0:
		E1_comp = None
		E2_comp = None
	elif extinction_competition == 1:
		E1_comp = 'universal'
		E2_comp = None
	elif extinction_competition == 2:
		E1_comp = None
		E2_comp = 'in_genera'
	
	scenar = simulation_name.split('_')[0]
	with open(f'{result_repertory}/{simulation_name}.par', 'w') as fileIn:
		fileIn.write(f'{scenar}\t{simulation_name}\t{E1}\t{E2}\t{E3}\t{pE2}\t{pE3}\t{D1}\t{D2}\t{C1}\t{nLocalities}\t{nGenerations}\t{maxNGenus}\t{nMaxLineages_in_genus}\n')

	living_genus = [] # list of genus which 1) were born and 2) not dead yet

	record = {}
	for genus in range(maxNGenus):
		record[genus] = {} 
		record[genus]['nLineages'] = 0
		record[genus]['lineage'] = [None] * nMaxDemes_in_genus
		record[genus]['oldest'] = [None] * nMaxDemes_in_genus
		record[genus]['youngest'] = [None] * nMaxDemes_in_genus
		record[genus]['geography'] = [None] * nMaxDemes_in_genus
		record[genus]['extinction_statu'] = [None] * nMaxDemes_in_genus

	# setup the list which will contain the lineages in a deme at each time.
	deme_species = [([None]*carrying_capacity[i]) for i in range(nLocalities)]

	living_genus.append(0)

	record[0]['nLineages'] +=1
	record[0]['lineage'][0] = 0
	record[0]['oldest'][0] = nGenerations
	record[0]['youngest'][0] = nGenerations
	record[0]['geography'][0] = 0
	record[0]['extinction_statu'][0] = 0

	deme_species[0][0] = (0, 0) # (genus, index of this species in record)

	record_extincted = {}
	record_extincted['genus'] = []
	record_extincted['lineage'] = []
	record_extincted['geography'] = []
	record_extincted['youngest'] = []	
	record_extincted['oldest'] = []
	record_extincted['extinction_statu'] = []

	record_extincted['genus'] = [None]*4000000
	record_extincted['lineage'] =[None]*4000000 
	record_extincted['geography'] = [None]*4000000
	record_extincted['youngest'] = 	[None]*4000000
	record_extincted['oldest'] = [None]*4000000
	record_extincted['extinction_statu'] = [None]*4000000
	nExtincted_lineages = 0

	mass_extinctions_genus = []
	mass_extinctions_times = []
	mass_extinctions_kind = []
	mass_extinctions_effect = []
	mass_extinctions_deme = []

	last_created_genus = 0
	counter = 1
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
			print(f'genus : {counter}')
		n_living_genus = 0
		for genus in living_genus:
			if record[genus]['lineage'].count(None)!=len(record[genus]['lineage']):
				n_living_genus += 1
		# loop over different genus for the origination of new lineage

		for i in range(len(living_genus)):
			genus = living_genus[i]
			# origination of a new genus
			if n_living_genus > 0:
				if record[genus]['lineage'].count(None)!=len(record[genus]['lineage']):
					if last_created_genus < list(record.keys())[-1]:
						new_geography = randint_pypy(low=0, high=nLocalities, size=1)[0]
						test_new_genus = speciation(new_geography=new_geography, species_in_deme=deme_species, k=k, carrying_capacity=carrying_capacity[new_geography],
										diversification=D2, competition=D2_compet)
						if test_new_genus:
							counter += 1
							last_created_genus += 1
							new_genus = last_created_genus
							living_genus.append(new_genus)
							record[new_genus]['nLineages'] = 1
							
							record[new_genus]['lineage'][0] = 0
							record[new_genus]['oldest'][0] = time
							record[new_genus]['youngest'][0] = time
							record[new_genus]['geography'][0] = new_geography
							record[new_genus]['extinction_statu'][0] = 0

							index_in_deme = deme_species[new_geography].index(None) # add the new genus to deme current species
							deme_species[new_geography][index_in_deme] = (new_genus, 0)
						else:
							pass
					else: #if last_created_genus >= len(record)-1:
						previous_lenght = len(living_genus)
						record, living_genus, n_living_genus = recreate_record(record=record, current_genus=living_genus, max_living_genus=maxNGenus, nMaxDemes_in_genus=nMaxDemes_in_genus)
						if len(living_genus) < previous_lenght:
							break

			else:				
				new_geography = randint_pypy(low=0, high=nLocalities, size=1)[0]
				test_new_genus = speciation(new_geography=new_geography, species_in_deme=deme_species, k=k, carrying_capacity=carrying_capacity[new_geography],
								diversification=D2, competition=D2_compet)
				if test_new_genus:
					if last_created_genus >= list(record.keys())[-1]:
						record, living_genus, n_living_genus = recreate_record(record=record, current_genus=living_genus, max_living_genus=maxNGenus, nMaxDemes_in_genus=nMaxDemes_in_genus)
					counter += 1
					n_living_genus += 1
					last_created_genus += 1
					new_genus = last_created_genus
					living_genus.append(new_genus)
					record[new_genus]['nLineages'] = 1
					
					record[new_genus]['lineage'][0] = 0
					record[new_genus]['oldest'][0] = time
					record[new_genus]['youngest'][0] = time
					record[new_genus]['geography'][0] = new_geography
					record[new_genus]['extinction_statu'][0] = 0

					index_in_deme = deme_species[new_geography].index(None) # add the new genus to deme current species
					deme_species[new_geography][index_in_deme] = (new_genus, 0)
					break
				else:
					pass
					
		
		for genus in living_genus:
			# only treats genus that are not fully extincted
			if record[genus]['lineage'].count(None)!=len(record[genus]['lineage']):
				# origination of a new lineage
				n_non_extincted_demes = len(record[genus]['extinction_statu']) - nonesum(record[genus]['extinction_statu']) - record[genus]['extinction_statu'].count(None)
				new_geography = randint_pypy(low=0, high=nLocalities, size=int(n_non_extincted_demes))

				for geography in new_geography:				
					n_new_lineages = speciation(new_geography=geography, species_in_deme=deme_species, k=k, carrying_capacity=carrying_capacity[geography],
								diversification=D1, competition=D1_compet)
					if n_new_lineages:
						# if the genus is not "full"
						if record[genus]['lineage'].count(None)>0:
							newPos = record[genus]['lineage'].index(None)
							record[genus]['nLineages'] += 1
							
							record[genus]['lineage'][newPos] = record[genus]['nLineages']-1
							record[genus]['oldest'][newPos] = time
							record[genus]['youngest'][newPos] = time
							
							record[genus]['geography'][newPos] = geography # the new lineage appears in a new geography
							record[genus]['extinction_statu'][newPos] = 0

							index_in_deme = deme_species[geography].index(None) # add the new genus to deme current species
							deme_species[geography][index_in_deme] = (genus, newPos)

				# colonization of a new geographic area
				n_non_extincted_demes = len(record[genus]['extinction_statu']) - nonesum(record[genus]['extinction_statu']) - record[genus]['extinction_statu'].count(None)
				n_colonizing_lineages = sum(binomial_pypy(n=1, p=C1, size=int(n_non_extincted_demes)))
				
				if n_colonizing_lineages>0:
					copied_deme = [ i for i in range(len(record[genus]['lineage'])) if isnone(record[genus]['lineage'][i])==False ]
					copied_deme = choice_pypy(a=copied_deme, size=n_colonizing_lineages, replace=True) # a species can colonize multiple new localities at the same time
					for colonizer in copied_deme:
						if record[genus]['lineage'].count(None)>0:
							new_geography = get_geography(record[genus], record[genus]['lineage'][colonizer], nLocalities, species_in_deme=deme_species)
							if new_geography != -1:
								newPos = record[genus]['lineage'].index(None)
								record[genus]['lineage'][newPos] = record[genus]['lineage'][colonizer]
								record[genus]['oldest'][newPos] = time
								record[genus]['youngest'][newPos] = time
								record[genus]['geography'][newPos] = new_geography
								record[genus]['extinction_statu'][newPos] = 0

								index_in_deme = deme_species[new_geography].index(None)  # add the new genus to deme current species
								deme_species[new_geography][index_in_deme] = (genus, newPos)

		for geography in range(len(deme_species)):
			# extinction of a given deme (lineage in a locality)
			current_species = [i for i in deme_species[geography] if not isnone(i)]
			current_species_index = [i for i in range(len(current_species))]
			if len(current_species)>0:
				n_non_extincted_species = len(current_species)
				n_lineages_to_kill = sum(extinction(new_geography=geography, species_in_deme=deme_species, carrying_capacity=carrying_capacity[geography],
									k=k, extinction=E1, competition=E1_comp, genera=None, nb_test=n_non_extincted_species))
				#n_lineages_to_kill = sum(binomial(n=1, p=E1, size=int(n_non_extincted_demes)))
				
				if n_lineages_to_kill>n_non_extincted_demes:
					n_lineages_to_kill=n_non_extincted_demes

				if n_lineages_to_kill>0:
						killed_deme = choice_pypy(a=current_species_index, size=int(n_lineages_to_kill), replace=False) # a deme cannot be extincted more that one time
						for i in killed_deme:
							genus = current_species[i][0]
							death = current_species[i][1]
							record[genus]['extinction_statu'][death] = 1
							record_extincted['genus'][nExtincted_lineages] = genus
							record_extincted['lineage'][nExtincted_lineages] = record[genus]['lineage'][death]
							record_extincted['geography'][nExtincted_lineages] = record[genus]['geography'][death]
							record_extincted['youngest'][nExtincted_lineages] = time
							record_extincted['oldest'][nExtincted_lineages] = record[genus]['oldest'][death]
							record_extincted['extinction_statu'][nExtincted_lineages] = record[genus]['extinction_statu'][death]
							nExtincted_lineages += 1
							index_in_deme = deme_species[record[genus]['geography'][death]].index((genus, death))  # delete the lineage to deme current species
							deme_species[record[genus]['geography'][death]][index_in_deme] = None

							record[genus]['lineage'][death] = None
							record[genus]['oldest'][death] = None
							record[genus]['youngest'][death] = None
							record[genus]['geography'][death] = None
							record[genus]['extinction_statu'][death] = None

		# loop over different genus for the extinction of a genus
		for genus in living_genus:
			if record[genus]['lineage'].count(None)!=len(record[genus]['lineage']):
				#species_in_genera = Counter(i[0] for i in species_in_deme[new_geography])[genera]
				demes_of_genus = []
				for i in range(nLocalities):
					demes_of_genus.append([])
				for i in range(len(deme_species)):
					for j in range(len(deme_species[i])):
						if not isnone(deme_species[i][j]) and deme_species[i][j][0] == genus:
							demes_of_genus[i].append(j)

				for index in range(len(demes_of_genus)):
					# number of non extincted lineages
					current_demes = demes_of_genus[index]
					if len(current_demes) > 0:
						n_non_extincted_demes = len(current_demes)
						# test extinction
						test_extinction_genus = extinction(new_geography=index, species_in_deme=demes_of_genus, carrying_capacity=carrying_capacity[index],
										k=k, extinction=E2, competition=E2_comp, genera=genus, nb_test=1)[0]
						
						#binomial(n=1, p=E2, size=1)[0]
						if test_extinction_genus:							
							nDemes_to_extinct = poisson_pypy(lam=pE2*n_non_extincted_demes, size=1)[0]
							
							if nDemes_to_extinct>n_non_extincted_demes:
								nDemes_to_extinct=n_non_extincted_demes

							if nDemes_to_extinct>0:
								mass_extinctions_genus.append(genus)
								mass_extinctions_deme.append(i)
								mass_extinctions_times.append(time)
								mass_extinctions_kind.append('E2')
								mass_extinctions_effect.append(nDemes_to_extinct)
								# list of demes to extinct
								demes_to_extinct = choice_pypy( a=current_demes, size=nDemes_to_extinct, replace=False)
								for i in demes_to_extinct:
									posi = deme_species[index][i][1]
									
									record[genus]['extinction_statu'][posi] = 1
									record[genus]['youngest'][posi] = time

									record_extincted['genus'][nExtincted_lineages] = genus
									record_extincted['lineage'][nExtincted_lineages] = record[genus]['lineage'][posi] 
									record_extincted['geography'][nExtincted_lineages] = record[genus]['geography'][posi] 
									record_extincted['youngest'][nExtincted_lineages] = record[genus]['youngest'][posi]
									record_extincted['oldest'][nExtincted_lineages] = record[genus]['oldest'][posi] 
									record_extincted['extinction_statu'][nExtincted_lineages] = record[genus]['extinction_statu'][posi] 
									nExtincted_lineages += 1

									index_in_deme = deme_species[record[genus]['geography'][posi]].index((genus, posi))  # delete the lineage to deme current species
									deme_species[record[genus]['geography'][posi]][index_in_deme] = None
									
									record[genus]['lineage'][posi] = None
									record[genus]['oldest'][posi] = None
									record[genus]['youngest'][posi] = None
									record[genus]['geography'][posi] = None
									record[genus]['extinction_statu'][posi] = None
		# check the demes where an extinction can take place
		not_empty_deme = []
		for i in range(len(deme_species)):
			if deme_species[i].count(None) != len(deme_species[i]):
				not_empty_deme.append(i)

		if len(not_empty_deme) != 0:
			# loop over different deme for the extinction of a deme
			for deme in not_empty_deme:
				test_extinction_deme = binomial_pypy(n=1, p=E3, size=1)[0]
				if test_extinction_deme:
					# number of non extincted lineages
					non_extincted_lineages = [i for i in range(len(deme_species[deme])) if not isnone(deme_species[deme][i]) ]
					n_non_extincted_lineages = len(non_extincted_lineages)
					# number of extinctions
					nLineages_to_extinct = poisson_pypy(lam=pE3 * n_non_extincted_lineages, size=1)[0]
					if nLineages_to_extinct > n_non_extincted_lineages:
						nLineages_to_extinct = n_non_extincted_lineages

					if nLineages_to_extinct > 0:
						# list of demes to extinct
						mass_extinctions_genus.append(-1)  # not adapted to E3
						mass_extinctions_deme.append(deme)
						mass_extinctions_times.append(time)
						mass_extinctions_kind.append('E3')
						mass_extinctions_effect.append(nLineages_to_extinct)
						# list of lineages to extinct
						lineages_to_extinct = choice_pypy(a=non_extincted_lineages, size=nLineages_to_extinct, replace=False)
						for i in lineages_to_extinct:
							posi = deme_species[deme][i][1]
							genus = deme_species[deme][i][0]
							record_extincted['genus'][nExtincted_lineages] = genus
							record_extincted['lineage'][nExtincted_lineages] = record[genus]['lineage'][posi]
							record_extincted['geography'][nExtincted_lineages] = deme
							record_extincted['youngest'][nExtincted_lineages] = time
							record_extincted['oldest'][nExtincted_lineages] = record[genus]['oldest'][posi]
							record_extincted['extinction_statu'][nExtincted_lineages] = 1
							nExtincted_lineages += 1

							deme_species[deme][i] = None

							record[genus]['lineage'][posi] = None
							record[genus]['oldest'][posi] = None
							record[genus]['youngest'][posi] = None
							record[genus]['geography'][posi] = None
							record[genus]['extinction_statu'][posi] = None

	# output file
	## fossils
	outfile = open(f'{result_repertory}/simulated_fossils_{simulation_name}.txt', 'w')
	outfile.write('simulation\tspecimen\tgenus\tspecies\tgeography\tmax_ma\tmin_ma\textinction_statu\n')

	i=0 ## In the case or nExtincted_lineages is egal to 0.
	for i in range(nExtincted_lineages):
		outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t1\n'.format(simulation_name, i, record_extincted['genus'][i], record_extincted['lineage'][i], record_extincted['geography'][i], record_extincted['oldest'][i],  record_extincted['youngest'][i]))

	cnt = i
	## current species
	for genus in record:
		if len(record[genus]['lineage'])!=record[genus]['lineage'].count(None):
			for i in range(len(record[genus]['lineage'])):
				if isnone(record[genus]['lineage'][i])==False:
					cnt += 1
					outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t0\n'.format(simulation_name, cnt, genus, record[genus]['lineage'][i], record[genus]['geography'][i], record[genus]['oldest'][i], time))
				

	outfile.close()

	## mass extinction events
	with open(f'{result_repertory}/simulated_mass_extinctions_{simulation_name}.txt', 'w') as outfile:
		outfile.write('simulation\ttime\taffected genus\taffected deme\tkind of extinction\tnumber of extincted demes\n')
		for i in range(len(mass_extinctions_genus)):
			outfile.write(f'{simulation_name}\t{mass_extinctions_times[i]}\t{mass_extinctions_genus[i]}\t{mass_extinctions_deme[i]}\t{mass_extinctions_kind[i]}\t{mass_extinctions_effect[i]}\n')

if __name__ == "__main__":
	parser= argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--proba_E1', '-e1', help='Rates of extinction E1', type=float, required=True, action='store')
	parser.add_argument('--proba_E2', '-e2', help='Rates of extinction E2', type=float, required=True, action='store')
	parser.add_argument('--proba_E3', '-e3', help='Rates of extinction E3', type=float, required=True, action='store')
	parser.add_argument('--proba_D1', '-d1', help='Rates of speciation D1', type=float, required=True, action='store')
	parser.add_argument('--proba_D2', '-d2', help='Rates of speciation D2', type=float, required=True, action='store')
	parser.add_argument('--proba_C1', '-c1', help='Rates of colonisation C1', type=float, required=True, action='store')
	parser.add_argument('--proportion_E2', '-pe2', help='Proportion of genus mass extincted', type=float,required=True, action='store')
	parser.add_argument('--proportion_E3', '-pe3', help='Proportion of demes mass extincted', type=float,required=True, action='store')
	parser.add_argument('--name', '-n', help='Name of simulation', type=str, required=True, action='store')
	parser.add_argument('--output_dir', '-o', help='Output directory', type=str, required=True, action='store')
	parser.add_argument('--k', '-k', help='Values of k the response to the density', type=int, required=True, action='store')
	parser.add_argument('--carrying_capacity', '-c', help='Values of carrying capacity', type=int, required=True, nargs='+')
	#parser.add_argument('--carrying_capacity', '-c', help='Values of carrying capacity', type=int, required=True, action='store')
	parser.add_argument('--competition_on_speciation', '-cs',help='Code to activate the density dependent speciation:\n0: don\'t use density dependence.'
					 '\n1: apply density dependence on D2.\n2: apply density dependence on D1 and D2.\nSee more information on the README.txt'
					,type=int, choices=range(0,3), default=0, action='store')
	parser.add_argument('--competition_on_extinction', '-ce', help='Code to activate the density dependent extinction:\n0: don\'t use density dependence.'
					 '\n1: apply density dependence on E1.\n2: apply density dependence on E2.\nSee more information on the README.txt',
					 type=int, choices=range(0,3), default=0, action='store')
	parser.add_argument('--nb_localities', '-nl', help='Number of localities', type=int, default=10, action='store')
	parser.add_argument('--nb_generations', '-ng', help='Number of generations', type=int, default=600, action='store')
	parser.add_argument('--nb_max_genus', '-nmg', help='Max number of genus', type=int, default=7000, action='store')
	parser.add_argument('--nb_max_lineages', '-ngl', help='Max number of lineages in genus', type=int, default=50, action='store')
	parser.add_argument('--verbose', '-v', help='Display the generations', action='store_true')
	args = parser.parse_args()

	rates = [args.proba_E1, args.proba_E2, args.proba_E3, args.proba_D1, args.proba_D2, args.proba_C1]
	proportions = [args.proportion_E2, args.proportion_E3]
	k = args.k
	simulation(rate=rates,
			proportion=proportions,
			simulation_name=args.name,
			result_repertory=args.output_dir,
            verbose=args.verbose, 
			nb_localities=args.nb_localities, 
			nb_generations=args.nb_generations,
            k=k,
			carrying_capacity=args.carrying_capacity,
			speciation_competition=args.competition_on_speciation,
			extinction_competition=args.competition_on_extinction,
			nb_max_genus=args.nb_max_genus, 
			nb_max_lineages_in_genus=args.nb_max_lineages)