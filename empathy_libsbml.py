from libsbml import *
import xmltodict, json

def fetch_compartment(model,id):
	compartment = model.getCompartment(id)
	return compartment.getName()

def fetch_species(model,id):
	species = model.getSpecies(id)
	return species.getName(),fetch_compartment(model,species.getCompartment())

def unpackBags(bags):
	resources = []
	try:
		if len(bags) == 1:
			type = 1
			resource = bags['rdf:Bag']['rdf:li'][0]['@rdf:resource']
			resource = resource.replace('http://identifiers.org/','')
			resource = resource.replace('/',':',1)
			resources.append(str(resource))			
		else:
			type = 2
			for bag in bags:
				resource = bag['rdf:Bag']['rdf:li'][0]['@rdf:resource']
				resource = resource.replace('http://identifiers.org/','')
				resource = resource.replace('/',':',1)
				resources.append(str(resource))
	except Exception, e:
		print 'Something went wrong:', str(e), type, resources
	return resources

def getAnnotations(sbmlObject):
	return {'No annotations':''}
	'''
	try:
		annotationString = sbmlObject.getAnnotationString()
		annotationString = annotationString.replace("\n","")
		annotationString = annotationString.replace("  ","")
		annotationDict = xmltodict.parse( annotationString )
		annotationDict = json.loads( json.dumps( annotationDict ) )
		#pp.pprint(annotationDict)
		annotations = {}
		for bqbiol in annotationDict['annotation']['rdf:RDF']['rdf:Description']:
			if 'bqbiol:' in bqbiol:
				whack = bqbiol.split(':')
				resources = unpackBags(annotationDict['annotation']['rdf:RDF']['rdf:Description'][bqbiol])
				annotations[whack[1]] = resources
		return annotations
	except Exception, e:
		print 'Exception:', str(e)
		return {}
	'''

def getNotes(sbmlObject):
	try:
		notesString = sbmlObject.getNotesString()
		notesString = notesString.replace("\n","")
		notesString = notesString.replace("  ","")
		notesDict = xmltodict.parse( notesString )
		notesDict = json.loads( json.dumps( notesDict ) )
		notes = []
		for html in notesDict['notes']['html:body']['html:p']:
			notes.append( str(html) )
		return notes
	except:
		return []

def parseSBML(sbml):	
	print 'Read SBML...'
	try:
		reader = SBMLReader()
		document = reader.readSBMLFromString(str(sbml))
		model = document.getModel()
		print '...complete'
	except:
		return 'Could not read SBML'

	#Localisations
	print 'Get compartments'
	try:
		compartments = []
		for comp in model.getListOfCompartments():
			compartment = {}
			compartment["id"] = "organelle_" + comp.getId()
			compartment["name"] = comp.getName()
			compartment["notes"] = getNotes(comp)
			compartment["SBO"] = comp.getSBOTermID()
			compartment["type"] = "organelle"
			compartment["source"] = "SBML"
			compartment["sourceId"] = comp.getId()
			compartment["synonyms"] = []
			compartment["inCompartment"] = comp.getOutside()
			compartment["property"] = []
			compartment["tags"] = []
			compartment["sandbox"] = []
			compartment["references"] = []
			compartment["notifications"] = []
			annotations = getAnnotations(comp)
			for ann in annotations:
				compartment[ann] = annotations[ann] 
			compartments.append(compartment)
		print '...complete'
	except:
		compartments = []
		print 'Could not read compartments'

	#Molecules
	print 'Get molecules...'
	molecules = []
	for mol in model.getListOfSpecies():
		molecule = {}
		molecule["id"] = mol.getId()
		molecule["name"] = mol.getName()
		molecule["SBO"] = mol.getSBOTermID()
		if mol.getSBOTermID() == 'SBO:0000247':
			molecule["tags"] = ['simple chemical']
		elif mol.getSBOTermID() == 'SBO:0000245':
			molecule["tags"] = ['macromolecule']
		else:
			molecule["tags"] = []
		molecule["notes"] = getNotes(mol)
		molecule["type"] = "molecule"
		molecule["source"] = "SBML"
		molecule["sourceId"] = mol.getId()
		molecule["synonyms"] = []
		molecule["inCompartment"] = "organelle_" + mol.getCompartment()
		molecule["property"] = []
		molecule["sandbox"] = []
		molecule["references"] = []
		molecule["notifications"] = []
		annotations = getAnnotations(mol)
		for ann in annotations:
			molecule[ann] = annotations[ann] 
		molecules.append(molecule)
	print '...complete'

	#Reactions
	print 'Get reactions...'
	try:	
		reactions = []
		for rxn in model.getListOfReactions():
			reaction = {}
			reaction["id"] = "organelle_" + rxn.getId()
			reaction["name"] = rxn.getName()
			reaction["SBO"] = rxn.getSBOTermID()
			reaction["tags"] = []
			reaction["notes"] = getNotes(rxn)
			reaction["type"] = "reaction"
			reaction["source"] = "SBML"
			reaction["sourceId"] = rxn.getId()
			reaction["synonyms"] = []
			reaction["property"] = []
			reaction["sandbox"] = []
			reaction["references"] = []
			reaction["notifications"] = []
			annotations = getAnnotations(rxn)
			for ann in annotations:
				reaction[ann] = annotations[ann] 

			listOfReactants = rxn.getListOfReactants()
			reactants = []
			for mol in listOfReactants:
				name,comp = fetch_species(model,mol.getSpecies())
				reactants.append({'id':mol.getSpecies(), 'stoichiometry':mol.getStoichiometry(), 'display':name, 'localisation':comp})
			reaction["reactants"] = reactants

			listOfModifiers = rxn.getListOfModifiers()
			modifiers = []
			for mol in listOfModifiers:
				name,comp = fetch_species(model,mol.getSpecies())
				modifiers.append({'id':mol.getSpecies(), 'display':name, 'localisation':comp})
			reaction["modifiers"] = modifiers

			listOfProducts = rxn.getListOfProducts()
			products = []
			for mol in listOfProducts:
				name,comp = fetch_species(model,mol.getSpecies())
				products.append({'id':mol.getSpecies(), 'stoichiometry':mol.getStoichiometry(), 'display':name, 'localisation':comp})
			reaction["products"] = products

		reactions.append(reaction)
		print '...complete'
	except Exception, e:
		reactions = []
		print 'Could not read reactions:', str(e)
	
	model = {"compartments":compartments, "molecules":molecules, "reactions":reactions }
	
	return model

'''
f = open('yeast_7.6_recon.xml', 'r')
sbml = f.read()
recon = parseSBML(sbml)
'''