from libsbml import *
import xmltodict, json, copy, uuid, requests

NEO_USERNAME = 'neo4j'
NEO_PASSWORD = 'donkey'

#Remapping required
nodeBlank = {"id":"", "type":"", "source":"", "sourceId":"", "name":"",	"synonyms":[], "inCompartment":"", "is":[], "is_a":[], "isDescribedBy":[], "isVersionOf":[], "property":[], "tags":[], "sandbox":[], "references":[], "notifications":[]	}

#Empty annotations
annotationsBlank = {'BQB_IS':[],'BQB_HAS_PART':[],'BQB_IS_PART_OF':[],'BQB_IS_VERSION_OF':[],'BQB_HAS_VERSION':[],'BQB_IS_HOMOLOG_TO':[],'BQB_IS_DESCRIBED_BY':[],'BQB_IS_ENCODED_BY':[],'BQB_ENCODES':[],'BQB_OCCURS_IN':[],'BQB_HAS_PROPERTY':[],'BQB_IS_PROPERTY_OF':[],'BQB_HAS_TAXON':[]}

#For qualifier index resolution
qualifiers = ['BQB_IS','BQB_HAS_PART','BQB_IS_PART_OF','BQB_IS_VERSION_OF','BQB_HAS_VERSION','BQB_IS_HOMOLOG_TO','BQB_IS_DESCRIBED_BY','BQB_IS_ENCODED_BY','BQB_ENCODES','BQB_OCCURS_IN','BQB_HAS_PROPERTY','BQB_IS_PROPERTY_OF','BQB_HAS_TAXON']

#Change SBO names to EMPATHY fields
def annotationMapper(properties,source):
	properties["is"] = source["BQB_IS"]
	properties["hasPart"] = source["BQB_HAS_PART"]
	properties["isPartOf"] = source["BQB_IS_PART_OF"]
	properties["isVersionOf"] = source["BQB_IS_VERSION_OF"]
	properties["hasVersion"] = source["BQB_HAS_VERSION"]
	properties["isHomologTo"] = source["BQB_IS_HOMOLOG_TO"]
	properties["isDescribedBy"] = source["BQB_IS_DESCRIBED_BY"]
	properties["isEncodedBy"] = source['BQB_IS_ENCODED_BY']
	properties["encodes"] = source['BQB_ENCODES']
	properties["occursIn"] = source['BQB_OCCURS_IN']
	properties["hasProperty"] = source['BQB_HAS_PROPERTY']
	properties["isPropertyOf"] = source['BQB_IS_PROPERTY_OF']
	properties["hasTaxon"] = source['BQB_HAS_TAXON']
	return properties

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
	annotations = copy.deepcopy(annotationsBlank)
	id = sbmlObject.getId()
	try:
		try:
			cvTerms = sbmlObject.getCVTerms()
			if cvTerms:
				for cvTerm in cvTerms:
					qual = qualifiers[cvTerm.getBiologicalQualifierType()]
					attributes = cvTerm.getResources()
					for i in range(0,attributes.getLength()):
						previous = annotations[qual]
						uri = attributes.getValue(i)
						uri = uri.replace("http://identifiers.org/","")
						key,val = uri.split('/')
						previous.append( json.dumps([key,val]) )
						annotations[qual] = previous
					return annotations
			else:
				return annotations
		except Exception, e:
			print 'Annotation exception', str(e), sbmlObject.getId()
			return annotations
	except Exception, e:
		print 'Exception:', str(e)
		return False

def getNotes(sbmlObject):
	try:
		notesString = sbmlObject.getNotesString()
		notesString = notesString.replace("\n","")
		notesString = notesString.replace("  ","")
		notesDict = xmltodict.parse( notesString )
		notesDict = json.loads( json.dumps( notesDict ) )
		'''
		notes = []
		for html in notesDict['notes']['html:body']['html:p']:
			notes.append( str(html) )
		return notes
		'''
		return json.dumps(notesDict)
	except:
		return ""

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
			compartment["type"] = "compartment"
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
			for qual in annotations.keys():		
				compartment[qual] = annotations[qual]
			compartments.append(compartment)
		print '...complete'
	except Exception, e:
		compartments = []
		print 'Could not read compartments', str(e)

	#Molecules
	print 'Get molecules...'
	try:
		molecules = []
		for mol in model.getListOfSpecies():
			molecule = {}
			molecule["id"] = mol.getId()
			molecule["name"] = mol.getName()
			molecule["SBO"] = mol.getSBOTermID()
			if mol.getSBOTermID() == 'SBO:0000247':
				molecule["type"] = 'simple chemical'
				molecule["tags"] = ['simple chemical']
			elif mol.getSBOTermID() == 'SBO:0000245':
				molecule["type"] = 'macromolecule'
				molecule["tags"] = ['macromolecule']
			elif mol.getSBOTermID() == 'SBO:0000253':
				molecule["tags"] = ['complex']				
				molecule["type"] = 'complex'
			else:
				molecule["type"] = 'None'
				molecule["tags"] = []
			molecule["charge"] = mol.getCharge()
			molecule["notes"] = getNotes(mol)
			molecule["source"] = "SBML"
			molecule["sourceId"] = mol.getId()
			molecule["synonyms"] = []
			molecule["inCompartment"] = "organelle_" + mol.getCompartment()
			molecule["property"] = []
			molecule["sandbox"] = []
			molecule["references"] = []
			molecule["notifications"] = []
			annotations = getAnnotations(mol)
			for qual in annotations.keys():		
				molecule[qual] = annotations[qual]
			molecules.append(molecule)
		print '...complete'
	except Exception, e:
		molecules = []
		print 'Could not read molecules', str(e)

	#Reactions
	print 'Get reactions...'
	try:	
		reactions = []
		for rxn in model.getListOfReactions():
			reaction = {}
			reaction["id"] = rxn.getId()
			reaction["name"] = rxn.getName()
			reaction["SBO"] = rxn.getSBOTermID()
			reaction["tags"] = []
			reaction["notes"] = getNotes(rxn)
			reaction["type"] = "reaction"
			reaction["source"] = "SBML"
			reaction["sourceId"] = rxn.getId()
			reaction["synonyms"] = []
			reaction["property"] = []
			reaction["references"] = []
			if rxn.getReversible():
				reaction["reversible"] = "true"
			else:
				reaction["reversible"] = "false"
			annotations = getAnnotations(rxn)
			for qual in annotations.keys():		
				reaction[qual] = annotations[qual]

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
		print '...complete', len(reactions), len( model.getListOfReactions() )
	except Exception, e:
		reactions = []
		print 'Could not read reactions:', str(e)
	
	model = {"compartments":compartments, "molecules":molecules, "reactions":reactions }

	return model

def collectCyphers(model):
	cyphers = []
	#Compartments
	print 'Collect compartments'
	compartments = model['compartments']
	for comp in compartments:
		#print comp
		properties = {}
		properties["id"] = str(uuid.uuid4())
		properties["type"] = "compartment"
		properties["name"] = comp['name']
		properties["sourceId"] = comp['id']	
		properties["source"] = "SBML"
		properties["synonyms"] = []
		properties["tags"] = ["compartment"]
		properties["inCompartment"] = comp['inCompartment']
		#Map SBO to EMPATHY
		properties = annotationMapper(properties,comp)	
		cyphers.append(['CREATE (n:compartment {props}) RETURN n.id', properties, 'compartment'])

	#Molecules
	print 'Collect molecules'
	mols = model['molecules']
	for mol in mols:
		#print mol
		properties = {}
		properties["id"] = str(uuid.uuid4())
		properties["type"] = "molecule"
		properties["name"] = mol['name']
		properties["sourceId"] = mol['id']	
		properties["source"] = "SBML"
		properties["synonyms"] = []
		properties["tags"] = ["molecule"]
		properties["inCompartment"] = mol['inCompartment']
		#Map SBO to EMPATHY
		properties = annotationMapper(properties,mol)	
		cyphers.append(['CREATE (n:molecule {props}) RETURN n.id', properties, 'molecule'])

	#Reactions
	print 'Collect reactions'
	rxns = model['reactions']
	for rxn in rxns:
		#print rxn
		properties = {}
		properties["id"] = str(uuid.uuid4())
		properties["type"] = "reaction"
		properties["name"] = rxn['name']
		properties["sourceId"] = rxn['id']	
		properties["source"] = "SBML"
		properties["synonyms"] = []
		properties["tags"] = ["reaction"]
		#Map SBO to EMPATHY
		properties = annotationMapper(properties,rxn)	
		cyphers.append(['CREATE (n:reaction {props}) RETURN n.id', properties, 'reaction'])

	#Reactions
	print 'Collect relations'
	rxns = model['reactions']
	for rxn in rxns:
		rxnId = rxn['id']

		if 'reactants' in rxn:
			for reactant in rxn['reactants']:
				molId = reactant['id']
				properties = {}
				stoichiometry = reactant['stoichiometry']
				cypher = 'MATCH (r:reaction {id:"' + rxnId + '"}), (m:molecule {sourceId:"' + molId + '"}) CREATE (r)-[s:hasReactant]->(m) SET s.stoichiometry="' + str(stoichiometry) + '" RETURN r.id'
				cyphers.append([cypher, properties, 'reactant'])

		if 'modifiers' in rxn:
			for modifier in rxn['modifiers']:
				molId = modifier['id']
				properties = {}
				cypher = 'MATCH (r:reaction {sourceId:"' + rxnId + '"}), (m:molecule {sourceId:"' + molId + '"}) CREATE (r)-[s:hasModifier]->(m) RETURN r.id'
				cyphers.append([cypher, properties, 'modifier'])

		if 'products' in rxn:
			for product in rxn['products']:
				molId = product['id']
				properties = {}
				stoichiometry = product['stoichiometry']
				cypher = 'MATCH (r:reaction {sourceId:"' + rxnId + '"}), (m:molecule {sourceId:"' + molId + '"}) CREATE (r)-[s:hasProduct]->(m) SET s.stoichiometry="' + str(stoichiometry) + '"  RETURN r.id'
				cyphers.append([cypher, properties, 'product'])

	return cyphers

#Master
def sbml2cyphers(sbml):
	model = parseSBML(sbml)
	cyphers = collectCyphers(model)
	return cyphers

'''
f = open('yeast_7.6_recon.xml', 'r')
sbml = f.read()
print uploadSBML(sbml,'7474')
'''