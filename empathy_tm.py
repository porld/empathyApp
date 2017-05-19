import requests, json

#Calliope triple
def buildTriple(type,word_form):
	triple = {}
	triple['type'] = type
	triple['id'] = []
	triple["word_form"] = word_form
	return triple

#Calliope response procesing
def processCalliope(hits):
	print 'Processing Calliope results...'
	results = []
	for hit in hits:
		for highlight in hit["highlight"]:
			processed = {}
			processed['id'] = hit['_id']
			processed['pmid'] = hit["_source"]["MedlineCitation"]["pmid"]["content"]
			processed['highlight'] = json.dumps(highlight)
			results.append(json.dumps (processed) ) #Stored in database as list of JSON 
	print 'Calliope hits', len(hits)
	return results

#Calliope query
def queryCalliope(triples):
	print 'Querying Calliope...'
	query = {
	"triples": triples,
	"highlight": {
	   "number_of_fragments":3,
		"fragment_size":100,
		"tag_schema":"styled",
		"fields":{"*":{"pre_tags":["<b>"],"post_tags":["</b>"]}}
	},
	"from": 0,
	"size": 10
	}
	try:
		url = 'http://nactem10.mib.man.ac.uk:5004/triplesAPI'
		content = requests.post(url,json=query )
		response = content.json()
		return response['hits']['hits']
	except Exception, e:
		return str(e)

#Reaction to Calliope format (json, list of strings)
def rxn2triples(record,organism):
	print 'Building triples...'
	triples = []

	try:
		#Get organism
		triple = buildTriple('Species',organism)
		triples.append(triple)

		#Gather chemical names
		molNames = []

		#Get reactants
		for mol in record['listOfReactants']:
			molNames.append( mol['name'] )

		#Get products
		for mol in record['listOfProducts']:
			molNames.append( mol['name'] )

		#Build chemical triple
		triple = buildTriple('Chemical',molNames)		
		triples.append(triple)

		return triples
	except Exception, e:
		print str(e)
		return False

#Pipeline
def calliopeCoordinator(rxn,org):
	print 'Hitting Calliope...'
	triples = rxn2triples(rxn,org)
	if triples:
		hits = queryCalliope(triples)
		return processCalliope(hits)
	else:
		print 'Malformed query'
		return False

'''
#Test
rxn = {}
rxn['listOfReactants'] = [{'name':'D-glucose','synonyms':[]},{'name':'ATP','synonyms':[]}]
rxn['listOfProducts'] = [{'name':'D-glucose 6-phosphate','synonyms':[]},{'name':'ADP','synonyms':[]}]
print calliopeCoordinator(rxn,['yeast','Saccharomyces','S. cerevisiae'])
'''