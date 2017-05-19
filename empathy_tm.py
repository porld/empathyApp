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
	results = []
	for hit in hits:
		for highlight in hit["highlight"]:
			processed = {}
			processed['id'] = hit['_id']
			processed['pmid'] = hit["_source"]["MedlineCitation"]["pmid"]["content"]
			processed['highlight'] = json.dumps(highlight)
			results.append(json.dumps (processed) ) #Stored in database as list of JSON 
	return results

#Calliope query
def queryCalliope(triples):
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
	triples = []

	try:
		#Get organism
		triple = buildTriple('Species',organism)
		triples.append(triple)

		#Get reactants
		for mol in record['listOfReactants']:
			molNames = [mol['name']]
			for syn in mol['synonyms']:
				if '?' in syn:
					syn = syn.replace('?','')
				molNames.append(syn)
			triple = buildTriple('Chemical',molNames)
			triples.append(triple)

		#Get products
		for mol in record['listOfProducts']:
			molNames = [mol['name']]
			for syn in mol['synonyms']:
				if '?' in syn:
					syn = syn.replace('?','')
				molNames.append(syn)
			triple = buildTriple('Chemical',molNames)		
			triples.append(triple)

		return triples
	except Exception, e:
		print str(e)
		return False

#Pipeline
def calliopeCoordinator(rxn,org):
	triples = rxn2triples(rxn,org)
	if triples:
		hits = queryCalliope(triples)
		return processCalliope(hits)
	else:
		print 'Malformed query'
		return False

#Test
rxn = {}
rxn['listOfReactants'] = [{'name':'glucose','synonyms':[]},{'name':'ATP','synonyms':['adenosine triphosphate']}]
rxn['listOfProducts'] = [{'name':'glucose-6-phosphate','synonyms':['G6P']},{'name':'ADP','synonyms':[]}]
for res in calliopeCoordinator(rxn,['yeast','Saccharomyces','S. cerevisiae']):
	print res
