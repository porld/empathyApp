import requests, json, re

#Some very common names (like water) confound the cooccurrence metric, so filter them out (BODGE!)
def naughtyCheck(name):
	naughty = ['water', 'proton', 'H+','H2O']
	if name not in naughty:
		return True
	else:
		return False

#Calliope triple
def buildTriple(type,word_form):
	triple = {}
	triple['type'] = type
	triple['id'] = []
	triple["word_form"] = word_form
	return triple

#Calliope response procesing
def processCalliope(hits):
	print 'Processing Calliope results...', hits
	results = []
	for hit in hits:
		for highlight in hit["highlight"]:
			processed = {}
			processed['id'] = hit['_id']
			processed['pmid'] = hit["_source"]["MedlineCitation"]["pmid"]["content"]
			processed['highlight'] = json.dumps(highlight)
			processed['curated'] = 'uncurated'
			results.append(json.dumps (processed) ) #Stored in database as list of JSON 
	print 'Calliope hits', len(results)
	return results

#Calliope query
def queryCalliope(triples):
	print 'Querying Calliope...'
	query = {
	"triples": triples,
	"highlight": {
	   "number_of_fragments":10,
		"fragment_size":200,
		"tag_schema":"styled",
		"fields":{"*":{"pre_tags":[""],"post_tags":[""]}}
	},
	"from": 0,
	"size": 10
	}
	print query
	try:
		url = 'http://nactem10.mib.man.ac.uk:5004/triplesAPI'
		content = requests.post(url,json=query )
		response = content.json()
		hits = response['hits']['hits']
		return hits, 'found ' + str(len(hits)) + '  hits'
	except Exception, e:
		print 'queryCalliope error', str(e)
		return str(e), 'query error'

#Reaction to Calliope format (json, list of strings)
def rxn2triples(record,organism):
	print 'Building triples...'
	triples = []

	try:
		#Get organism
		triple = buildTriple('Species',organism)
		triples.append(triple)

		#Get reactants
		for mol in record['listOfReactants']:
			name = mol['name']
			name = re.sub(r"^\s+", "", name)
			name = re.sub(r"\s+$", "", name)
			if naughtyCheck(name):
				triple = buildTriple('Chemical', [name])		
				triples.append(triple)

		#Get products
		for mol in record['listOfProducts']:
			name = mol['name']
			name = re.sub(r"^\s+", "", name)
			name = re.sub(r"\s+$", "", name)
			if naughtyCheck(name):
				triple = buildTriple('Chemical', [name])		
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
		hits, calliope_msg = queryCalliope(triples)
		return processCalliope(hits), calliope_msg
	else:
		print 'Malformed query'
		return False, 'query failure'

'''
#Test
rxn = {}
rxn['listOfReactants'] = [{'name':'glucose'},{'name':'ATP'}]
rxn['listOfProducts'] = [{'name':'glucose 6-phosphate'},{'name':'ADP'},{'name':'H+ '}]
response, hitCount = calliopeCoordinator(rxn,['yeast','Saccharomyces','S. cerevisiae'])
print response
'''