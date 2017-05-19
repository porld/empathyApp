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
	   "number_of_fragments":5,
		"fragment_size":150,
		"tag_schema":"styled",
		"fields":{"*":{"pre_tags":["<b>"],"post_tags":["</b>"]}}
	},
	"from": 0,
	"size": 10
	}
	print query
	try:
		url = 'http://nactem10.mib.man.ac.uk:5004/triplesAPI'
		content = requests.post(url,json=query )
		response = content.json()
		print 'Calliope server says', response
		hits = response['hits']['hits']
		return hits, 'found ' + len(hits) + '  hits'
	except Exception, e:
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
			molNames = [ mol['name'] ]
			triple = buildTriple('Chemical',molNames)		
			triples.append(triple)

		#Get products
		for mol in record['listOfProducts']:
			molNames = [ mol['name'] ]
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
		hits, calliope_msg = queryCalliope(triples)
		return processCalliope(hits), calliope_msg
	else:
		print 'Malformed query'
		return False, 'query failure'


#Test
rxn = {}
rxn['listOfReactants'] = [{'name':'glucose'},{'name':'ATP'}]
rxn['listOfProducts'] = [{'name':'glucose 6-phosphate'},{'name':'ADP'}]
print calliopeCoordinator(rxn,['yeast','Saccharomyces','S. cerevisiae'])
