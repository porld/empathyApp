import requests

def pdbIdentifier(keyword,organism):
	url = 'http://www.uniprot.org/uniprot/?query=organism:' + str(organism) + '+AND+database:pdb+AND+' + keyword + '&format=tab&columns=database(PDB)'
	content = requests.get(url)
	content = content.text
	lines = content.split('\n')
	
	matches = []
	for line in lines[1:(len(lines)-1)]:
		cols = line.split(';')
		for pdb in cols:
			if len(pdb) == 4:
				matches.append({"pdb":pdb+'?'})
	return matches[0:2]

def proteinIdentifier(keyword,organism):
	url = 'http://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:' + str(organism) + '+AND+' + keyword + '&format=tab&columns=id&limit=10'
	content = requests.get(url)
	content = content.text
	lines = content.split('\n')
	
	matches = []
	for line in lines[1:(len(lines)-1)]:
		cols = line.split('\t')
		matches.append({"uniprot":cols[0]+'?'})
	return matches

def cactvs_keyword2smiles(keyword):
	url = 'https://cactus.nci.nih.gov/chemical/structure/' + keyword + '/smiles'
	content = requests.get(url)
	if('Page not found (404)' in content.text):
		return ''
	return content.text

def cactvs_keyword2inchi(keyword):
	url = 'https://cactus.nci.nih.gov/chemical/structure/' + keyword + '/stdinchi'
	content = requests.get(url)
	if('Page not found (404)' in content.text):
		return ''
	return content.text

def smallMoleculeSynonyms(keyword):
	url = 'https://cactus.nci.nih.gov/chemical/structure/' + keyword + '/names'
	content = requests.get(url)
	if('Page not found (404)' in content.text):
		return []
	content = content.text
	lines = content.split('\n')
	possibles = []
	for line in lines[1:(len(lines)-1)]:
		possibles.append(line+'?')
	return possibles[0:9]

def smallMoleculeIdentifier(keyword):
	matches = []
	matches.append(["smiles",cactvs_keyword2smiles(keyword)+'?'])
	matches.append(["InChI",cactvs_keyword2inchi(keyword)+'?'])
	return matches