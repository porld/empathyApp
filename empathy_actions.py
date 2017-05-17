from rdkit import Chem
import requests, json, time, xmltodict, libchebipy

#-----------------------------------------------------------------------------------------
#Chemical converters
def pickInchi(inchis):
	if 'CHEBI' in inchis:
		return inchis['CHEBI'], 'CHEBI'
	elif 'KEGG' in inchis:
		return inchis['KEGG'], 'KEGG'
	elif 'PUBCHEM' in inchis:
		return inchis['PUBCHEM'], 'PUBCHEM'
	elif 'MetaCyc' in inchis:
		return inchis['MetaCyc'], 'MetaCyc'
	else:
		return '', ''

#Mol string to RDKit mol
def molString2mol(molString):
	mol = Chem.MolFromMolBlock(molString)
	if mol is None:
		print 'Mol block failure'
		return ''
	else:
		return mol

#InChI to smiles		
def inchi2smiles(inchi):
	try:
		mol = Chem.MolFromInchi(inchi)
		smiles = Chem.MolToSmiles(mol)
		return smiles
	except Exception, e:
		print "Error", str(e)
		return False
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
#PROTEIN NAME TO IDENTITY
#Keyword, organism to PDB
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
				matches.append({"pdb":'?'+pdb})
	return matches[0:2]

#Keyword, organism to Uniprot
def proteinIdentifier(keyword,organism):
	url = 'http://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:' + str(organism) + '+AND+' + keyword + '&format=tab&columns=id&limit=10'
	content = requests.get(url)
	content = content.text
	lines = content.split('\n')
	
	matches = []
	for line in lines[1:(len(lines)-1)]:
		cols = line.split('\t')
		matches.append({"uniprot":'?'+cols[0]})
	return matches
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
##CHEMICAL NAME TO IDENTITY
#Keyword to smiles
def cactvs_keyword2smiles(keyword):
	url = 'https://cactus.nci.nih.gov/chemical/structure/' + keyword + '/smiles'
	content = requests.get(url)
	if('Page not found (404)' in content.text):
		return ''
	return content.text

#Keyword to InChI
def cactvs_keyword2inchi(keyword):
	url = 'https://cactus.nci.nih.gov/chemical/structure/' + keyword + '/stdinchi'
	content = requests.get(url)
	if('Page not found (404)' in content.text):
		return ''
	return content.text
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
#CHEMICAL SYNONYMS
def smallMoleculeSynonyms(keyword):
	url = 'https://cactus.nci.nih.gov/chemical/structure/' + keyword + '/names'
	content = requests.get(url)
	if('Page not found (404)' in content.text):
		return []
	content = content.text
	lines = content.split('\n')
	possibles = []
	for line in lines[1:(len(lines)-1)]:
		possibles.append('?'+line)
	return possibles[0:9]
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
#CHEMICAL ID TO IDENTITY
#MetaCyc to smiles
def metacyc2smiles(metacyc):
	try:
		url = 'https://metacyc.org/getxml?id=META:'+metacyc
		content = requests.get(url)
		answer = content.text
		metadict = xmltodict.parse(answer)
		time.sleep(2)
		inchi = str( metadict['ptools-xml']['Compound']['inchi']['#text'] )
		return inchi2smiles(inchi)
	except:
		return False

#ChEBI to smiles
def chebi2smiles(chebi):
	chebi_entity = libchebipy.ChebiEntity(chebi)
	inchi = chebi_entity.get_inchi()
	if inchi is None:
		print 'No InChI:', inchi
		return False
	else:
		print 'Found InChI:', inchi
		return inchi2smiles(inchi)

#KEGG to smiles
def kegg2smiles(keggId):
	try:
		url = 'http://www.genome.jp/dbget-bin/www_bget?-f+m+compound+'+keggId
		content = requests.get(url)
		mol = molString2mol(content.text)
		smiles = Chem.MolToSmiles(mol)
		return smiles
	except:
		return False

#PubChem to smiles
def pubchem2smiles(pubchem):
	try:
		url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'+pubchem+'/record/SDF/?record_type=3d&response_type=display'
		content = requests.get(url)
		mol = molString2mol(content.text)
		smiles = Chem.MolToSmiles(mol)
		return smiles
	except:
		return False
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
#Jamboree overseers
def chemical2structure(keyword,isList):
	print 'chemical2structure', keyword, isList, type(isList)
	message = 'Search for structure: none found'
	for entry in isList:
		known = entry[0]
		print known
		#Chase up ChEBI
		if 'chebi' in known:
			print "Search ChEBI", entry[1]
			smiles = chebi2smiles(entry[1])
			if smiles:
				message = 'Retrieved structure from ChEBI'
				return json.dumps(["smiles",smiles]), message
		elif known is 'kegg':
			smiles = kegg2smiles(entry[1])
			if smiles:
				message = 'Retrieved structure from KEGG'
				return json.dumps(["smiles",smiles]), message
		elif known is 'metacyc':
			smiles = metacyc2smiles(entry[1])
			if smiles:
				message = 'Retrieved structure from MetaCyc'
				return json.dumps(["smiles",smiles]), message
		elif known is 'pubchem':
			smiles = pubchem2smiles(entry[1])
			if smiles:
				message = 'Retrieved structure from PubChem'
				return json.dumps(["smiles",smiles]), message
		else:
			pass
	#Last resort is CACTUS
	if len(matches) == 0:
		smiles = cactvs_keyword2smiles(keyword)
		message = 'Retrieved structure by "' + keyword + '" from CACTUS'
		return json.dumps(["smiles",'?'+smiles]), message
	else:
		return False, 'None found'
#-----------------------------------------------------------------------------------------

'''
print chemical2structure('glucose',[["pubchem","5793"]])
'''