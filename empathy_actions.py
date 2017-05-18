from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.rdmolops import GetFormalCharge
import requests, json, time, xmltodict, libchebipy

#-----------------------------------------------------------------------------------------
#Chemical converters

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
	try:
		url = 'https://cactus.nci.nih.gov/chemical/structure/' + keyword + '/names'
		content = requests.get(url)
		if('Page not found (404)' in content.text):
			return False, 'no synonym found'
		content = content.text
		lines = content.split('\n')
		possibles = []
		for line in lines[1:(len(lines)-1)]:
			possibles.append('?'+line)
		return possibles[0:9], 'got synonyms from CACTUS'
	except:
		return False, 'synonym search error'
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
#CHEMICAL PROPERTIES
def syncProperties(smiles):
	try:
		mol = Chem.MolFromSmiles(smiles)
		formula = CalcMolFormula(mol)
		charge = GetFormalCharge(mol)
		return formula, charge, 'calculated properties from structure'
	except:
		return False, False, 'property calculation error'
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

#Coordinator
def chemical2structure(keyword,isList):
	print 'ACTION chemical2structure', keyword, isList, type(isList)
	message = 'Search for structure: none found'
	for entry in isList:
		entry[0] = str(entry[0])
		entry[1] = str(entry[1])
		known = entry[0]
		#Chase up ChEBI
		if 'smiles' in known:
			print "ACTION Already got SMILES", entry[1]
			message = 'already had structure'
			return False, message
		elif 'InChI' in known:
			print "ACTION Convert InChI", entry[1]
			inchi = entry[1]
			if '?' in inchi:
				inchi = inchi.replace('?','')
			smiles = inchi2smiles(inchi)
			if smiles:
				message = 'built structure from InChI'
				return json.dumps(["smiles",smiles]), message
		elif 'chebi' in known:
			print "ACTION Search ChEBI", entry[1]
			smiles = chebi2smiles(entry[1])
			if smiles:
				message = 'got structure from ChEBI'
				return json.dumps(["smiles",smiles]), message
		elif 'kegg.compound' in known:
			print "ACTION Search KEGG", entry[1]
			smiles = kegg2smiles(entry[1])
			if smiles:
				message = 'got structure from KEGG'
				return json.dumps(["smiles",smiles]), message
		elif 'metacyc' in known:
			print "ACTION Search MetaCyc", entry[1]
			smiles = metacyc2smiles(entry[1])
			if smiles:
				message = 'got structure from MetaCyc'
				return ["smiles",smiles], message
		elif 'pubchem' in known:
			print "ACTION Search Pubchem", entry[1]
			smiles = pubchem2smiles(entry[1])
			if smiles:
				message = 'got structure from PubChem'
				return json.dumps(["smiles",smiles]), message
		else:
			print "ACTION Cannot use", entry
	#Last resort is CACTUS
	smiles = cactvs_keyword2smiles(keyword)
	if smiles:
		print "ACTION Search CACTUS by name", keyword
		message = 'got structure from CACTUS'
		return json.dumps(["smiles",'?'+smiles]), message
	else:
		print "ACTION Could not find", keyword
		return False, 'no structure found'
#-----------------------------------------------------------------------------------------

#Tests
#print chemical2structure("glucose", [ [u'chebi', u'CHEBI:37447']])
#print syncProperties("C(C(C(=O)[O-])N)C(=O)[O-]")