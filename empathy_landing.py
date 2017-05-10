from flask import Flask, make_response, url_for, request
from flask_restful import Resource
from flask_httpauth import HTTPBasicAuth
from passlib.hash import sha256_crypt
from flask_mail import Mail, Message
import sys, subprocess, uuid, os, json, requests, socket, time, pickle, passwordmeter, copy, datetime
import empathy_actions as actions
import empathy_libsbml as libsbml

HOST = 'www.metabolicjamboree.co.uk'
PORT = 8080


#Initialise app
app = Flask(__name__, static_url_path='')
app.config['SECRET_KEY'] = 'big_secret'

if len(sys.argv) != 2:
	print 'Missing arguments (landing.py <email password> <port>)'

mail_username = 'empathy.manchester@gmail.com'
mail_password = sys.argv[1]

#-----------------------------------------------------------------------------------------
#Mail config
app.config.update(
	DEBUG=True,
	MAIL_SERVER='smtp.gmail.com',
	MAIL_PORT=465,
	MAIL_USE_SSL=True,
	MAIL_USERNAME = mail_username,
	MAIL_PASSWORD = mail_password
	)
mail = Mail(app)
#Collaborator tokens
tokens = {}
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
#Neo4J Credentialas
NEO_USERNAME = 'neo4j'
NEO_PASSWORD = 'empathy_donkey'
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
#Authorisation
auth = HTTPBasicAuth()

def relaunch_dockers(user_credentials):
	#First collect a list of unique Docker names
	dockers = []
	for user in user_credentials:
		#print 'User:', user
		reconstructions = user_credentials[user]['reconstructions']
		for recon in reconstructions:
			print recon['fsid']
			dockers.append(recon['fsid'])
	#Unique
	dockers = list(set(dockers))
	successes = []
	failures = []
	for d in dockers:
		try:
			subprocess.Popen(["docker", "start", d ])
			print d, 'success'
			successes.append(d)
		except:
			print d, 'failed'
			failures.append(d)
	return successes, failures

#Load previous credentials
try:
	user_credentials = pickle.load( open( "users.pickle", "rb" ) )
	print user_credentials
	print 'Resurrecting previous Dockers...'
	successes, failures = relaunch_dockers(user_credentials)
	print 'Relaunched Dockers:', successes
	print 'Failed to relaunch:', failures
except:
	user_credentials  = {}
	pickle.dump( user_credentials, open( "users.pickle", "wb" ) )

#What's this for?
#user_credentials = {}

@auth.verify_password
def verify_password(username, password):
    if username in user_credentials:
    	user = user_credentials[username]
        return sha256_crypt.verify(password, user["password"])
    return False

@auth.get_password
def get_password(username):
    if username in user_credentials:
    	user = user_credentials[username]
    	password = user['password']
    	return password
    return None

@auth.error_handler
def unauthorized():
    return make_response(json.dumps({'error': 'Unauthorized access'}), 401)
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
def flatten_parameters(parameters):
	temp = {}
	for key in parameters:
		value = parameters[key]
		if type(value) is list:
			value = list(set(value))
		temp[key] = value
	return temp

#Neo4J functions
def send_cypher(cypher,parameters,port):
	headers = {'content-type': 'application/json'}
	url = 'http://localhost:' + port + '/db/data/transaction/commit'
	parameters = flatten_parameters(parameters)
	query = { "statements" : [ { "statement" : cypher, "parameters" : parameters } ] }
	##print url, NEO_USERNAME, NEO_PASSWORD
	return requests.post(url, auth=(NEO_USERNAME,NEO_PASSWORD), data=json.dumps(query), headers=headers)

#Socket  function
def send_message(handle,message,port):
	headers = {'content-type': 'application/json'}
	url = 'https://' + HOST + ':8082/socket'
	try:
		#print 'Hitting socket:', url
		requests.post(url, data=json.dumps({"handle":handle,"message":message,"port":port}), headers=headers)
		return True
	except:
		#print 'Could not hit socket:', url
		return False
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
##USER HANDLING
#Register a new user
@app.route('/registerUser', methods=['POST'])
def register_user():
	json_data = request.get_json(force=True)
	username = json_data['username']
	password = json_data['password']
	#print username, password
	#Add to in-memory list for checking
	#print 'Check username is novel'
	if username in user_credentials:
		return json.dumps([False, 'Select different username'])
	#print 'Check username is long'
	if len(username) < 8:
		return json.dumps([False, 'Username too short'])
	#print 'Check password strength'
	strength, _ = passwordmeter.test(password)
	if strength < 0.3:
		return json.dumps([False, 'Password too weak'])
	#print 'Encrypt password'
	hash_pwd = sha256_crypt.encrypt(password)
	#print 'Create user'
	user_credentials[username] = {	"password": hash_pwd, "reconstructions": []	}
	#Update the pickle
	#print 'Update pickle'
	pickle.dump( user_credentials, open( "users.pickle", "wb" ) )
	return json.dumps([True, username])

@app.route('/signIn', methods=['POST'])
def sign_in():
	json_data = request.get_json(force=True)
	username = json_data['username']
	password = json_data['password']
	try:
		user = user_credentials[username]
		if sha256_crypt.verify(password, user["password"]):
			return json.dumps([True, 'Successfully signed in'])
		else:
			return json.dumps([False, 'Could not sign in'])
	except:
		return json.dumps([False, 'Could not sign in'])

'''
@app.route('/changePwd', methods=['POST'])
def change_pwd():
	json_data = request.get_json(force=True)
	username = json_data['username']
	password = json_data['password']
	new_password = json_data['new_password']
	#Verify current user
	if username not in user_credentials:
		return json.dumps([False, 'Unknown user'])
	else:
		user = user_credentials[username]
	if sha256_crypt.verify(password, user["password"]):
		strength, _ = passwordmeter.test(new_password)
		if strength < 0.3:
			return json.dumps([False, 'New password too weak'])
		reconstructions = user['reconstructions']
		hash_pwd = sha256_crypt.encrypt(new_password)
		user_credentials[username] = {	"password": hash_pwd, "reconstructions": reconstructions	}
		#Update the pickle
		pickle.dump( user_credentials, open( "users.pickle", "wb" ) )
		return json.dumps([True, username])
	else:
		return json.dumps([False, 'Could not validate user'])
'''

@app.route('/pwdStrength/<pwd>', methods=['GET'])
def pwdStrength(pwd):
	#print 'Password:', pwd, 
	strength, _ = passwordmeter.test(pwd)
	return str(strength)
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
#Docker launch
def PickUnusedPort():
	s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	s.bind(('localhost', 0))
	addr, port = s.getsockname()
	s.close()
	#print 'Found free:', port
	return str(port)

#Bring up a unique docker instance for a reconstruction
@app.route('/docker', methods=['POST'])
@auth.login_required
def launch_neo4j_docker():
	try:
		json_data = request.get_json(force=True)
		username = json_data['username']
		recon_name = json_data['recon_name']
		if len(recon_name) == 0:
			recon_name = 'nemo'
		fsid = "EMPATHY_" + str(uuid.uuid4())
		##print 'DOCKER: Create Neo4j-Docker instance:', fsid
		freePort = PickUnusedPort()
		##print 'DOCKER: found port', freePort
		port_bind = "--publish=" + freePort + ":7474"
		docker = ["docker", "run", port_bind, "--volume=/Users/paul/Fire/FS/" + fsid + ":/data", "--name", fsid, "neo4j:3.0.0"]
		try:
			#print 'Start subprocess'
			subprocess.Popen(docker)
			print 'DOCKER: Created on port:', freePort, fsid, 'for', username
			user = user_credentials[username]
			reconstructions = user['reconstructions']
			#print 'Reconstructions (before)', reconstructions
			recon = {"name":recon_name,"port":freePort,"fsid":fsid}
			#print 'Recon', recon
			reconstructions.append(recon)
			#print 'Reconstructions (after)', reconstructions
			user["reconstructions"] = reconstructions
			#print 'User', user
			user_credentials[username] = user
			#print 'Credentials', user_credentials
			#Update the pickle
			pickle.dump( user_credentials, open( "users.pickle", "wb" ) )
			return json.dumps([freePort, fsid])
		except:
			return json.dumps(['error', False])
	except:
		return json.dumps(['unspecified', False])

#Check the Docker instance on <port> is alive and kicking
@app.route('/checkLive', methods=['POST'])
@auth.login_required
def checkLive():
	print 'Check live'
	json_data = request.get_json(force=True)
	username = json_data['username']
	recon_name = json_data['recon_name']
	port = json_data['port']
	if json_data['notes']:
		notes = json_data['notes']
	else:
		notes = 'None provided'

	url = 'http://localhost:' + port + '/db/data/transaction/commit'
	print 'Hitting', url
	i = 1
	live = [False, i]
	while not live[0] and i <= 60:
		try:
			response = requests.post(url)
			#print i, 'Live on', port
			live = [True,i]
		except:
			#print i, 'Not live on', port
			live = [False,i]
		time.sleep(1)
		i = i + 1
	if live[0]:
		print 'Alive'
		#Change default password to EMPATHY password
		#print 'Changing default password'
		headers = {'content-type': 'application/json'}
		pwd_response = requests.post('http://localhost:' + port + '/user/neo4j/password', auth=('neo4j','neo4j'), data=json.dumps({"password":NEO_PASSWORD}), headers=headers)
		
		#Create recon node
		time.sleep(3) #Just wait a minute
		#print 'Create reconstruction node'
		node_id = str(uuid.uuid4())
		properties = {"id": node_id, "founder":username,"name":recon_name,"notes":notes}
		#print properties
		cypher = "CREATE (n:recon {props}) RETURN n"
		parameters = {"props": properties}
		response = send_cypher(cypher,parameters,port)
		#print response.json()
		
		return json.dumps([True,'Created ' + recon_name + ' on port ' + port])
	else:
		print 'Given up'
		return json.dumps([False,'Could not create ' + recon_name + ' on port ' + port])

#Fetch the user's list of recons
@app.route('/userRecons', methods=['POST'])
@auth.login_required
def userRecons():
	print 'Fetch user reconstructions'
	json_data = request.get_json(force=True)
	username = json_data['username']
	#print 'User recons username:', username
	if username in user_credentials:
		user = user_credentials[username]
		reconstructions = user["reconstructions"]
		return json.dumps([True,reconstructions])
	else:
		return json.dumps([False,'Could not find user'])

#Add user to recon (part 1)
@app.route('/msgNewUser', methods=['POST'])
@auth.login_required
def msgNewUser():
	print 'Add user to reconstruction 1'
	json_data = request.get_json(force=True)
	username = json_data["username"]
	recon = json_data["recon"]
	port = json_data["port"]
	fsid = json_data['fsid']
	user_email = json_data['user_email']
	new_token = str(uuid.uuid4()) + '_' + str(uuid.uuid4())
	tokens[new_token] = {"recon":recon,"port":port,"fsid":fsid}
	try:
		msg = Message("You have been invited to collaborate on a metabolic reconstruction for " + recon, sender=("EMPATHY simple metabolic networks","empathy.manchester@gmail.com"),recipients=[user_email])
		msg.html = 	'<b>EMPATHY Metabolic Network Reconstruction</b><br><br>You have been invited by ' + username + ' to collaborate on a reconstruction of ' + recon + ' metabolism.<br><br>Follow the link below to start collaborating.<br><br>https://www.metabolicjamboree.co.uk/index.html?token=' + new_token
		mail.send(msg)
		time.sleep(2)
		return json.dumps([True,new_token])
	except:
		time.sleep(2)
		return json.dumps([False,"Mail not sent"])

#Add user to recon (part 2)
@app.route('/addNewUser', methods=['POST'])
@auth.login_required
def addNewUser():
	print 'Add user to reconstruction 2'
	json_data = request.get_json(force=True)
	username = json_data["username"]
	token = json_data["token"]
	#print username, token
	user = user_credentials[username]
	recons = user["reconstructions"]
	try:
		#Look up token	
		print 'Resolve token'
		recon_details = tokens[token]
		print 'Resolve recon details', recon_details
		recons.append(recon_details)
		print 'Reconstructions', recons
		user["reconstructions"] = recons
		print 'Updated user', user
		user_credentials[username] = user
		#Write to pickle
		#print 'Write update credentials'
		pickle.dump( user_credentials, open( "users.pickle", "wb" ) )
		#Remove token
		#print 'Remove token', token
		del tokens[token]
		return json.dumps([True,recon_details])
	except:
		return json.dumps([False,'Could not add collaborator'])
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
#Node and relationship blanks
nodeBlank = {"id":"", "type":"", "source":"", "sourceId":"", "name":"",	"synonyms":[], "inCompartment":"", "is":[], "is_a":[], "isDescribedBy":[], "isVersionOf":[], "property":[], "tags":[], "sandbox":[], "references":[], "notifications":[]	}
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
def fetchList(label,port):
	#Don't look for compartment (no such thing for reactions)	
	print 'Fetch list:', label
	new_list = []
	if 'molecule' in label:
		print 'Fetch ', label
		cypher = "MATCH (n:" + label + ") RETURN n.id AS id, n.inCompartment AS compartment, n.name AS name, n.tags AS tags ORDER BY name"
		response = send_cypher(cypher,{},port)
		cypher_response = response.json()
		for row in cypher_response["results"][0]["data"]:
			id = row["row"][0]
			compartment = str(row["row"][1])
			name = row["row"][2]
			tags = row["row"][3]
			name = name + '_[' + compartment + ']'
			new_list.append({"id":id,"name":name,"tags":tags})
	else:
		print 'Fetch', label
		cypher = "MATCH (n:" + label + ") RETURN n.id AS id, n.name AS name, n.tags AS tags ORDER BY name"
		response = send_cypher(cypher,{},port)
		cypher_response = response.json()
		for row in cypher_response["results"][0]["data"]:
			id = row["row"][0]
			name = row["row"][1]
			tags = row["row"][2]
			new_list.append({"id":id,"name":name,"tags":tags})
	print 'New_list:', len(new_list)
	return new_list

@app.route('/makeConnection', methods=['POST'])
@auth.login_required
def makeConnection():
	json_data = request.get_json(force=True)
	port = json_data['port']
	record_handle = json_data['record_handle'] #<port>_<record id>
	reaction = json_data['reaction']
	molecule = json_data['molecule']
	role = json_data['role']

	#Create connection
	cypher = "MATCH (r:reaction {id:'" + reaction + "'}), (e:molecule {id:'" + molecule + "'}) CREATE (r)-[:" + role + "]->(e)"
	parameters = {}
	response = send_cypher(cypher,parameters,port)
	#print response.json()
		
	#Emit updated record
	#print 'BROADCAST:', record_handle
	send_message(record_handle, {"key":"","value":"","id":reaction}, port)
	return json.dumps(True)

'''
@app.route('/fetchSYNBIOCHEM', methods=['POST'])
@auth.login_required
def fetchFromSYNBIOCHEM():
	json_data = request.get_json(force=True)
	port = json_data['port']
	ncbi = json_data['ncbi']
	print 'Fetch from SYNBIOCHEM-DB:', ncbi

	#--------------------------------------------
	#ENZYMES
	#Get enzymes
	#print 'Get enzymes'
	cypher = "MATCH (n:Organism)<-[*]-(o:Organism)-[]-(e:Enzyme) WHERE (n.taxonomy='" + ncbi + "') RETURN DISTINCT(e) AS enzyme"
	headers = {'content-type': 'application/json'}
	query = { "statements" : [ { "statement" : cypher } ] }
	synbiochem = requests.post('http://db.synbiochem.co.uk/db/data/transaction/commit', data=json.dumps(query), headers=headers)

	#print 'Put enzymes'
	is_keys = ["SYNBIOCHEM-DB"]
	is_a_keys = []
	property_keys = []
	isDescribedBy_keys = ["uniprot"]
	isVersionOf_keys = []
	enzymes = synbiochem.json()
	for enz in enzymes["results"][0]["data"]:
		enzyme = enz["row"][0]
		
		#Create new enzyme
		properties = copy.deepcopy(nodeBlank)
		properties["id"] = str(uuid.uuid4())
		properties["type"] = "macromolecule"
		properties["source"] = "SYNBIOCHEM-DB"
		properties["sourceId"] = enzyme["entry"]
	 	properties["tags"] = ["macromolecule"]

		#Source into 'is'
		enzyme["SYNBIOCHEM-DB"] = enzyme["entry"]
		
		#Name
		if "name" in enzyme:
			properties["name"] = enzyme["name"]
		else:
			properties["name"] = enzyme["entry"]

		#Synonyms
		if "names" in enzyme:
			synonyms = enzyme["names"]
			synonyms = list( set( synonyms.split(";") ) )
			new_synonyms = []
			for syn in synonyms:
				new_synonyms.append('?'+syn)
			properties["synonyms"] = new_synonyms

		#Keys
		for key in enzyme.keys():
			if key in is_keys:
				properties["is"].append( json.dumps([key,'?'+enzyme[key]]) )
			elif key in isDescribedBy_keys:
				properties["isDescribedBy"].append( json.dumps([key,'?'+enzyme[key]]) )
			elif key in isVersionOf_keys:
				properties["isVersionOf"].append( json.dumps([key,'?'+enzyme[key]]) )
			elif key in property_keys:
				properties["property"].append( json.dumps([key,enzyme[key]]) )
			elif key in is_a_keys:
				properties["is_a"].append( json.dumps([key,'?'+enzyme[key]]) )
			else:
				pass
		cypher = "CREATE (n:molecule {props}) RETURN n"
		parameters = {"props": properties}
		response = send_cypher(cypher,parameters,port)
		
	#Broadcast list of molecules
	#print 'Broadcast enzymes (1)'
	cypher2 = "MATCH (n:molecule) RETURN n.id AS id, n.name AS name, n.inCompartment AS compartment, n.tags AS tags"
	response2 = send_cypher(cypher2,{},port)
	cypher_response = response2.json()
	new_list = []
	for row in cypher_response["results"][0]["data"]:
		id = row["row"][0]
		name = row["row"][1]
		compartment = str( row["row"][2] )
		tags = row["row"][3]
		name = name + '_[' + compartment[0:2] + ']'
		new_list.append({"id":id,"name":name,"tags":tags})

	#print 'Broadcast enzymes (2)', port+"_molecule"
	send_message(port+"_molecule", new_list,  port)
	#--------------------------------------------

	#--------------------------------------------
	#CHEMICALS
	#Get chemicals
	#print 'Get chemicals'
	cypher = "MATCH (n:Organism)<-[*]-(o:Organism)-[]-(e:Enzyme)-[]-(r:Reaction)-[t]-(c:Chemical) WHERE (n.taxonomy='" + ncbi + "') RETURN DISTINCT(c) AS chemical"
	headers = {'content-type': 'application/json'}
	query = { "statements" : [ { "statement" : cypher } ] }
	synbiochem = requests.post('http://db.synbiochem.co.uk/db/data/transaction/commit', data=json.dumps(query), headers=headers)
	chemicals = synbiochem.json()

	#print 'Put chemicals'
	is_keys = []
	isDescribedBy_keys = ["mnx","metacyc","knapsack","chebi","inchi","smiles","kegg.compound","lipidmaps","hmdb","cas","wikipedia.en"]
	isVersionOf_keys = []
	property_keys = ["monoisotopic_mass","formula","charge"]
	for chem in chemicals["results"][0]["data"]:
		chemical = chem["row"][0]

		#Create new chemical
		properties = copy.deepcopy(nodeBlank)
		properties["id"] = str(uuid.uuid4())
		properties["type"] = "simple chemical"
		properties["source"] = "SYNBIOCHEM-DB"
		properties["sourceId"] = chemical["id"]
	 	properties["tags"] = ["simple chemical"]
		
		#Source into 'is'
		chemical["SYNBIOCHEM-DB"] = chemical["id"]

		#Name		
		if "name" in chemical:
			properties["name"] = chemical["name"]
		else:
			properties["name"] = chemical["id"]
		
		#Synonyms
		if "names" in chemical:
			synonyms = chemical["names"]
			synonyms = list( set( synonyms ) )
			new_synonyms = []
			for syn in synonyms:
				new_synonyms.append('?'+syn)
			properties["synonyms"] = new_synonyms
	
		for key in chemical.keys():
			if key in is_keys:
				properties["is"].append( json.dumps([key,'?'+chemical[key]]) )
			elif key in isDescribedBy_keys:
				properties["isDescribedBy"].append( json.dumps([key,'?'+chemical[key]]) )
			elif key in isVersionOf_keys:
				properties["isVersionOf"].append( json.dumps([key,'?'+chemical[key]]) )
			elif key in property_keys:
				properties["property"].append( json.dumps([key,chemical[key]]) )
			elif key in is_a_keys:
				properties["is_a"].append( json.dumps([key,'?'+chemical[key]]) )
			else:
				pass
		cypher = "CREATE (n:molecule {props}) RETURN n"
		parameters = {"props": properties}
		response = send_cypher(cypher,parameters,port)
	
	#Broadcast list of molecules
	#print 'Broadcast chemicals (1)'
	cypher2 = "MATCH (n:molecule) RETURN n.id AS id, n.name AS name, n.inCompartment AS compartment, n.tags AS tags"
	response2 = send_cypher(cypher2,{},port)
	cypher_response = response2.json()
	new_list = []
	for row in cypher_response["results"][0]["data"]:
		id = row["row"][0]
		name = row["row"][1]
		compartment = str(row["row"][2])
		tags = row["row"][3]
		name = name + '_[' + compartment[0:2] + ']'
		new_list.append({"id":id,"name":name,"tags":tags})
	##print new_list

	#print 'Broadcast chemicals (2)', port+"_molecule"
	send_message(port+"_molecule", new_list,  port)
	#--------------------------------------------

	#--------------------------------------------
	#REACTIONS
	#Get reactions
	#print 'Get reactions'
	cypher = "MATCH (n:Organism)<-[*]-(o:Organism)-[]-(e:Enzyme)-[]-(r:Reaction) WHERE (n.taxonomy='" + ncbi + "') RETURN DISTINCT(r) AS reaction"
	headers = {'content-type': 'application/json'}
	query = { "statements" : [ { "statement" : cypher } ] }
	synbiochem = requests.post('http://db.synbiochem.co.uk/db/data/transaction/commit', data=json.dumps(query), headers=headers)
	reactions = synbiochem.json()

	#print 'Put reactions'
	is_keys = []
	isDescribedBy_keys = ["rhea","kegg.reaction","ec","seed","metacyc","reactome","bigg.reaction"]
	isVersionOf_keys = []
	property_keys = ["balance"]
	for rxn in reactions["results"][0]["data"]:
		reaction = rxn["row"][0]

		#Create new reaction
		properties = copy.deepcopy(nodeBlank)
		properties["id"] = str(uuid.uuid4())
		properties["type"] = "reaction"
		properties["source"] = "SYNBIOCHEM-DB"
		properties["sourceId"] = reaction["id"]
	 	properties["tags"] = ["reaction"]
	 	properties["reversible"] = "true";
		
		#Source into 'is'
		reaction["SYNBIOCHEM-DB"] = reaction["id"]

		if "name" in reaction:
			properties["name"] = reaction["name"]
		else:
			properties["name"] = reaction["id"]

		for key in reaction.keys():
			if key in is_keys:
				properties["is"].append( json.dumps([key,'?'+reaction[key]]) )
			elif key in isDescribedBy_keys:
				properties["isDescribedBy"].append( json.dumps([key,'?'+reaction[key]]) )
			elif key in isVersionOf_keys:
				properties["isVersionOf"].append( json.dumps([key,'?'+reaction[key]]) )
			elif key in property_keys:
				properties["property"].append( json.dumps([key,reaction[key]]) )
			elif key in is_a_keys:
				properties["is_a"].append( json.dumps([key,'?'+reaction[key]]) )
			else:
				pass
		cypher = "CREATE (n:reaction {props}) RETURN n"
		parameters = {"props": properties}
		response = send_cypher(cypher,parameters,port)
	#--------------------------------------------

	#--------------------------------------------
	#MODIFIERS
	#Build enzyme-reaction links
	#print 'Get modifiers'
	cypher = "MATCH (n:Organism)<-[*]-(o:Organism)-[]-(e:Enzyme)-[l]-(r:Reaction) WHERE (n.taxonomy='" + ncbi + "') RETURN COLLECT(l.source) AS link, r.id AS reactionId, e.entry AS enzymeId"
	headers = {'content-type': 'application/json'}
	query = { "statements" : [ { "statement" : cypher } ] }
	synbiochem = requests.post('http://db.synbiochem.co.uk/db/data/transaction/commit', data=json.dumps(query), headers=headers)
	links = synbiochem.json()

	#print 'Put modifiers'
	for l in links["results"][0]["data"]:
		properties = l["row"][0]
		reactionId = l["row"][1]
		enzymeId = l["row"][2]
		##print 'LINKS', properties, reactionId, enzymeId
		cypher = "MATCH (r:reaction {sourceId:'" + reactionId + "'}), (e:molecule {sourceId:'" + enzymeId + "'}) CREATE (r)-[:hasModifier {props}]->(e)"
		parameters = {"props": {"source":properties} }
		response = send_cypher(cypher,parameters,port)
	#--------------------------------------------
	
	#--------------------------------------------
	#REACTANTS/PRODUCTS
	#Build reaction-metabolite links
	#print 'Get reactants and products'
	cypher = "MATCH (n:Organism)<-[*]-(o:Organism)-[]-(e:Enzyme)-[]-(r:Reaction)-[l]-(c:Chemical) WHERE (n.taxonomy='" + ncbi + "') RETURN l AS link, TYPE(l) AS label, r.id AS reactionId, c.id AS chemicalId"
	headers = {'content-type': 'application/json'}
	query = { "statements" : [ { "statement" : cypher } ] }
	synbiochem = requests.post('http://db.synbiochem.co.uk/db/data/transaction/commit', data=json.dumps(query), headers=headers)
	links = synbiochem.json()
	##print links

	#print 'Put reactants and products'
	#Some weird flattening required
	uniq = []
	for l in links["results"][0]["data"]:
		properties = l["row"][0]
		stoichiometry = l["row"][0]["stoichiometry"]
		if stoichiometry > 0:
			label = 'hasProduct'
		else:
			label = 'hasReactant'
		reactionId = l["row"][2]
		chemicalId = l["row"][3]
		uniq.append( json.dumps([properties, label, reactionId, chemicalId]) )
	uniq = list(set(uniq))
	
	for rel in uniq:
		r = json.loads(rel)
		properties = r[0]
		label = r[1]		
		reactionId = r[2]
		chemicalId = r[3]				
		cypher = "MATCH (r:reaction {sourceId:'" + reactionId + "'}), (c:molecule {sourceId:'" + chemicalId + "'}) CREATE (r)-[:" + label + " {props}]->(c)"
		parameters = {"props": properties}
		response = send_cypher(cypher,parameters,port)
	#--------------------------------------------

	#--------------------------------------------
	#CAUTION! ROGUE REACTIONS! (Reactions without participants)
	#print 'Caution! Rogue robots! Removing reactions without reactants or products'

	#Remove links on dangling reactions
	cypher = "MATCH (r:reaction) WHERE NOT (r)-[:hasReactant|hasProduct]-() MATCH (r)-[l]-() DELETE(l)"
	parameters = {}
	response = send_cypher(cypher,parameters,port)

	#Remove dangling reactions	
	cypher = "MATCH (r:reaction) WHERE NOT (r)-[:hasReactant|hasProduct]-() DELETE(r)"
	parameters = {}
	response = send_cypher(cypher,parameters,port)
	#--------------------------------------------

	#--------------------------------------------
	#Broadcast list of reactions
	#print 'Broadcast reactions (1)'
	cypher2 = "MATCH (n:reaction) RETURN n.id AS id, n.name AS name, n.name_level AS name_level, n.tags AS tags, n.tags_level AS tags_level"
	response2 = send_cypher(cypher2,{},port)
	cypher_response = response2.json()
	new_list = []
	for row in cypher_response["results"][0]["data"]:
		id = row["row"][0]
		name = row["row"][1]
		name_level = row["row"][2]
		tags = row["row"][3]
		tags_level = row["row"][4]
		new_list.append({"id":id,"name":name,"tags":tags})
	##print new_list
	#print 'Broadcast reactions (2)', port+"_reaction"	
	send_message(port+"_reaction", new_list,  port)
	#--------------------------------------------

	print 'SYNBIOCHEM-DB import complete'
	return json.dumps(True)
'''

@app.route('/importSBML', methods=['POST'])
@auth.login_required
def importSBML():
	json_data = request.get_json(force=True)
	port = json_data['port']
	sbml = json_data['sbml']

	#Build handle for passing progress messages
	sbml_handle = port + '_sbml'

	print 'Import SBML...'
	try:
		send_message(sbml_handle, 'Parsing SBML',  port)
		cyphers = libsbml.sbml2cyphers(sbml)
		send_message(sbml_handle, 'SBML parsed',  port)
		print 'libSBML parsing successful'
		#Push cyphers
		print 'Pushing cyphers to', port, '...'
		try:
			i = 1
			j = 1
			total = len(cyphers)
			for cypher in cyphers:
				parameters = {"props": cypher[1]}
				response = send_cypher(cypher[0],parameters,port)
				if j == 100:
					what_are_we_adding = cypher[2]
					sbml_message = 'Adding to database: ' + str(i) + '/' + str(total) + ' (' + what_are_we_adding + 's)'
					send_message(sbml_handle, sbml_message,  port)
					j = 1
				else:
					j = j + 1
				i = i + 1
		except Exception, e:
			send_message(sbml_handle,'Upload error',  port)
			print 'Could not push cyphers'		

		##Push notifications to socket
		print 'Disseminating lists via web socket...'
		#Compartments
		try:
			cypher2 = "MATCH (n:compartment) RETURN n.id AS id, n.name AS name, n.tags AS tags"
			response2 = send_cypher(cypher2,{},port)
			cypher_response = response2.json()
			new_list = []
			for row in cypher_response["results"][0]["data"]:
				id = row["row"][0]
				name = row["row"][1]
				tags = row["row"][2]
				new_list.append({"id":id,"name":name,"tags":tags})
			send_message(port+"_compartment", new_list,  port)				
		except Exception, e:
			print 'Could not push compartments', str(e)

		#Molecules
		try:
			cypher2 = "MATCH (n:molecule) RETURN n.id AS id, n.name AS name, n.inCompartment AS compartment, n.tags AS tags"
			response2 = send_cypher(cypher2,{},port)
			cypher_response = response2.json()
			new_list = []
			for row in cypher_response["results"][0]["data"]:
				id = row["row"][0]
				name = row["row"][1]
				compartment = str( row["row"][2] )
				tags = row["row"][3]
				name = name + '_[' + compartment + ']'
				new_list.append({"id":id,"name":name,"tags":tags})
			send_message(port+"_molecule", new_list,  port)							
		except Exception, e:
			print 'Could not push molecules', str(e)

		#Reactions
		try:
			cypher2 = "MATCH (n:reaction) RETURN n.id AS id, n.name AS name, n.tags AS tags"
			response2 = send_cypher(cypher2,{},port)
			cypher_response = response2.json()
			new_list = []
			for row in cypher_response["results"][0]["data"]:
				id = row["row"][0]
				name = row["row"][1]
				tags = row["row"][2]
				new_list.append({"id":id,"name":name,"tags":tags})							
			send_message(port+"_reaction", new_list,  port)
		except Exception, e:
			print 'Could not push compartments', str(e)

		send_message(sbml_handle, 'Complete',  port)
		return json.dumps(True)
	except:
		send_message(sbml_handle, 'Error in SBML upload',  port)
		print 'Could not convert to libSBML object'
		return json.dumps(False)

	print 'SBML import complete'
	return json.dumps(True)

@app.route('/createNode', methods=['POST'])
@auth.login_required
def createNode():
	json_data = request.get_json(force=True)
	port = json_data['port']
	message_handle = json_data['message_handle'] #<port>_<type>
	label = json_data['label']
	print 'Create node', message_handle

	#Create new node
	properties = copy.deepcopy(nodeBlank)
	properties["id"] = str(uuid.uuid4())
	properties["type"] = label
	properties["source"] = "curator"
	properties["sourceId"] = "curator"
	properties["tags"] = [label]
	properties["name"] = "new "+label
	
	cypher = "CREATE (n:" + label + " {props}) RETURN n"
	parameters = {"props": properties}
	response = send_cypher(cypher,parameters,port)
	#print response.json()
	
	#Disseminate new list of nodes
	new_list = fetchList(label,port)
	send_message(message_handle, new_list,  port)
	return properties["id"]

@app.route('/subcell', methods=['POST'])
@auth.login_required
def subCell():
	json_data = request.get_json(force=True)
	port = json_data['port']
	cell = json_data['cell'] #Use this to build in logic for different cell types
	message_handle = json_data['message_handle'] #<port>_<type>
	#print 'Create cellular organisation', message_handle

	cells = {}
	cells["bacterium"] = {	"extracellular":{"name":"extracellular space","inCompartment":"null"}, "membrane_cell":{"name":"cell membrane","inCompartment":"extracellular"}, "cytosol":{"name":"cytosol","inCompartment":"membrane_cell"}	}
	cells["eukaryote"] = {	"extracellular":{"name":"extracellular","inCompartment":"null"}, "membrane_cell":{"name":"cell membrane","inCompartment":"extracellular"}, "cytosol":{"name":"cytosol","inCompartment":"membrane_cell"}	}
	cells["system"] = {	"extracellular":{"name":"extracellular","inCompartment":"null"}, "membrane_cell":{"name":"cell membrane","inCompartment":"extracellular"}, "cytosol":{"name":"cytosol","inCompartment":"membrane_cell"}	}

	try:
		organelles = cells[cell]
	except:
		return "No such cell type within generic patterns"

	for org in organelles.keys():
		props = organelles[org]
		#Create new node
		properties = copy.deepcopy(nodeBlank)
		properties["id"] = org
		properties["type"] = "organelle"
		properties["source"] = "autobot"
		properties["sourceId"] = "autobot_"+org
		properties["tags"] = []
		properties["name"] = props["name"]
		properties["inCompartment"] = props["inCompartment"]
	
		cypher = "CREATE (n:compartment {props}) RETURN n"
		parameters = {"props": properties}
		response = send_cypher(cypher,parameters,port)
		#print response.json()
		
	#Disseminate new list of nodes
	new_list = fetchList('compartment',port)
	send_message(message_handle, new_list,  port)
	return json.dumps(True)

#Fetch pane 1 list of type <label>
@app.route('/listNode', methods=['POST'])
@auth.login_required
def listNode():
	json_data = request.get_json(force=True)
	label = json_data['label']
	message_handle = json_data['message_handle'] #<port>_<type>
	port = json_data['port']
	print 'listNode', message_handle, port

	#Get new list of nodes
	new_list = fetchList(label,port)

	if message_handle != '':
		send_message(message_handle, new_list, port)
	return json.dumps(new_list)

@app.route('/destroyNode', methods=['POST'])
@auth.login_required
def destroyNode():
	json_data = request.get_json(force=True)
	label = json_data['label']
	target = json_data['target']
	message_handle = json_data['message_handle'] #<port>_<type>
	record_handle = json_data['record_handle'] #<port>_<record_id>
	port = json_data['port']
	
	try:
		#Destroy node's relations
		rel_cypher = "MATCH (n {id:'" + target + "'})-[r]-() DELETE r"
		print 'DESTROY:', rel_cypher
		parameters = {}
		response = send_cypher(rel_cypher,parameters,port)
		#Destroy node
		node_cypher = "MATCH (n {id:'" + target + "'}) DELETE n"
		print 'DESTROY:', node_cypher
		parameters = {}
		response = send_cypher(node_cypher,parameters,port)

		#Push updated list to subscribers
		new_list = fetchList(label,port)
		send_message(message_handle, new_list,  port)		
		return json.dumps(True)
	except:
		return json.dumps(False)

#Use this to take participants out of reactions
@app.route('/destroyEdge', methods=['POST'])
@auth.login_required
def destroyEdge():
	json_data = request.get_json(force=True)
	targetA = json_data['targetA']
	targetB = json_data['targetB']
	record_handle = json_data['record_handle'] #<port>_<record_id>
	port = json_data['port']

	try:
		#Destroy relations between A (reaction) and B
		rel_cypher = "MATCH (n {id:'" + targetA + "'})-[r]-(m {id:'" + targetB + "'}) DELETE r"
		print 'DESTROY:', rel_cypher
		parameters = {}
		response = send_cypher(rel_cypher,parameters,port)
		
		##Get updated version of node A (reaction)		
		#Fetch core of node
		cypher = "MATCH (n) WHERE n.id='" + targetA + "' RETURN n"
		response = send_cypher(cypher,{},port)
		data = response.json()
		record = data["results"][0]["data"][0]["row"][0]

		#Signify change to subscribers so they can fetch whole new record
		print 'Broadcast change to', record_handle
		send_message(record_handle, {'key':'','id': targetA},  port)
		
		return json.dumps(True)
	except:
		return json.dumps(False)

#Retrieve a node
@app.route('/fetchSelection', methods=['POST'])
@auth.login_required
def fetchSelection():
	json_data = request.get_json(force=True)
	port = json_data['port']
	selection = json_data['selection']
	label = json_data['label']
	#print 'Fetch selection:', port, selection, label

	#Fetch core of node
	cypher = "MATCH (n) WHERE n.id='" + selection + "' RETURN n"
	response = send_cypher(cypher,{},port)
	data = response.json()
	record = data["results"][0]["data"][0]["row"][0]

	if label == 'reaction':
		#Fetch edges
		cypher = "MATCH (r:reaction)-[l]->(m:molecule) WHERE r.id='" + selection + "' RETURN l,TYPE(l) AS type,m.id AS moleculeId, m.name AS moleculeName, m.inCompartment AS compartment"
		edgeResponse = send_cypher(cypher,{},port)
		edgeData = edgeResponse.json()

		reactants = []
		modifiers = []
		products = []
		for row in edgeData['results'][0]['data']:
			edge = row['row']
			print 'Reaction:', edge
			edgeProperties = edge[0]
			type = edge[1]
			moleculeId = edge[2]
			moleculeName = edge[3]
			moleculeCompartment = str( edge[4] )
			if type == 'hasReactant':
				reactants.append({"id":moleculeId,"name":moleculeName,"properties":edgeProperties,"compartment":moleculeCompartment})
			elif type == 'hasModifier':
				modifiers.append({"id":moleculeId,"name":moleculeName,"properties":edgeProperties,"compartment":moleculeCompartment})
			elif type == 'hasProduct':
				products.append({"id":moleculeId,"name":moleculeName,"properties":edgeProperties,"compartment":moleculeCompartment})
			else:
				pass
		record["listOfReactants"] = reactants
		record["listOfModifiers"] = modifiers
		record["listOfProducts"] = products
	elif label == 'molecule':
		#Fetch edges
		cypher = "MATCH (m:molecule)<-[l:hasReactant|hasProduct|hasModifier]-(r:reaction) WHERE m.id='" + selection + "' RETURN TYPE(l), r.id, r.name"
		edgeResponse = send_cypher(cypher,{},port)
		edgeData = edgeResponse.json()

		reactants = []
		modifiers = []
		products = []
		for row in edgeData['results'][0]['data']:
			edge = row['row']
			type = edge[0]
			rxnId = edge[1]
			rxnName = edge[2]
			if type == 'hasReactant':
				reactants.append({"id":rxnId,"name":rxnName})
			elif type == 'hasModifier':
				modifiers.append({"id":rxnId,"name":rxnName})
			elif type == 'hasProduct':
				products.append({"id":rxnId,"name":rxnName})
			else:
				pass
		reactions = {}
		reactions["asReactant"] = reactants
		reactions["asModifier"] = modifiers
		reactions["asProduct"] = products
		record["reactions"] = reactions

	#print record
	return json.dumps(record)

@app.route('/updateText', methods=['POST'])
@auth.login_required
def updateText():
	json_data = request.get_json(force=True)
	port = json_data['port']
	record_handle = json_data['record_handle'] #<port>_<id>
	id = json_data['id']
	key = json_data['key']
	value = json_data['value']
	#print 'UPDATE\tid:', id, 'key:', key, 'value:', value
	
	#Update node
	cypher = "MATCH (n) WHERE n.id='" + id + "' SET n." + key + "={value} RETURN n"
	parameters = {"value":value}
	response = send_cypher(cypher,parameters,port)
	#print 'BROADCAST:', record_handle, key, value
	send_message(record_handle, {'key':key,'value':value, 'id': id},  port)
	return json.dumps(True)

@app.route('/listPop', methods=['POST'])
@auth.login_required
def listPop():
	json_data = request.get_json(force=True)
	port = json_data['port']
	record_handle = json_data['record_handle'] #<port>_<id>
	id = json_data['id']
	key = json_data['key']
	value = json_data['value']
	#value = list(set(value))
	
	#Update node
	cypher = 'MATCH (n) WHERE n.id="' + id + '" SET n.' + key + '={value} RETURN n'
	parameters = {"value":value}
	response = send_cypher(cypher,parameters,port)
	#print 'BROADCAST:', record_handle, key, value
	send_message(record_handle, {'key':key,'value':value, 'id': id},  port)
	return json.dumps(True)

@app.route('/listPush', methods=['POST'])
@auth.login_required
def listPush():
	json_data = request.get_json(force=True)
	port = json_data['port']
	record_handle = json_data['record_handle'] #<port>_<id>
	id = json_data['id']
	key = json_data['key']
	value = json_data['value']
	#value = list(set(value))
	#print 'listPush:', key, value
	
	#Update node
	cypher = 'MATCH (n) WHERE n.id="' + id + '" SET n.' + key + '={value} RETURN n'
	parameters = {"value":value}
	response = send_cypher(cypher,parameters,port)
	#print 'BROADCAST:', record_handle, key, value
	send_message(record_handle, {'key':key,'value':value, 'id': id},  port)
	return json.dumps(True)

#For finding reactions in the reaction builder
@app.route('/queryMolecules', methods=['POST'])
@auth.login_required
def queryMolecules():
	json_data = request.get_json(force=True)
	port = json_data['port']
	query = json_data['query']
	
	#Update node
	cypher = "MATCH (n:molecule) WHERE n.name =~ '.*" + query + ".*' RETURN n.id AS id, n.name AS name, n.inCompartment AS  compartment"
	parameters = {}
	response = send_cypher(cypher,parameters,port)
	results = response.json()
	##print results
	matches = [];
	for row in results["results"][0]["data"]:
		match = row["row"]
		matches.append(match)
	return json.dumps(matches)
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
#Action stations

def actionPush(port,record_handle,record,key,value,username,password,message):
	#Push out new list to record
	#print 'Update record'
	headers = {'content-type': 'application/json'}
	url = 'https://' + username + ':' + password + '@' + HOST + ':' + str(PORT) + '/listPush'
	bundle = {"port":port,"record_handle":record_handle,"id":record['id'],"key":key,"value":value}
	requests.post(url, data=json.dumps(bundle), headers=headers)

	#Push out message notification
	#print 'Send notification message'
	url = 'https://' + username + ':' + password + '@' + HOST + ':' + str(PORT) + '/listPush'
	notifications = record['notifications']
	notifications.append(message)
	bundle = {"port":port,"record_handle":record_handle,"id":record['id'],"key":"notifications","value":notifications}
	requests.post(url, data=json.dumps(bundle), headers=headers)	

@app.route('/actionMolecule', methods=['POST'])
@auth.login_required
def actionMolecule():
	json_data = request.get_json(force=True)
	username = json_data['username']
	password = json_data['password']
	record_handle = json_data['record_handle'] #<port>_<record_id>
	record = json_data['record']
	action = json_data['action']
	port = json_data['port']

	if action == "smallMoleculeIdentifier":
		#First run action
		try:
			name = record['name']
			result = actions.smallMoleculeIdentifier(name)
			#print 'ACTION:', result
		except:
			return json.dumps(False)
		
		#Then post result (fetch current then add updates)
		oldList = []
		for row in record['isDescribedBy']:
			oldList.append(json.dumps(row))
		for row in result:
			oldList.append(json.dumps(row))
		newList = list(set(oldList))
		
		#Push out new list and push notification
		actionPush(port,record_handle,record,"isDescribedBy",newList,username,password,"New molecule identifiers added (isDescribedBy)")			
		return json.dumps(True)

	elif action == "smallMoleculeSynonyms":
		#First run action
		try:
			name = record['name']
			result = actions.smallMoleculeSynonyms(name)
			#print 'ACTION:', result
		except:
			return json.dumps(False)
		
		#Then post result (fetch current then add updates)
		oldList = record['synonyms']
		#print "oldList:", oldList
		for row in result:
			oldList.append(row)
		#print "newList 1:", oldList
		newList = list(set(oldList))
		#print "newList 2:", newList
		#Push out new list and push notification
		actionPush(port,record_handle,record,"synonyms",newList,username,password,"New synonyms added")			
		return json.dumps(True)

	return json.dumps(True)
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
##STATIC PAGES
#Landing page
@app.route('/')
def index():
	return app.send_static_file('index.html')
#-----------------------------------------------------------------------------------------

# Run the app.
if __name__ == '__main__':
	#app.run(debug=True, port=PORT, threaded=True)
	#https://github.com/enginebai/PyMessager
	context = ('certs/fullchain.pem', 'certs/privkey.pem')
	app.run(host='0.0.0.0',debug=False, ssl_context=context, port=PORT, threaded=True)
	#app.run(host='0.0.0.0',debug=False, port=PORT, threaded=True)
