from flask import abort, Flask, url_for, request, make_response, send_file, Markup
from flask_restful import Resource, Api, reqparse
import json, logging, requests, xmltodict, uuid, time
from requests.auth import HTTPBasicAuth
from flask_cors import CORS
from rdkit import Chem
from rdkit.Chem import AllChem


'''
NB: Neo4j authorization disabled (authorization implemented since 2.1.7)

See - http://stackoverflow.com/questions/29096616/how-to-disable-basic-auth-on-neo4j-2-2-0-rc01
particularly if neo4j installed using homebrew

Services live at 'http://localhost:5000/*', where * is the path defined adding the resource below
'''

app = Flask(__name__)

#Should enable service access despite same-origin policy (cf. JSONP)
cors = CORS(app)
api = Api(app)


#-----------------------------------------------------------------------------------------
#CHEMISTRY
def mol_from_smiles(smiles):
	return Chem.MolFromSmiles(smiles)

def mol_to_2D(mol):
	AllChem.Compute2DCoords(mol)
	return Chem.MolToPDBBlock(mol)

def mol_to_3D(mol):
	mol = Chem.AddHs(mol)
	AllChem.EmbedMolecule(mol)
	AllChem.UFFOptimizeMolecule(mol)
	mol = Chem.RemoveHs(mol)
	return Chem.MolToPDBBlock(mol)


class chemistry_services(Resource): 
	def get(self):
		listing = []
		listing.append('mol_to_2D: convert a SMILES string to a 2D MOL file')
		listing.append('mol_to_3D: convert a SMILES string to a 3D MOL file')
		return json.dumps(listing)


class smiles_2D(Resource): 
	def get(self,smiles):
		try:
			smiles = smiles.decode('utf8')
			mol = mol_from_smiles(smiles)
			mol = mol_to_2D(mol)
			return make_response(mol)
		except:
			return 'Could not convert SMILES'


class smiles_3D(Resource): 
	def get(self,smiles):
		try:
			smiles = smiles.decode('utf8')
			mol = mol_from_smiles(smiles)
			mol = mol_to_3D(mol)
			return make_response(mol)
		except:
			return 'Could not convert SMILES'


class smiles_post(Resource): 
	def post(self,dim):
		try:
			json_data = request.get_json(force=True)
			smiles = json_data['smiles']
			smiles = smiles.decode('utf8')
			print smiles, dim
			mol = mol_from_smiles(smiles)
			if dim == "3D":
				mol = mol_to_3D(mol)
			else:
				mol = mol_to_2D(mol)
			return make_response(mol)
		except:
			return 'Could not convert SMILES'
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
#SBGN

class svg_unknown(Resource): 
	def get(self,label,hex): #Remember security on label (check for / or ..)
		try:
			svg_content = '<?xml version="1.0"?><svg xmlns="http://www.w3.org/2000/svg"><ellipse cx="26" cy="13" rx="25" ry="12" style="fill:#' + hex + ';stroke:black;stroke-width:1"/><text x="10" y="15" font-family="sans-serif" font-size="12px" fill="black">' + label + '</text></svg>'
			response = make_response(svg_content)
			#response.headers['Access-Control-Allow-Origin'] = '*'
			response.mimetype = 'image/svg+xml'
			return response	
		except:
			return 'Error'
			

class svg_reaction(Resource): 
	def get(self,label,hex): #Remember security on label (check for / or ..)
		try:
			svg_content = '<?xml version="1.0"?><svg xmlns="http://www.w3.org/2000/svg"><g transform="translate(15,10) rotate(0)"><rect width="15" height="15" style="fill:#' + hex + ';stroke-width:2;stroke:rgb(0,0,0)" /></g></svg>'
			response = make_response(svg_content)
			#response.headers['Access-Control-Allow-Origin'] = '*'
			response.mimetype = 'image/svg+xml'
			return response	
		except:
			return 'Error'


class svg_simple_chemical(Resource): 
	def get(self,label,hex): #Remember security on label (check for / or ..)
		try:
			svg_content = '<?xml version="1.0"?><svg xmlns="http://www.w3.org/2000/svg"><ellipse cx="42" cy="42" rx="40" ry="40" style="fill:#' + hex + ';stroke:black;stroke-width:2"/><text x="10" y="45" font-family="sans-serif" font-size="16px" fill="black">' + label + '</text></svg>'
			response = make_response(svg_content)
			#response.headers['Access-Control-Allow-Origin'] = '*'
			response.mimetype = 'image/svg+xml'
			return response	
		except:
			return 'Error'
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
#ENDPOINT DEFINITION
api.add_resource(chemistry_services,'/chemistry/listing')
api.add_resource(smiles_2D,'/chemistry/smiles_2D/<smiles>')
api.add_resource(smiles_3D,'/chemistry/smiles_3D/<smiles>')
api.add_resource(smiles_post,'/chemistry/smiles_post/<dim>')
api.add_resource(svg_unknown,'/SBGN/unknown/<label>/<hex>')
api.add_resource(svg_reaction,'/SBGN/reaction/<label>/<hex>')
api.add_resource(svg_simple_chemical,'/SBGN/simple_chemical/<label>/<hex>')
#-----------------------------------------------------------------------------------------

if __name__ == '__main__':
#	Run private
#	app.run(debug=True,threaded=True,port=5001)
#	Run public
	context = ('certs/fullchain.pem', 'certs/privkey.pem')
	app.run(host='0.0.0.0',port=8081,debug=False,ssl_context=context,threaded=True)