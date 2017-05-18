var landingApp = angular.module('landingApp', ['ngTouch','ngSanitize', 'ui.grid', 'ui.grid.selection', 'ui.grid.autoResize']);

landingApp.controller('landingCtrl', ['$scope', '$http', '$rootScope', '$window', '$filter', 'uiGridConstants', function ($scope, $http, $rootScope, $window, $filter, uiGridConstants) {

	//-------------------------------------------------
	//ENV
	$scope.static_url = 'www.metabolicjamboree.co.uk';
	$scope.socket_static_url = 'www.metabolicjamboree.co.uk:8082';
	$scope.external_static_url = 'www.metabolicjamboree.co.uk:8081';

	//For server to push messages to client
	$scope.general_messages = [];

	//Base64 encoding
	function Base64(username,password) {
		var keyStr = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/=';
		input = username + ":" + password;
		var output = "";
		var chr1, chr2, chr3 = "";
		var enc1, enc2, enc3, enc4 = "";
		var i = 0;
		do {
			chr1 = input.charCodeAt(i++);
			chr2 = input.charCodeAt(i++);
			chr3 = input.charCodeAt(i++);
		
			enc1 = chr1 >> 2;
			enc2 = ((chr1 & 3) << 4) | (chr2 >> 4);
			enc3 = ((chr2 & 15) << 2) | (chr3 >> 6);
			enc4 = chr3 & 63;

			if (isNaN(chr2)) {
				enc3 = enc4 = 64;
				}
			else if (isNaN(chr3)) {
				enc4 = 64;
				}
		
			output = output + keyStr.charAt(enc1) + keyStr.charAt(enc2) + keyStr.charAt(enc3) + keyStr.charAt(enc4);
			chr1 = chr2 = chr3 = "";
			enc1 = enc2 = enc3 = enc4 = "";
			} while (i < input.length);
		return 'Basic ' + output;
		};

	//Connect to broadcast server
	console.log('SOCKET Connecting to broadcast server');
	var socket = io.connect('https://' + $scope.socket_static_url + '/mq', {reconnection: true})

	//Connect to message socket	
	$scope.socketId = 'No connection';
	$scope.socket_status = '';
	socket.on('connect', function() {
		console.log('SOCKET Connected to socket:', socket.id);
		console.log('SOCKET', socket);
		$scope.socketId = socket.id;
		$scope.socket = true;
		$scope.socket_status = 'connected';
		$scope.$apply();
		});

	//Socket error
	socket.on('error', function() {
		console.log('Socket error', socket);
		$scope.socket_status = 'error';
		});

	//Socket reconnect
	socket.on('reconnect', function() {
		console.log('Reconnecting socket', socket);
		$scope.socket_status = 'reconnecting';
		});

	//Socket reconnection failure
	socket.on('reconnect_failed', function() {
		console.log('Failed to reconnect socket', socket);
		$scope.socket_status = 'Failed to reconnect';
		});

	$scope.disconnect = function() {
		console.log('Disconnecting');
		socket.disconnect();
		console.log(socket);
		};
	//-------------------------------------------------

	//-------------------------------------------------
	//Actions
	$scope.moleculeActions = [	{"name":"Suggest structure","description":"Use name to structure services to suggest a molecular structure","run":"chemical2structure"},
								{"name":"Find synonyms","description":"Search remote databases for alternative names for this molecule","run":"smallMoleculeSynonyms"},
								{"name":"Cross-reference","description":"Find this molecule in other biological databases","run":"smallMoleculeCrossReference"},
								{"name":"Find literature","description":"Find papers demonstrating this molecule in this organism","run":"smallMoleculeCooccurrence"},
								{"name":"Convert to SMILES","description":"Make InChI into SMILES","run":"smallMoleculeInChIToSmiles"}	];

	$scope.reactionActions = [	{"name":"Auto-balance","description":"Try to tweak the reaction into balance","run":"reactionBalancer"},
								{"name":"Find synonyms","description":"Search remote databases for alternatives names for this reaction","run":"reactionSynonyms"},
								{"name":"Cross-reference","description":"Find this reaction in other biological databases","run":"reactionCrossReference"},
								{"name":"Find literature","description":"Find papers demonstrating this molecule in this organism","run":"smallMoleculeCooccurrence"}	];

	$scope.compartmentActions = [	{"name":"Find synonyms","description":"Search remote databases for alternative names for this compartment","run":"compartmentSynonyms"},
									{"name":"Cross-reference","description":"Find this molecule in other biological databases","run":"smallMoleculeCrossReference"}	];

	//-------------------------------------------------


	//-------------------------------------------------	
	//REGISTRATION
	$scope.sign = true;
	$scope.username = '';
	$scope.password = '';
	$scope.message = '';
	$scope.logged_in = false;
	//For registering a new user
	$scope.passwords = {}
	$scope.passwords.a = '';
	$scope.strength = 'gray';
	$scope.passwords.b = '';
	$scope.agree = false;
	//Socket connection status
	$scope.socket = false;
	//Switch between sign in and registration
	$scope.signToReg = function() {
		$scope.sign = !$scope.sign;
		};
	//Change functions
	$scope.username_update = function(u) {
		$scope.username = u;
		};
	$scope.password_update = function(p) {
		$scope.password = p;
		};
	//Password strength test
	$scope.pwdStrength = function() {
		if($scope.passwords.a.length > 0) {
			response = $http.get('https://' + $scope.static_url + '/pwdStrength/'+$scope.passwords.a)
				.success(function(data) {
					if(data < 0.2) {
						$scope.strength = 'red';
						}
					else if (data < 0.3) {
						$scope.strength = 'orange';
						}
					else {
						$scope.strength ='green';
						}
					})
				.error(function (data, status) {
					console.log('Error',status,data);
					return status;
					});
			}
		else {
			return 0;
			}
		};

	//Register user
	$scope.signUp = function() {
		bundle = {"username":$scope.username, "password":$scope.passwords.a};
		console.log('Registration bundle:', bundle);
		url = 'https://' + $scope.static_url + '/registerUser';
		response = $http.post(url, angular.toJson(bundle) )
			.success(function(data) {
				if(data[0] == true) {
					console.log('Registered');
					$scope.message = 'You just registered!';
					$scope.sign = true;
					}
				else {
					console.log(data[1]);
					$scope.message = data[1];
					}
				})
			.error(function (data, status) {
				console.log('Error',status,data);
				$scope.message = 'Error';
				return false;
				});
		};

	//Sign in user
	$scope.signIn = function() {
		bundle = {"username":$scope.username, "password":$scope.password};
		console.log('Bundle:', $scope.username, $scope.password);
		
		//Generate Basic auth header
		$scope.credentials = Base64($scope.username,$scope.password) ;
		console.log('Base64 credentials:', $scope.credentials);
		$http.defaults.headers.common.Authorization = $scope.credentials;
		
		url = 'https://' + $scope.static_url + '/signIn';
		response = $http.post(url, angular.toJson(bundle) )
			.success(function(data) {
				console.log('Sign in response', data);
				if(data[0] == true) {
					$scope.message = 'Signing in';
					}
				else {
					$scope.message = data[1];
					}
				})
			.error(function (data, status) {
				console.log('Error',status,data);
				$scope.message = 'Error';
				return false;
				})
			.then(function (data) {
					console.log('Then',data.data[0]);
					if(data.data[0] == true) {
						console.log('Look for token')
						port = window.location.search.substr('token');
						token = port.replace('?token=','');
						if(token != '') {
							console.log('Found token', token);
							//Resolve token
							url = 'https://' + $scope.static_url + '/addNewUser'
							response = $http.post(url, angular.toJson({"username":$scope.username,"token":token}))
								.success(function(data) {
									console.log('Resolved token', data);
									})
								.error(function (data, status) {
									console.log('Error', status);
									})
								.then(function(data) {
									console.log(data.data);
									$scope.recon = data.data[1]["recon"];
									$scope.fsid = data.data[1]["fsid"];
									$scope.port = data.data[1]["port"];
									console.log('Token resolution', $scope.recon,$scope.fsid,$scope.port);
									run_app();
									})
							}
						else {
							console.log('No token');
							}
						$scope.logged_in = true;
						$scope.user_recons();
						}
					else {
						$scope.logged_in = false;
						};
				});
		};

	//Sign out
	$scope.signOut = function() {
		console.log('Sign out');
		$scope.sign = true;
		$scope.username = '';
		$scope.password = '';
		$scope.credentials = '';
		$http.defaults.headers.common.Authorization = '';
		$scope.message = '';
		$scope.logged_in = false;
		$scope.passwords = {};
		$scope.agree = false;
		$scope.recon = '';
		$scope.port = '';
		$scope.fsid = '';
		$scope.record = {};
		$scope.selection = '';
		$scope.reconstruction_list = [];
		$scope.label = 'molecule';
		$scope.query = '';
		};
	//-------------------------------------------------	

	//-------------------------------------------------
	//RECONSTRUCTION FUNCTIONS
	initialise = function() {
		console.log('Initialise');
		$scope.spinner = false;				//spinner activity
		$scope.sbml_spinner = false;		//SBML spinner activity
		$scope.sbml_upload_message = '';	//SBML upload message
		$scope.recon = '';					//reconstruction name
		$scope.port = '';					//DB port
		$scope.fsid = '';					//filesystem id
		$scope.label = 'molecule';			//label
		$scope.listLabel = 'molecule';		//listLabel (tracks the list type - see https://github.com/porld/empathyApp/issues/28)
		$scope.colours = ['gray','gray','white'];
		$scope.selection = '';				//Short handle to current selection
		$scope.record = {};					//Full object for current selection
		$scope.path = [];
		$scope.query = '';					//Query
		jsonLists = ['is','isDescribedBy','isVersionOf','property',"is_a","references"]; //This is for processing lists made of JSONified objects (see unpacking updates within $watch on $scope.selection)
		//$scope.reconstruction_list = []; 	//List reconstructions owned by user
		
		//Panel 1 grid
		$scope.gridOptions = {
			enableRowSelection: true,
			enableRowHeaderSelection: false,
			exporterMenuCsv: false,
			enableColumnMenus: false,
			enableGridMenu: false,
			gridMenuShowHideColumns: false,
			enableSelectAll: false,
			selectionRowHeaderWidth: 35,
			rowHeight: 35,
			showGridFooter:false,
			multiSelect:false,
			modifierKeysToMultiSelect: false,
			noUnselect: false,
			enableFiltering: false,
			enableHiding: false,
			enableSorting: false,
			showHeader: false,
			enableVerticalScrollbar: 0
			};

		//Panel 1 column definitions
		$scope.gridOptions.columnDefs = [
			{ name: 'id', visible: false },
			{ name: 'name', displayName:'Name'},
			{ name: 'tags', visible: false} //displayName:'Tags', cellTemplate: '<span style="line-height: 35px; height:100%; overflow-x:scroll !important;" ng-repeat="tag in row.entity.tags"><small><a href="#" ng-click="grid.appScope.filterData(tag)" style="margin: 2px 2px 2px 2px; background: #e4e4e4; color: #666666" class="round label">{{tag}}</a></small></span>'}
			];

		//Register grid API
		$scope.gridOptions.onRegisterApi = function(gridApi){
			console.log('Register grid API');
			//set gridApi on scope
			$scope.gridApi = gridApi;
			gridApi.selection.on.rowSelectionChanged($scope,function(row){
				console.log('Selected',row.entity.id, row.entity.name);
				$scope.label = $scope.listLabel; //Make sure our switches are applied when moving between records (https://github.com/porld/empathyApp/issues/28)
				$scope.selection = row.entity.id;
				//$scope.fetchSelection(row.entity.id); //No longer required because of $watch on selection
				});
			};
		console.log('initialised');
		}
	initialise();

	//Create new reconstruction
	$scope.launch_docker = function(recon_name,notes) {
		$scope.recon_spinner = true;
		initialise();
		console.log('Launch Docker', $scope.recon_spinner);
		url = 'https://' + $scope.static_url + '/docker'
		$http.post(url, angular.toJson({"username":$scope.username,"recon_name":recon_name}) )
			.success(function(data) {
				console.log('Triggered: launch Docker', data);
				$scope.port = data[0];
				$scope.fsid = data[1];
				})
			.error(function (data, status) {
				console.log('Error', status, data);
				})
			.then(function(data, status) {
				console.log('Port:', $scope.port);
				console.log('Check live and activate');
				url = 'https://' + $scope.static_url + '/checkLive';
				$http.post(url,angular.toJson({"username":$scope.username,"port":$scope.port,"recon_name":recon_name,"notes":notes}))
					.success(function(data) {
						console.log('Response', data);
						})
					.error(function (data, status) {
						console.log('Error', status, data);
						})
					.then(function() {
						console.log('Update list of user reconstructions');
						$scope.user_recons();
						$scope.recon_spinner = false;
						});
					$scope.spinner = false;
					});
			};
	
	//Load generic cellular organisation
	$scope.load_cell = function(cell) {
		$scope.recon_spinner = true;
		console.log('Add generic cellular organisation', cell);
		message_handle = $scope.port + '_compartment';
		url = 'https://' + $scope.static_url + '/subcell'
		$http.post(url, angular.toJson({"port":$scope.port,"message_handle":message_handle,"cell":cell}) )
			.success(function(data) {
				console.log('Triggered: load_cell', data);
				})
			.error(function (data, status) {
				console.log('Error', status, data);
				})
			.then(function(data, status) {
				console.log('Created cellular organisation:', cell, 'on', message_handle);
				$scope.recon_spinner = false;
				});
		};

	//Get existing reconstructions
	$scope.user_recons = function() {
		console.log('Fetch user reconstructions');
		url = 'https://' + $scope.static_url + '/userRecons'
		$http.post(url, angular.toJson({"username":$scope.username}) )
			.success(function(data) {
				$scope.reconstruction_list = data[1];
				console.log('Reconstruction list:', $scope.reconstruction_list);
				})
			.error(function (data, status) {
				console.log('Error', status, data);
				});
		};
	
	//Load a reconstruction
	$scope.load_reconstruction = function(reconstruction) {
		console.log('Loading reconstruction', reconstruction, $scope.label);
		$scope.label = 'molecule';
		initialise();
		updateListSocket();
		$scope.recon = reconstruction.name;
		$scope.port = reconstruction.port;
		$scope.fsid = reconstruction.fsid;
		$scope.general_messages = [];

		//Attach general message socket
		socket.on($scope.port + '_' + $scope.credentials, function(message) {
			message = angular.fromJson(message);
			console.log('General message', message);
			message.read = false;
			$scope.unseen = true;
			$scope.general_messages.push(message);
			console.log('General messages', $scope.general_messages);
			});

		$scope.recon_spinner = false;
		run_app();
		};

	//Fetch reconstruction from SYNBIOCHEM-DB
	$scope.synbiochem_import = function(ncbi) {
		console.log('Importing from SYNBIOCHEM-DB', ncbi);
		$scope.recon_spinner = true;
		url = 'https://' + $scope.static_url + '/fetchSYNBIOCHEM'
		$http.post(url, angular.toJson({"port":$scope.port,"ncbi":ncbi}) )
			.success(function(data) {
				console.log('Initiated');
				})
			.error(function (data, status) {
				console.log('Error', status, data);
				})
			.then(function () {
				console.log('Running');
				$scope.recon_spinner = false;			
				});
		};

	//Import from SBML file
	function SBMLparser(sbml) {
		$scope.sbml_upload_message = 'uploading...';
		console.log('Parsing SBML...');
		url = 'https://' + $scope.static_url + '/importSBML'
		$http.post(url, angular.toJson({"port":$scope.port,"sbml":sbml}) )
			.success(function(data) {
				console.log('Pushing SBML to server');
				})
			.error(function (data, status) {
				console.log('Error', status, data);
				})
			.then(function () {
				console.log('Converting on server');
				$scope.sbml_spinner = false;
				$scope.sbml_upload_message = '';
				});
		};
	
	//Filereader to fetch SBML string
	$scope.importSBML = function() {
		$scope.sbml_spinner = true;

		//Create socket for SBML upload messages
		console.log('SOCKET Create SBML callback',socket);
		$scope.sbml_message = 'reading...';
		$scope.sbml_handle = $scope.port + '_sbml';
		socket.on($scope.sbml_handle, function(sbml_message) {
			console.log('SBML message:', sbml_message);
			$scope.sbml_message = sbml_message;
			$scope.$apply();
			});

		console.log('SBML spinner:', $scope.sbml_spinner);
		var file = document.getElementById("SBMLinput").files[0];
		console.log('File:', file);
		if(file) {
			console.log('Importing from SBML');
			var sbmlReader = new FileReader();
			sbmlReader.readAsText(file, "UTF-8");
			//Load
			sbmlReader.onload = function (evt) {
                document.getElementById("SBMLinput").innerHTML = evt.target.result;
                $scope.sbmlString = sbmlReader.result;
                $scope.fileName = document.getElementById("SBMLinput").files[0].name;
                $scope.fileSize = document.getElementById("SBMLinput").files[0].size;
	            console.log('SBML filename:', $scope.fileName, $scope.fileSize);
	            SBMLparser($scope.sbmlString);
				}
			//Progress bar
			sbmlReader.onprogress = function(data) {
				if (data.lengthComputable) {                                            
					var progress = parseInt( ((data.loaded / data.total) * 100), 10 );
					$scope.sbml_message = 'Reading: ' + progress + '%';
					}
				}
			//Error
			sbmlReader.onerror = function (evt) {
				console.log('SBML read error');
				$scope.sbml_spinner = false;
        		}
			}
		};

	//Add collaborator to reconstruction
	$scope.collaborator_update = function(collaborator) {
		$scope.collaborator = collaborator;
		}; 
	$scope.invite = function() {
		$scope.invitation = 'Sending invitation';
		console.log('Inviting', $scope.collaborator);
		bundle = {"username": $scope.username, "recon": $scope.recon, "fsid": $scope.fsid, "port": $scope.port, "user_email": $scope.collaborator };		
		url = 'https://' + $scope.static_url + '/msgNewUser'
		$http.post(url, angular.toJson(bundle) )
			.success(function(data) {
				console.log('Invitation sent');
				$scope.invitation = 'Invitation sent';
				})
			.error(function (data, status) {
				console.log('Error: invitation not sent', status, data);
				$scope.invitation = 'Invitation sent';
				});
		$scope.invitation = '';
		};
	//-------------------------------------------------	

	//-------------------------------------------------	
	//Socket handlers		
	//Pick up broadcasts about list
	function updateListSocket(new_label) {
		console.log('updateListSocket. Switching from', $scope.label, 'to', new_label);
		//Create new callback
		console.log('SOCKET Create new callback',socket);
		$scope.message_handle = $scope.port + '_' + new_label;
		socket.on($scope.message_handle, function(latest_list) {
			console.log('SOCKET', socket);
			//console.log('List update', $scope.message_handle, '>>>', latest_list);
			$scope.gridOptions.data = latest_list;
			$scope.grid_data = latest_list;
			$scope.count = latest_list.length;
			$scope.$apply();
			$scope.spinner = false;
			});

		//Remove old callback
		console.log('SOCKET remove old callback', socket);
		socket.off($scope.port + '_' + $scope.label);

		//Update label
		$scope.label = new_label;		//Change label
		};
	
	//Note the record socket is attached by the $watch on $scope.selection	
	//-------------------------------------------------	

	//-------------------------------------------------	
	//In-app functions	
	//Fetch list corresponding to $scope.list_message_handle ($scope.port_$scope.label)
	$scope.initialiseList = function() {
		$scope.spinner = true;
		$scope.list_message_handle = $scope.port + '_' + $scope.label;
		console.log('Fetch list', $scope.list_message_handle);
		bundle = {"label":$scope.label, "message_handle": $scope.list_message_handle, "port": $scope.port};
		url = 'https://' + $scope.static_url + '/listNode'
		console.log('Hitting url:', url, 'with message handle', $scope.list_message_handle);
		$http.post(url, angular.toJson(bundle) )
			.success(function(data) {
				console.log('Triggered: initialise list');
				})
			.error(function (data, status) {
				console.log('Error',status,data);
				})
			.then(function (data) {
				//console.log('initialiseList RESTful response', data.data);
				$scope.gridOptions.data = data.data; //Commented out because this is the live-filtered list
				$scope.grid_data = data.data; 
				$scope.count = data.data.length;
				$scope.spinner = false;
				});
		};

	//Create a new node of type $scope.label
	$scope.createNode = function() {
		$scope.spinner = true;
		//createNode($scope.label,$scope.username,$scope.password,$scope.port,$scope.static_url,$http);
		console.log('Create new node');
		message_handle = $scope.port + '_' + $scope.label;
		bundle = {"label":$scope.label, "message_handle": message_handle, "port":$scope.port};
		url = 'https://' + $scope.static_url + '/createNode'
		$http.post(url, angular.toJson(bundle) )
			.success(function(data) {
				console.log('Triggered: create node', data);
				$scope.selection = data;
				//$scope.fetchSelection(data);
				})
			.error(function (data, status) {
				console.log('Error',status,data);
				});
		fetchCompartmentList();
		$scope.spinner = false;
		};

	//Destroy node from $scope.selection
	$scope.destroyNode = function() {
		$scope.spinner = true;
		console.log('Destroy', $scope.selection);
		record_handle = $scope.port + '_' + $scope.selection;
		message_handle = $scope.port + '_' + $scope.label;
		bundle = {"target":$scope.selection, "label":$scope.label, "record_handle": record_handle, "message_handle": message_handle, "port":$scope.port};
		url = 'https://' + $scope.static_url + '/destroyNode';
		$http.post(url, angular.toJson(bundle) )
			.success(function(data) {
				console.log('Triggered: destroy node', data);
				})
			.error(function (data, status) {
				console.log('Error',status,data);
				})
			.then(function() {
				console.log('Destroyed', $scope.selection);
				//Detach socket callback
				socket.off(record_handle);
				$scope.selection = '';
				$scope.record = '';
				$scope.spinner = false;
				});
		};

	//Destroy edge (take molecular species out of reaction)
	$scope.destroyEdge = function(rxn_id,part_id) {
		console.log('Remove',part_id,'from',rxn_id);
		$scope.spinner = true;
		record_handle = $scope.port + '_' + rxn_id;
		bundle = {"targetA":rxn_id, "targetB":part_id, "label":$scope.label, "record_handle": record_handle, "port":$scope.port};
		url = 'https://' + $scope.static_url + '/destroyEdge';
		$http.post(url, angular.toJson(bundle) )
			.success(function(data) {
				console.log('Triggered: removing', part_id, 'from', rxn_id);
				})
			.error(function (data, status) {
				console.log('Error',status,data);
				})
			.then(function() {
				console.log('Removed');
				$scope.spinner = false;
				});
		};

	//Change label
	$scope.change_label = function(label) {
		console.log('Change label');
		colours = ['gray','gray','gray'];
		if(label == 'molecule') {
			$scope.colours = ['gray','gray','white'];
			updateListSocket(label);		//Update broadcast handle
			$scope.listLabel = label;
			$scope.initialiseList();			//Fetch list
			}
		else if(label == 'reaction') {
			$scope.colours = ['gray','white','gray'];
			$scope.listLabel = label;
			updateListSocket(label);		//Update broadcast handle
			$scope.initialiseList(label);
			}
		else {
			$scope.colours = ['white','gray','gray'];
			$scope.listLabel = label;
			updateListSocket(label);		//Update broadcast handle
			$scope.initialiseList();
			};
		$scope.selection = '';
		console.log('LABEL AND LISTLABEL:', $scope.label,$scope.listLabel);
		return colours;
		};
	
	//Fetch the node by REST
	function fetchNode(selection) {
		if(selection != '') {
			console.log('Fetching node', selection);
			bundle = {"port":$scope.port, "selection":selection, "label":$scope.label};
			//console.log('SURELY THIS DOES NOT WORK!?');
			//Reset a few selections to enable selection of same value as previous
			$scope.typeChoice = '';
			$scope.compChoice = '';
			url = 'https://' + $scope.static_url + '/fetchSelection';
			$http.post(url, angular.toJson(bundle) )
				.success(function(data) {
					console.log('Triggered: fetch node (success)');
					})
				.error(function (data, status) {
					console.log('Error',status,data);
					$scope.record = {};
					})
				.then(function (data) {
					console.log('Record returned (then)',data);
					$scope.record = data.data;
					console.log('RECORD',$scope.record);
					//Reset tickers and associates
					$scope.notTicker = 1;
					$scope.refTicker = 1;
					$scope.notification = '';
					$scope.reference = '';

					//Unpack JSONified lists
					$scope.record.is = unpackJson($scope.record.is);
					$scope.record.hasPart = unpackJson($scope.record.hasPart);
					$scope.record.isPartOf = unpackJson($scope.record.isPartOf);
					$scope.record.isVersionOf = unpackJson($scope.record.isVersionOf);
					$scope.record.hasVersion = unpackJson($scope.record.hasVersion);
					$scope.record.isHomologTo = unpackJson($scope.record.isHomologTo);
					$scope.record.isDescribedBy = unpackJson($scope.record.isDescribedBy);
					$scope.record.isEncodedBy = unpackJson($scope.record.isEncodedBy);
					$scope.record.encodes = unpackJson($scope.record.encodes);
					$scope.record.occursIn = unpackJson($scope.record.occursIn);
					$scope.record.hasProperty = unpackJson($scope.record.hasProperty);
					$scope.record.isPropertyOf = unpackJson($scope.record.isPropertyOf);
					$scope.record.hasTaxon = unpackJson($scope.record.hasTaxon);
					//$scope.record.is_a = unpackJson($scope.record.is_a);
					//$scope.record.references = unpackJson($scope.record.references);

					//Unpack notes
					$scope.record.notes = angular.fromJson($scope.record.notes);
					
					console.log('Unpacked record',$scope.record);

					//$scope.record.notifications = unpackJson($scope.record.notifications);					
					//Trigger ancillary functions
					if($scope.label == 'molecule') { //Compartments and molecules have parent compartments
						fetchCompartmentList();
						}
					//smilesFromRecord(); //Disabled as smilesViewer takes care of this now
					//$scope.reference = $scope.record.references[0];
					//$scope.notification = $scope.record.notifications[0];
					$scope.spinner = false;
					});
			}
		else {
			$scope.record = {};
			}
		};

	//$watch for changes to selection to track in $scope.path
	console.log('Attaching $watch to selection');
	$scope.$watch('selection', function(newValue, oldValue) {
		if($scope.logged_in) { //Only if we're logged in
			$scope.spinner = true;
			console.log('Selection just changed:', oldValue, newValue);

			//Remove old record callback
			console.log('SOCKET remove old record callbacks', $scope.port+'_'+oldValue);
			socket.off($scope.port + '_' + oldValue);
			console.log('SOCKET removed old record callback', socket);

			//Create new callback (also includes handler for processing messages)
			$scope.record_handle = $scope.port + '_' + newValue;
			console.log('Create new record callback', newValue, socket);
			socket.on($scope.record_handle, function(record_update) {
				console.log('SOCKET hit socket handle:', $scope.record_handle, socket);
				console.log('SOCKET record_update content:', record_update);
				id = record_update['id'];
				key = record_update['key'];
				value = record_update['value'];
				//Look for blank updates (this means refresh the whole record)
				if(key === '') {
					console.log('Blank key means refresh whole record', id);
					//RESTful trigger an update to the whole of $scope.record
					fetchNode(id);
					}
				//Non-blank updates involving poking things into the current record
				else {
					//Poke into $scope.record
					console.log('RECORD UPDATE\tid:',id,'key:',key,'value:',value);
					if(jsonLists.indexOf(key) > -1) { //jsonLists holds the names of lists of JSONified objects
						//Unpack JSONified lists
						console.log('SOCKET: this is a JSON list', key);
						value = unpackJson(value);
						};

					//Update field
					console.log('SOCKET: poke key', key, 'into record', $scope.record, 'with new value', value);
					$scope.record[key] = value;
					$scope.$apply();
			
					//Update grid if name or tags have changed
					if( (key == 'name') || (key == 'tags') ) {
						console.log('Update grid');
						$scope.initialiseList();
						};
					
					}
				$scope.spinner = false;
				});	

			if($scope.label == 'molecule') {
				//Fetch latest compartment list
				fetchCompartmentList();
				console.log('Fetched compartments:', $scope.compartmentList);
				}
			
			//Fetch record on selection
			fetchNode(newValue);
	
			//Path handling
			$scope.path.push(newValue);
			$scope.path = $filter('limitTo')($scope.path, -2, 0);
			console.log('Path:',$scope.path);
			}
		});
	
	//-------------------------------------------------	
	//ACTIONS UPON TEXT FIELDS
	//Update node property
	$scope.updateText = function(key,value) {
		console.log('Update text', $scope.selection, 'with', key, ':', value);
		if(key === 'notes') {		//Wrap notes in JSON (unwrapper expects this)
			console.log('Update text wrapping to JSON');
			value = angular.toJson(value);
			}
		bundle = {"port": $scope.port, "record_handle": $scope.record_handle, "id":$scope.selection, "key":key, "value":value};
		url = 'https://' + $scope.static_url + '/updateText';
		$http.post(url, angular.toJson(bundle) )
			.success(function(data) {
				console.log('Triggered: update node property', data);
				})
			.error(function (data, status) {
				console.log('Error',status,data);
				});
		};
	//-------------------------------------------------	

	//-------------------------------------------------	
	//ACTIONS UPON LISTS
	//Add element to list
	$scope.listPush = function(key,value) {
		if(typeof value === 'undefined') {
			console.log('No value defined');
			}
		else if(value.length == 0) {
			console.log('Zero length element', value, value.length);
			}
		else if($scope.record[key].includes(value) ) {
			console.log('Duplicate element', value, 'in', $scope.record[key]);
			}
		else {
			console.log('Push', value, 'to', key, 'in', $scope.selection);
			newList = $scope.record[key];
			newList.push(value);
			bundle = {"port": $scope.port, "record_handle": $scope.record_handle, "id":$scope.selection, "key":key, "value":newList};
			url = 'https://' + $scope.static_url + '/listPush';
			$http.post(url, angular.toJson(bundle) )
				.success(function(data) {
					console.log('Triggered: list push');
					})
				.error(function (data, status) {
					console.log('Error',status,data);
					});
			}
		};

	//Change curation status of subfield in JSON list of type field (restrict to value in case of subfields with same names)
	$scope.curateList = function(field,value) {
		console.log('Curation swap', value, 'in', field, 'of', $scope.selection);
		newList = [];
		flag = true;
		for (i = 0; i < $scope.record[field].length; i++) {
			currentValue = $scope.record[field][i];
			//Check if we're looking at the right value
			if(currentValue === value) {
				if(currentValue[0] == '?') {
					currentValueEnd = currentValue.length - 1;
					currentValue = currentValue.substr(-1*currentValueEnd);
					}
				else{
					currentValue = '?' + currentValue;
					}
				flag = false;
				newList.push( currentValue );
				}
			else {
				newList.push( currentValue );
				}
			};
		if(flag) { //If the property isn't there create a new one
			newList.push( currentValue );
			}
		bundle = {"port": $scope.port, "record_handle": $scope.record_handle, "id":$scope.selection, "key":field, "value":newList};
		url = 'https://' + $scope.static_url + '/listPush';
		$http.post(url, angular.toJson(bundle) )
			.success(function(data) {
				console.log('Triggered: update json text', data);
				})
			.error(function (data, status) {
				console.log('Error',status,data);
				});
		};
	
	//Remove element from list
	$scope.listPop = function(key,value) {
		console.log('Pop', value, 'from', key, 'in', $scope.selection);
		currentList = $scope.record[key];
		//console.log('Current list:', currentList);
		newList = [];
		for (i = 0; i < currentList.length; i++) {
			if(currentList[i] != value) {
				newList.push(currentList[i]);
				}
			};
		//console.log('New list:', newList);
		bundle = {"port": $scope.port, "record_handle": $scope.record_handle, "id":$scope.selection, "key":key, "value":newList};
		url = 'https://' + $scope.static_url + '/listPop';
		$http.post(url, angular.toJson(bundle) )
			.success(function(data) {
				console.log('Triggered: list pop');
				})
			.error(function (data, status) {
				console.log('Error',status,data);
				});
		};
	//-------------------------------------------------	

	//-------------------------------------------------	
	//ACTIONS UPON JSON LISTS
	//Add element to JSON list
	$scope.jsonListPush = function(field,value) {
		source = value["source"];
		id = value["id"];
		console.log('JSON push', field, value, source, id);
		if(typeof id === 'undefined') {
			console.log('No value defined');
			}
		else if(source.length == 0) {
			console.log('Zero length element (source)', source);
			}
		else if(id.length == 0) {
			console.log('Zero length element (id)', id);
			}
		else if($scope.record[field].indexOf([source,id]) > -1) {
			console.log('Duplicate element', [source,id], 'in', $scope.record[field]);
			}
		else {
			console.log('Push', [source,id], 'to', field, 'in', $scope.selection);
			newList = [];
			for (i = 0; i < $scope.record[field].length; i++) {
				newList.push( angular.toJson($scope.record[field][i]) );
				};
			newList.push( angular.toJson([source,id]) );
			bundle = {"port": $scope.port, "record_handle": $scope.record_handle, "id":$scope.selection, "key":field, "value":newList};
			url = 'https://' + $scope.static_url + '/listPush';
			$http.post(url, angular.toJson(bundle) )
				.success(function(data) {
					console.log('Triggered: json list push', data);
					})
				.error(function (data, status) {
					console.log('Error',status,data);
					});
			}
		};

	//Change curation status of subfield in JSON list of type field (restrict to value in case of subfields with same names)
	$scope.curateJsonText = function(field,subfield,value) {
		console.log('Curation swap', subfield, '(', value, ')', 'in', field, 'of', $scope.selection);
		newList = [];
		//JSONify each subfield in field
		flag = true;
		for (i = 0; i < $scope.record[field].length; i++) {
			currentPair = $scope.record[field][i];
			//Check if we're looking at the right subfield
			if((currentPair[0] === subfield) && (currentPair[1] === value)) {
				newValue = currentPair[1];
				if(newValue[0] == '?') {
					newValueEnd = newValue.length - 1;
					newValue = newValue.substr(-1*newValueEnd);
					}
				else{
					newValue = '?' + newValue;
					}
				flag = false;
				newList.push( angular.toJson([currentPair[0],newValue]) );
				}
			else {
				newList.push( angular.toJson(currentPair) );
				}
			};
		if(flag) { //If the property isn't there create a new one
			newList.push( angular.toJson([subfield,value]) );
			}
		bundle = {"port": $scope.port, "record_handle": $scope.record_handle, "id":$scope.selection, "key":field, "value":newList};
		url = 'https://' + $scope.static_url + '/listPush';
		$http.post(url, angular.toJson(bundle) )
			.success(function(data) {
				console.log('Triggered: update json text', data);
				})
			.error(function (data, status) {
				console.log('Error',status,data);
				});
		};

	//Update element in JSON list
	$scope.updateJsonText = function(field,subfield,value) {
		console.log('JSON push', field, subfield, value);
		if(typeof value === 'undefined') {
			console.log('No value defined');
			}
		else if(value.length == 0) {
			console.log('Zero length element (source)', source);
			}
		else {
			console.log('Edit', [subfield,value], 'in', field, 'of', $scope.selection);
			newList = [];
			//JSONify each subfield in field
			flag = true;
			for (i = 0; i < $scope.record[field].length; i++) {
				currentPair = $scope.record[field][i];
				//Check if we're looking at the right subfield
				console.log('Current pair:',currentPair);
				if(currentPair[0] === subfield) {
					currentPair[1] = value;
					flag = false;
					}
				newList.push( angular.toJson(currentPair) );
				};
			if(flag) { //If the property isn't there
				newList.push( angular.toJson([subfield,value]) );
				}
			bundle = {"port": $scope.port, "record_handle": $scope.record_handle, "id":$scope.selection, "key":field, "value":newList};
			url = 'https://' + $scope.static_url + '/listPush';
			$http.post(url, angular.toJson(bundle) )
				.success(function(data) {
					console.log('Triggered: update json text', data);
					})
				.error(function (data, status) {
					console.log('Error',status,data);
					});
			}
		};

	//Remove element from list of Json entries
	$scope.jsonListPop = function(field,jsonObject) {
		console.log('Json pop', jsonObject, 'from', field, 'in', $scope.selection);
		currentList = $scope.record[field];
		newList = [];
		for (i = 0; i < currentList.length; i++) {
			//console.log(currentList[i]);
			if(currentList[i] != jsonObject) {
				newList.push(angular.toJson(currentList[i])); //JSONify so it can go into the database (we go back to the client via the socket)
				}
			};
		//console.log('New list:', newList);
		bundle = {"port": $scope.port, "record_handle": $scope.record_handle, "id":$scope.selection, "key":field, "value":newList};
		url = 'https://' + $scope.static_url + '/listPop';
		$http.post(url, angular.toJson(bundle) )
			.success(function(data) {
				console.log('Triggered: json list pop', data);
				})
			.error(function (data, status) {
				console.log('Error',status,data);
				});
		};
	//-------------------------------------------------		
	
	//Render record
	$scope.renderRecord = function(record) {
		if($scope.label == 'molecule') {
			$scope.render = renderMolecule(record);
			}
		else if($scope.label == 'reaction') {
			$scope.render = renderReaction(record);
			}
		else if($scope.label == 'compartment') {
			$scope.render = renderCompartment(record);
			}
		else {
			$scope.render = 'What was that and where did you get it from!?';
			};
		};
	//-------------------------------------------------	

	//-------------------------------------------------
	//Filter
	//Grid filter function for searching
	$scope.filterData = function(query) {
		console.log('Filtering...', query);
		query = query.toLowerCase();
  		temp = []
  		for (i = 0; i < $scope.grid_data.length; i++) {
  			//Search against name string
		  	rowString = $scope.grid_data[i]["name"];
		  	rowString = rowString.toLowerCase();
		  	caughtFlag = false
			if (rowString.includes(query)) {
		    	temp.push($scope.grid_data[i]);
		    	caughtFlag = true;
		    	}
			if(!caughtFlag) {
				//Search against tags
				tagList = $scope.grid_data[i]["tags"];
				for (j = 0; j < tagList.length; j++) {
					tag = tagList[j];
					tag = tag.toLowerCase();
					if (tag.includes(query)) {
						temp.push($scope.grid_data[i]);
						}
					};
				}
			};
		$scope.gridOptions.data = temp;
		};
	//-------------------------------------------------

	//-------------------------------------------------
	//Application
	run_app = function() {
		console.log('Run app...');
		$scope.spinner = false;
		$scope.recon_spinner = false;
		$scope.change_label('molecule');

		};
	//End run_app
	//-------------------------------------------------

	//-------------------------------------------------
	//Ancillary functions
	
	//Change selection from within ng-repeats (for hopping between reactions and molecular species, for example)
	$scope.makeSelection = function(id,label) {
		$scope.selection = id;
		$scope.label = label;
		};
	
	//Molecular species types
	$scope.molecular_species_types = ['simple chemical','macromolecule','complex','other'];
	$scope.type_flicker = function(flick) {
		if(flick < ($scope.molecular_species_types.length - 1) ) {
			flick = flick + 1;
			}
		else {
			flick = 0;
			}
		$scope.typeChoice = $scope.molecular_species_types[flick];
		return flick;
		};
	
	//Look up compartment name on record load (loaded when fetchCompartmentList fulfils promise)
	function compartmentNameFromList() {
		//console.log('Looking for', $scope.record.inCompartment, 'in', $scope.compartmentList);
		$scope.compName = '';
		for (k = 0; k < $scope.compartmentList.length; k++) {
			if ($scope.compartmentList[k]["id"] === $scope.record.inCompartment) {
				$scope.compName = $scope.compartmentList[k]["name"];
				}
			};
		if($scope.compName === '') {
			$scope.record.inCompartment = '';
			$scope.compName = 'Unknown';
			}
		};

	//$watch for changes to inCompartment
	$scope.$watch('record.inCompartment', function(newValue, oldValue) {
		//console.log('$watch picked up change:', oldValue, '>', newValue);
		compartmentNameFromList();
		});

	//Fetch compartment list (for specifying molecules, parent compartments, etc.). Note this is a direct RESTful load.
	$scope.compartmentList = [];
	function fetchCompartmentList() {
		console.log('Fetch compartment list');
		bundle = {"label":"compartment", "message_handle": "", "port": $scope.port};
		url = 'https://' + $scope.static_url + '/listNode';
		console.log('Hitting url:', url, 'without message handle');
		$http.post(url, angular.toJson(bundle) )
			.success(function(data) {
				console.log('Triggered: fetch compartment list');
				})
			.error(function (data, status) {
				console.log('Error',status,data);
				})
			.then(function(data) {
				console.log('Fetched compartment list', data.data);
				$scope.compartmentList = data.data;
				})
		};
		
	//Pull smiles (DO NOT USE - SEE valueFromJsonList)
	$scope.smiles = '';
	function smilesFromRecord() {
		console.log('Pull SMILES from record');
		$scope.smiles = '';
		for (i = 0; i < $scope.record.isDescribedBy.length; i++) {
			entry = $scope.record.isDescribedBy[i];
			if(entry[0] === 'smiles') {
				$scope.smiles = entry[1];
				}
			};
		};

	//Pull value from Json list
	function valueFromJsonList(field,subfield) {
		console.log('Pull (valueFromJsonList)', subfield, 'from', field);
		for (i = 0; i < $scope.record[field].length; i++) {
			entry = $scope.record[field][i];
			if(entry[0] === subfield) {
				
				return entry[1];
				}
			};
		};
		
	//Split key:value Json string into object
	function unpackJson(listObject) {
		//Unpack Jsonified properties
		unpackedList = [];
		for (i = 0; i < listObject.length; i++) {
			unpackedList.push(angular.fromJson(listObject[i]));
			};
		return unpackedList;
		};
		
	//Molecule search
	$scope.searchMolecules = function(query) {
		if(query.length > 2) {
			console.log('Query:', query);
			$scope.queryMatches = [];
			bundle = {"query":query,"port":$scope.port};
			url = 'https://' + $scope.static_url + '/queryMolecules';
			$http.post(url, angular.toJson(bundle) )
				.success(function(data) {
					console.log('Triggered: search molecules', data);
					})
				.error(function (data, status) {
					console.log('Error', status, data);
					})
				.then(function(data, status) {
					$scope.queryMatches = data.data;
					});
			}
		};
		
	//Pick up reaction participation type
	$scope.getRole = function(role) {
		$scope.participant = role;
		};
	
	//Reference ticker
	$scope.refTicker = 1;
	$scope.refTick = function(way,max) {
		console.log('refTick', $scope.refTicker);
		if(way == 'up') {
			if($scope.refTicker == max ) {
				$scope.refTicker = 1;
				}
			else {
				$scope.refTicker = $scope.refTicker + 1;
				}
			$scope.reference = $scope.record.references[$scope.refTicker-1];
			}
		else if(way == 'down') {
			if($scope.refTicker == 1) {
				$scope.refTicker = max;
				}
			else {
				$scope.refTicker = $scope.refTicker - 1;
				}
			$scope.reference = $scope.record.references[$scope.refTicker-1];
			}
		else {
			$scope.reference =	$scope.record.references[0];
			}
		console.log('Reference', $scope.reference);
		};

	//Notification ticker
	$scope.notTicker = 1;
	$scope.notTick = function(way,max) {
		console.log('notTick', $scope.notTicker);
		if(way == 'up') {
			if($scope.notTicker == max ) {
				$scope.notTicker = 1;
				}
			else {
				$scope.notTicker = $scope.notTicker + 1;
				}
			$scope.notification = $scope.record.notifications[$scope.notTicker-1];
			}
		else if(way == 'down') {
			if($scope.notTicker == 1) {
				$scope.notTicker = max;
				}
			else {
				$scope.notTicker = $scope.notTicker - 1;
				}
			$scope.notification = $scope.record.notifications[$scope.notTicker-1];
			}
		else {
			$scope.notification = $scope.record.notifications[0];
			}
		console.log('Notification', $scope.notification);
		};

	//REACTION BUILDER
	$scope.reactionBuilder = function(participant) {
		var naughty = true;
		if($scope.participant == 'reactant') {
			hasRole = 'hasReactant';
			}
		else if($scope.participant == 'modifier') {
			hasRole = 'hasModifier';
			}
		else if($scope.participant == 'product') {
			hasRole = 'hasProduct';
			}
		else {
			naughty = false;
			}
			
		if(naughty) {
			console.log('Add', $scope.participant, 'with id', participant, 'to', $scope.selection);
			record_handle = $scope.port + '_' + $scope.selection;
			bundle = {"port":$scope.port,"record_handle":record_handle,"reaction":$scope.selection,"molecule":participant,"role":hasRole};
			url = 'https://' + $scope.static_url + '/makeConnection';
			$http.post(url, angular.toJson(bundle) )
				.success(function(data) {
					console.log('Making connection', data);
					})
				.error(function (data, status) {
					console.log('Error', status, data);
					})
				.then(function(data, status) {
					console.log('Connection made');
					});
			}
		};
	//-------------------------------------------------


	//-------------------------------------------------
	//Actions
	$scope.moleculeAction = function(action) {
		console.log('Run molecule action:', action);
		url = 'https://' + $scope.static_url + '/actionMolecule'
		$http.post(url, angular.toJson({"port":$scope.port,"record":$scope.record,"action":action.run, "credentials":$scope.credentials, "record_handle":$scope.port+'_'+$scope.selection}) )
			.success(function(data) {
				console.log('Triggered:', action.run);
				})
			.error(function (data, status) {
				console.log('Error', status, data);
				});
		};

	//-------------------------------------------------

	$scope.$watch('selection', function(newValue, oldValue) {
		console.log('SOCKET NOW',socket);
		});

	//End controller
	}]);


//----------------------------------------------------------------------------------------
//DIRECTIVES

//SMILES
landingApp.directive('smilesViewer', function ($parse, $http) {
    var molObject = {
		restrict: 'E',
		replace: true,
		scope: {record: '=',
				dim: '=',
				rotate: '=',
				ngid: '=',
				external: '='},
		link: function (scope, element, attrs) {
			scope.$watch('record', function(record) {
				//Look for SMILES in record
				if(record.is) {
					console.log('record.is:', scope.record.is);
					smiles = '';
					for (i = 0; i < scope.record.is.length; i++) {
						if( scope.record.is[i][0] == 'smiles' ) {
							smiles = scope.record.is[i][1];
							console.log('smilesViewer found:', smiles);
							}
						};

					console.log('PreSMILES:', smiles);
					console.log('End SMILES:', smiles.substr(smiles.length-1,1) );
					if( smiles[0] == '?') {
						smiles = smiles.substr(1,smiles.length);
						console.log('SMILES chop:', smiles);
						}

					//Pick up element and clean out
					var clear_id = '#' + scope.ngid
					clear_element = angular.element(document.querySelector(clear_id));
					clear_element.empty();
		
					if(smiles.length > 0 ) {
						console.log('smiles viewer:', smiles, scope.dim, scope.rotate);
						//Load PDB into Protein Viewer
						function loadPDB(pdb,rotate) {
							console.log('loadPDB');
							var structure = pv.io.pdb(pdb);
							var ligand = structure.select({rnames : ['UNL']});
							//Attach protein viewer to element
							var viewer = pv.Viewer(document.getElementById(scope.ngid), { quality : 'high', width: 'auto', height : 'auto', antialias : true, outline : false, });
							viewer.on('viewerReady', function() {
								viewer.clear();
								viewer.ballsAndSticks('ligand', ligand);
								viewer.centerOn(ligand);
								viewer.autoZoom();
								//viewer.fitParent();
								viewer.spin(rotate);
								});
							};

						//Send for structure
						url = 'https://' + scope.external + '/chemistry/smiles_post/' + scope.dim;
						$http.post(url, angular.toJson({"smiles":smiles}) )
							.success(function(data) {
								console.log('Triggered: smiles2PDB');
								})
							.error(function (data, status) {
								console.log('Error',status,data);
								})
							.then(function(data) {
								pdb = data.data;
								//console.log(pdb);
								loadPDB(pdb,scope.rotate);
								});
						}
					}

				},true);//End watch

    		},

    	}
		return molObject;
	});

//PDB
landingApp.directive('pdbViewer', function ($parse, $http) {
    var molObject = {
		restrict: 'E',
		replace: true,
		scope: {pdb: '=',
				dim: '=',
				rotate: '='},
		link: function (scope, element, attrs) {
			scope.$watch('pdb', function(pdb) {
				clear_element = angular.element(document.querySelector('#threepdb'));
				clear_element.empty();
				
				var viewer = pv.Viewer(document.getElementById('threepdb'), { quality : 'high', width: 'auto', height : 'auto', antialias : true, outline : false});
				var structure;
				function pdb_preset() {
					viewer.clear();
					var ligand = structure.select({rnames : []});
					viewer.ballsAndSticks('ligand', ligand);
					viewer.cartoon('protein', structure, { color: color.ssSuccession() });
					}
				function load_pdb(pdb_id) {
					var xhr = new XMLHttpRequest();
					xhr.open('GET', 'http://files.rcsb.org/view/'+pdb_id+'.pdb');
					xhr.setRequestHeader('Content-type', 'application/x-pdb');
					xhr.onreadystatechange = function() {
						if (xhr.readyState == 4) {
							if (xhr.status == 200) {
								structure = pv.io.pdb(xhr.responseText);
								pdb_preset();
								viewer.centerOn(structure);
								viewer.autoZoom();
								viewer.fitParent();
								viewer.spin(scope.rotate);
								}
							else {
								document.getElementById("threepdb").innerHTML = '<small>Could not load PDB structure with identifier:' + pdb_id + '</small>';
								}
							}
						}
					xhr.send();
					}
				
				window.onresize = function(event) {
					viewer.fitParent();
					}
				window.ondblclick = function(event) {
					if(viewer.spin() === true) {
						viewer.spin(false);
						}
					else {
						viewer.spin(true);
						}
					}
				load_pdb(scope.pdb);

				});
    		},
    	}
		return molObject;
	});
