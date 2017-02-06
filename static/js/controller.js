var empathyApp = angular.module('empathyApp', ['ngTouch','ngSanitize']);

empathyApp.config(['$httpProvider', function($httpProvider) {
	$httpProvider.defaults.useXDomain = true;	
   delete $httpProvider.defaults.headers.common['X-Requested-With'];	
}]);

empathyApp.controller('PanelCtrl', ['$scope', '$http', '$location', '$interval', function ($scope, $http, $location, $interval) {

	//Base url for the static files
	$scope.static_url = "http://localhost:5000";

	//Pick up port from URL
	port = window.location.search.substr('port');
	$scope.port = port.replace('?port=','');
	console.log('Waiting for Neo4j running in Docker on port', $scope.port);
	url = $scope.static_url + '/checkLive/' + $scope.port;
	console.log(url);
	$scope.docker_live = [false,0];
	$scope.docker_error = '';
	$scope.change_pwd = false;
	response = $http.get(url)
		.success(function(data) {
			console.log('Checking Docker status', data);
			$scope.docker_live = data;
			})
		.error(function (data, status) {
			console.log('Docker error',status,data);
			$scope.docker_live = false;
			})
		.then(function (data) {
			console.log('Docker status on port', $scope.port, ':', $scope.docker_live);
			if($scope.docker_live[0]) {
				console.log('Attempt sign in');
				$scope.change_pwd = true;
				}
			else {
				console.log('Could not find Docker on', $scope.port, 'after', $scope.docker_live[1], 'attempts');
				$scope.docker_error = 'Could not locate your database on port ' + $scope.port;
				}
			});
		
	//Log in
	$scope.login = false;
	$scope.username = '';
	$scope.password = '';

	$scope.verify = function(username,password) {
		$scope.passport = {"username":username,"password":password,"port":$scope.port}; //Change this port when Docker is up

		//Sign in
		console.log('Trying to sign in...');
		bundle = {"passport":$scope.passport};
		url = $scope.static_url + '/signIn';
		response = $http.post(url, angular.toJson(bundle) )
			.success(function(data) {
				console.log('Sign in response', data);
				$scope.login = data;
				})
			.error(function (data, status) {
				console.log('Error',status,data);
				return false;
				});
		console.log('Sign in over', response);

		if($scope.login) {
			console.log('Signed in');
			run_app();
			return true;
			}
		else {
			console.log('Could not sign in');
			return false;
			}
		}

	run_app = function() {
		$scope.login = true;
		//Initialise some variables
		console.log('Logged in');
		$scope.label = 'molecule';
		$scope.message_handle = '';
		$scope.list = {};
		$scope.selection = '';
		$scope.record = {};
		$scope.alive = true;
	
		//Connect to ZeroMQ broadcast server
		var socket = io.connect('http://' + document.domain + ':' + location.port + '/mq')
		//Connect to message socket	
		socket.on('connect', function() {
			console.log('Connected to message socket');
			});

		//Pick up broadcasts about list
		socket.on($scope.port + '_' + $scope.label, function(latest_list) {
			console.log('List update');
			$scope.list = latest_list;
			console.log($scope.list);
			$scope.$apply();
			});

		//Wait for things to be in place before fetching
		initialiseList($scope.label,$scope.passport,$scope.port + '_' + $scope.label,$scope.static_url,$http)

		} //End else

	//FUNCTIONS

	//Click function for selection a record
	$scope.makeSelection = function(id) {
		$scope.selection = id;
		$scope.fetchSelection(id,$scope.passport,static_url,$http);
		//Update socket
		$scope.record_socket();
		};

	//Update record property (REST fire, returns via socket.io)
	$scope.updateNode = function(key,value) {
		console.log('Update', $scope.selection, 'with', key, ':', value);
		bundle = {"passport":$scope.passport, "message_handle": $scope.port+'_'+$scope.selection, "id":$scope.selection, "key":key, "value":value};
		url = static_url + '/updateNode';
		$http.post(url, angular.toJson(bundle) )
			.success(function(data) {
				console.log('Triggered', data);
				})
			.error(function (data, status) {
				console.log('Error',status,data);
				});
		};
		
	$scope.changeFunction = function() {
		$scope.change = true;
		};

	//Pick up broadcasts about current record properties
	$scope.record_socket = function() {
		console.log('Refresh socket:', $scope.port, '_', $scope.selection);
		socket.on($scope.port + '_' + $scope.selection, function(record_update) {
			key = record_update['key'];
			value = record_update['value']
			$scope.record[key] = value;
			initialiseList($scope.label,$scope.passport,$scope.port + '_' + $scope.label,static_url,$http);
			$scope.$apply();
			});
		};

	//Node create and destroy functions
	$scope.createNode = function(label) {
		console.log('Creating new', label, 'node...');
		createNode(label,$scope.passport,$scope.port + '_' + $scope.label,$scope.static_url,$http);
		};
	$scope.destroyNode = function() {
		console.log('Destroy', $scope.selection);
		destroyNode($scope.selection,$scope.label,$scope.passport,$scope.port + '_' + $scope.label,$scope.static_url,$http);
		$scope.selection = '';
		$scope.record = '';
		$scope.message_handle = '';
		};

	//Fetch selection
	$scope.fetchSelection = function(selection,passport,static_url,$http) {
		console.log('Fetching node', selection);
		bundle = {"passport":passport, "selection":selection};
		url = static_url + '/fetchSelection';
		$http.post(url, angular.toJson(bundle) )
			.success(function(data) {
				//console.log('Triggered', data);
				$scope.record = data;
				console.log($scope.record);
				})
			.error(function (data, status) {
				console.log('Error',status,data);
				$scope.record = 'error';
				});
		};

	}]);



//Initialise list
initialiseList = function(label,passport,message_handle,static_url,$http) {
	console.log('Fetch initial list');
	bundle = {"passport":passport, "label":label, "message_handle": message_handle};
	url = static_url + '/listNode';
	$http.post(url, angular.toJson(bundle) )
		.success(function(data) {
			console.log('Triggered', data);
			})
		.error(function (data, status) {
			console.log('Error',status,data);
			});
	};
	
//Create new list entry
createNode = function(label,passport,message_handle,static_url,$http) {
	console.log('Fetch initial list');
	bundle = {"passport":passport, "label":label, "message_handle": message_handle};
	url = static_url + '/createNode';
	$http.post(url, angular.toJson(bundle) )
		.success(function(data) {
			console.log('Triggered', data);
			})
		.error(function (data, status) {
			console.log('Error',status,data);
			});
	};

//Destroy list entry
destroyNode = function(target,label,passport,message_handle,static_url,$http) {
	console.log('Destroy', target);
	bundle = {"passport":passport, "label":label, "message_handle": message_handle, "target":target};
	url = static_url + '/destroyNode';
	$http.post(url, angular.toJson(bundle) )
		.success(function(data) {
			console.log('Triggered', data);
			})
		.error(function (data, status) {
			console.log('Error',status,data);
			});
	};
