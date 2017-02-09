from flask import Flask, make_response, url_for, request
from flask_restful import Resource
from flask_socketio import SocketIO, emit

#Initialise app
app = Flask(__name__)
app.config['SECRET_KEY'] = 'big_secret'

print 'Socket handler for messaging and push'

#-----------------------------------------------------------------------------------------
#SOCKET APP
socketio = SocketIO(app)
@socketio.on('connect', namespace='/mq')
def test_connect():
    print 'Client connected'

@app.route('/socket', methods=['POST'])
def socketHandler():
	json_data = request.get_json(force=True)
	handle = json_data['handle'] #<port>_<id>
	message = json_data['message']
	port = json_data['port'] #Use this to define a port-specific namespace
	print 'EMIT:', handle, message, port
	try:
		socketio.emit(handle, message,  namespace='/mq_' + port)
		return str(True)
	except:
		return str(False)
#-----------------------------------------------------------------------------------------


# Run the app.
if __name__ == '__main__':
	socketio.run(app,host='0.0.0.0',debug=False, port=8082)