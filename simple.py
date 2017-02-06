from flask import Flask, make_response, url_for, request

#Initialise app
app = Flask(__name__, static_url_path='')
app.config['SECRET_KEY'] = 'big_secret'

#Landing page
@app.route('/')
def index():
	return 'Hello world'

# Run the app.
if __name__ == '__main__':
	#app.run(debug=True, port=5000)
	app.run(host='0.0.0.0',debug=False, port=80)
