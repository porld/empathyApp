#!/usr/bin/env bash

#Fetch Python modules
sudo pip install --upgrade pip
sudo pip install flask, flask-socketio, flask-httpauth, flask-restful, flask-mail, passlib, passwordmeter, gevent

#Run app
python landing.py