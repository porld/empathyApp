source activate cheminformatics

#http://werkzeug.pocoo.org/docs/0.10/serving/
#http://stackoverflow.com/questions/7340784/easy-install-pyopenssl-error

#Set up port redirect
sudo /sbin/iptables -t nat -I PREROUTING -p tcp --dport 80 -j REDIRECT --to-port 8080

python empathy_socket.py &
python api_external.py &
python landing.py $1 &