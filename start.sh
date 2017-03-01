source activate cheminformatics

#http://werkzeug.pocoo.org/docs/0.10/serving/
#http://stackoverflow.com/questions/7340784/easy-install-pyopenssl-error

#Set up port redirect
echo 'Redirect ports'
sudo /sbin/iptables -t nat -I PREROUTING -p tcp --dport 80 -j REDIRECT --to-port 8080

echo 'Run socket'
python empathy_socket.py &
#echo 'Run external functions'
#python api_external.py &
echo 'Run core'
python landing.py $1 &