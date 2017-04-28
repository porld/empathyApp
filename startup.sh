#Jamboree app startup script

echo 'Reroute IP tables'
sudo /sbin/iptables -t nat -I PREROUTING -p tcp --dport 443 -j REDIRECT --to-port 8080

echo 'Launch ancillary services'
python empathy_ancillary.py &

echo 'Launch Socket'
python empathy_socket.py &

echo 'Launch core'
python empathy_landing.py $1 &