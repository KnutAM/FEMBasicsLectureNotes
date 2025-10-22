# This script stops and restarts the service, and must be run with sudo rights
systemctl daemon-reload
echo "Stopping the pluto server and waiting 2 s"
systemctl stop pluto-server
sleep 2
echo "Restarting the pluto server"
systemctl start pluto-server
systemctl enable pluto-server
echo "Restart complete"
