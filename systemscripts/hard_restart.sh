# This script stops and restarts the service, and must be run with sudo rights
systemctl daemon-reload
echo "Stopping the pluto server and waiting 2 s"
systemctl stop plutosliderserver.service
sleep 2
echo "Restarting the pluto server"
systemctl start plutosliderserver.service
systemctl enable --now plutosliderserver.service
nginx -t
systemctl enable --now nginx
echo "Restart complete"
