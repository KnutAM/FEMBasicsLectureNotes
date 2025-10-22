#!/usr/bin/env bash

SCRIPT_DIR="/home/pluto/sites/FEMBasicsLectureNotes/systemscripts"

# Allow web server to make outbound connections 
setsebool -P httpd_can_network_connect 1

cp -p "$SCRIPT_DIR/plutosliderserver.service" "/etc/systemd/system/plutosliderserver.service"

cp -p "$SCRIPT_DIR/pluto.conf" "/etc/nginx/conf.d/pluto.conf"

cp -p "$SCRIPT_DIR/certbot-renew.service" /etc/systemd/system/certbot-renew.service

cp -p "$SCRIPT_DIR/certbot-renew.timer" /etc/systemd/system/certbot-renew.timer


systemctl daemon-reload
systemctl enable --now plutosliderserver.service
systemctl enable --now certbot-renew.timer

nginx -t
systemctl enable --now nginx