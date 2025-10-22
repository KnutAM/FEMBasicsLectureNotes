#!/usr/bin/env bash

SCRIPT_DIR="/home/pluto/sites/FEMBasicsLectureNotes/systemscripts"

# Allow web server to make outbound connections 
setsebool -P httpd_can_network_connect 1

cp -p "$SCRIPT_DIR/plutosliderserver.service" "/etc/systemd/system/plutosliderserver.service"

cp -p "$SCRIPT_DIR/pluto.conf" "/etc/nginx/conf.d/pluto.conf"
