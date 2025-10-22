#!/usr/bin/env bash

SCRIPT_DIR="/home/pluto/sites/FEMBasicsLectureNotes/systemscripts"

cp -p "$SCRIPT_DIR/plutosliderserver.service" "/etc/systemd/system/plutosliderserver.service"

cp -p "$SCRIPT_DIR/pluto.conf" "/etc/nginx/conf.d/pluto.conf"
