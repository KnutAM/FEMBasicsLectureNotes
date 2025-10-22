#!/usr/bin/env bash
# Quick status check for PlutoSliderServer + NGINX + HTTPS setup on RHEL

set -e

echo "=== PlutoSliderServer status ==="
if systemctl is-active --quiet plutosliderserver; then
    echo "✅ PlutoSliderServer is active"
else
    echo "❌ PlutoSliderServer is NOT running"
    systemctl --no-pager -l status plutosliderserver || true
fi

echo
echo "=== NGINX status ==="
if systemctl is-active --quiet nginx; then
    echo "✅ NGINX is active"
else
    echo "❌ NGINX is NOT running"
    systemctl --no-pager -l status nginx || true
fi

echo
echo "=== Port bindings ==="
sudo ss -ltnp | egrep '(:80 |:443 |:8080)' || echo "⚠️  Expected ports not found"

echo
echo "=== Firewall zones & allowed services ==="
sudo firewall-cmd --get-active-zones
echo
for z in $(sudo firewall-cmd --get-active-zones | awk 'NR%2==1'); do
  echo "--- Zone: $z ---"
  sudo firewall-cmd --zone="$z" --list-services
done

echo
echo "=== Certificate expiry ==="
if sudo test -f /etc/letsencrypt/live/teaching.ims.chalmers.se/fullchain.pem; then
    sudo openssl x509 -enddate -noout -in /etc/letsencrypt/live/teaching.ims.chalmers.se/fullchain.pem
else
    echo "⚠️  Certificate not found at /etc/letsencrypt/live/teaching.ims.chalmers.se/"
fi

echo
echo "=== HTTPS test (local) ==="
curl -ksI https://localhost/ | head -n 1 || echo "⚠️  HTTPS not responding"

echo
echo "=== Summary ==="
echo "If all checks show ✅ and curl returns HTTP/1.1 200 or 301, everything is OK."