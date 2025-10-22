# FEMBasicsLectureNotes
This repository contains notebook lecture notes for the course *FEM Basics* (VSM167) at Chalmers University of Technology. 

To reduce the memory requirement, using too many notebooks should be avoided, it is more efficient to bunch multiple parts into a single notebook. 

*Do not make any manual changes to the git repo on the virtual webserver directly*: Instead, make changes to the [remote](https://github.com/KnutAM/FEMBasicsLectureNotes) either directly or via pull requests / pushing. 

## Details on the setup
**NOTE:** This information is not yet updated to the Chalmers system

* The directory `pluto-slider-server-environment` contains the configurations for the pluto slider server.
* The service file is here: `/etc/systemd/system/pluto-server.service`
* The startup script is here: `/usr/local/bin/pluto-slider-server.sh`


### Restarting the service
Normally, notebooks are updated automatically when they change due to a change in the git repository. However, in case something goes wrong and the service needs to be restarted, one use `devutils/restart_hard.sh` (requires sudo rights). 

### Starting up
To start up the service (first time on a new computer), run 
```
sudo chmod 744 /usr/local/bin/pluto-slider-server.sh
sudo chmod 664 /etc/systemd/system/pluto-server.service
sudo systemctl daemon-reload
sudo systemctl start pluto-server
sudo systemctl enable pluto-server
```
