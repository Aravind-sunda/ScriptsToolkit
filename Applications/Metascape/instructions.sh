# download the zip file from metascape 
# copy the zip file to the directory that you want to install the docker env
# follow the instructions in the following link:
# https://metascape.org/gp/index.html#/menu/msbio

#gitignore will have ignore the directory
unzip msbio_v3.5.20260201.zip  -d metascape/
cd metascape/

# using docker: if running locally on personal computer
bin/install.sh

# using singularity: if running on a cluster
bin/install_sms.sh


# _BACKGROUND add this to the multigene list file, so that metascape can recognize the background gene s