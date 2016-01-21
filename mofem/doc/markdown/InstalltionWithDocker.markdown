##Installation with Docker##

Docker is an open platform that allows for the distribution and deployment of
applications across different systems. In the context of MoFEM it allows for
the distribution of pre-compiled libraries for the Docker platform that can
then run on any system. The Docker platform works by using facilities provided
by the Linux kernel to provide lightweight containers and thereby avoiding the
need to run costly virtual machines. It’s through the use of containers that
MoFEM is compiled and ran.

In Mac OS X a lightweight Linux distribution is virtualized to run the Docker
containers in.

##1. Download and Install Docker##

First install docker as per the instructions here: https://docs.docker.com/installation/#installation

##2. Setup Aliases##

The setup of Docker and the use of MoFEM within it will be through the command
line interface. To ease the installation and continual usage aliases have been
set up for common Docker commands. The aliases have shortened the commands to
more memorable words.

If you wish to use aliases then add the text below to your shell profile, most
likely bash. Note the name of this file will be different depending on the type
and version of Mac OS X you are running and is likely to be a hidden file. Here
are common locations:

Mac OS X: ~/.bash_profile or ~/.profile

Throughout this guide the Docker commands are given in their alias form but
their full form can be found below.

Note the image from which the containers are based upon may be updated. To use
the new version change the dockpull alias command to the correct version.

~~~~~~
#Docker Aliases
#--------------

#Set my hostname
DOCK_HOSTNAME="mofem_host"

#Initialise docker variables in the current shell
alias dockshellinit='`boot2docker shellinit 2>/dev/null | grep export`'

#List all containers
alias docklist='docker ps -a'

#Pull the ubuntu_mofem image from remote Docker server. Note this version number may change
alias dockpull='docker pull likask/ubuntu_mofem:v0.6'

#Update docker time: After the system waking from sleep Docker’s time is wrong and this can affect building/compilation.
alias docktime='/usr/local/bin/boot2docker ssh sudo ntpclient -s -h pool.ntp.org'

#Create and start a temporary container, and remove it on exit - all changes will be lost
alias dockcreatetmp='docker run -it --rm --hostname $DOCK_HOSTNAME -v $HOME:/mnt/home --volumes-from mofem_build mofem_host/ubuntu_mofem:v0.6 /bin/bash'

#Create and start container with the name "sshd" and mount build and data volumes
alias dockcreatesshd='docker run -d -P -p 2222:22 --name sshd --hostname $DOCK_HOSTNAME -v $HOME:/mnt/home --volumes-from mofem_build mofem_host/ubuntu_mofem:v0.6 /usr/sbin/sshd -D'
#start sshd container before connecting
alias dockstartsshd='docker start sshd'
#Login into sshd container using ssh and xserver (the command below finds the IP address)
alias sshdock='ssh -X root@$(echo $(boot2docker ip 2>&1 | sed "s/^.*: //g")) -p 2222'
#Remove container “sshd”
alias dockstopsshd='docker rm -f sshd'
#Transfer your public key to Docker
alias dockauthorize='scp -P 2222 ~/.ssh/id_rsa.pub root@$(echo $(boot2docker ip 2>&1 | sed "s/^.*: //g")):.ssh/authorized_keys'
~~~~~~
-----------------------------

##3. Setting Up Docker##

There are two ways of setting up Docker: launching the boot2docker app or setting up a terminal. The instructions below are for the terminal only. Note if you use the boot2docker app the above aliases are not loaded in the terminal which is automatically opened.


###3.1 Initialise boot2docker###

This sets up the virtualization layer and so is only required the very first time.

~~~~~~
boot2docker init
~~~~~~

###3.2 Start boot2docker###

This starts the virtual machine for the lightweight Linux OS. This needs to be done each time Docker is stopped, such as after a system reboot.

~~~~~~
boot2docker start
~~~~~~

To ascertain the current status of the virtual machine run:
~~~~~~
boot2docker status
~~~~~~

###4.2 Initialise Required Variables###

This sets up the neccesary environment variables for boot2docker to work with Docker.

~~~~~~
dockshellinit
~~~~~~

-----------------------------
##4. Download MoFEM Docker Image##

Download the ubuntu_mofem image which contains the pre-compiled libraries for MoFEM and also use this command to update to the most recent release image.

~~~~~~
dockpull
~~~~~~

##6. Setup MoFEM in Docker##

Setting up MoFEM in Docker requires the creation of 2 containers from 2
different images and mounting the host OS in one of them. The diagram below
illustrates how this works with the relevant commands given thereafter.

![Docker_ImagesContainersVolumes](DockerImagesContainersVolumes.png)

###4.1 Create mofem_build Container###

The user generated data is kept on 2 volumes mounted on this container. This
generated data could be the compiled MoFEM code.

To create the mofem_build container from the Ubuntu image with the 2 volumes
run:

~~~~~~
docker run --name mofem_build  -v /build -v /data ubuntu /bin/bash
~~~~~~

###5.2 Create sshd Container for SSH Access###

To access the 2 volumes in the mofem_build container they are mounted into an
sshd container. This container has the precompiled libraries necessary for
MoFEM  and can be accessed via ssh. This keeps the data volumes separate from
the precompiled libraries and so when the libraries are updated the container
can be replaced with the new one and the user’s data is not lost.

To create the sshd container run:

~~~~~~
dockcreatesshd
~~~~~~

###6.3 SSH Access to SSHD Container###

Start sshd container before connecting:
~~~~~~
dockstartsshd

~~~~~~

Connect to ubuntu_mofem container using ssh:

~~~~~~
sshdock

~~~~~~

Password: mofem

To stop the need to enter a password every time you log into the container you
can generate a public-private key pair and copy the public key to the sshd
container. To generate a public key in Mac OS X or Linux:

~~~~~~
ssh-keygen -t rsa
~~~~~~

This will prompt you for a password which you should leave blank. Then copy
your new public key to your sshd container:

~~~~~~
dockauthorize
~~~~~~

You should now be able to log-in without having to enter the password.

###6.4 Time De-Sync Problem###

When awaking the host system from sleep the Docker clock may de-sync with the
host clock. This can create compilation warnings and problems. To solve this
use the command below in the host OS.

~~~~~~
docktime
~~~~~~

##6. Building MoFEM in Docker##

###6.1 Building MoFEM v0.2###

In v0.2 the core MoFEM libraries are built separately from the user’s modules.
This adds further complications during the building process which will be
explained below. The reason for this is it allows the MoFEM core to be used as
a library for other projects.

Building should be done through the sshd container into the ''$MOFEM_INSTALL_DIR'' mounted
volume and the source code on the host system can be found in /mnt/home

Please note that the location of the MoFEM source code needs to be specified by
the user due to differing locations. The placeholder for this location in the
cmake command is: `$MOFEM_SOURCE_CODE_DIR` and should be changed to the
directory in Mac OS X from the home directory.

First we need to clone MoFEM source directory into your ```$HOME``` on your
native system (not docer). You could be able to edit source code using your
favorite editor on your native system and manage *Git* repository.

Cloning MoFEM sourcecode (note that following line you should execute in your native operating system):
~~~~~~
cd $HOME
git clone https://bitbucket.org/likask/mofem-cephas.git mofem-cephas
~~~~~~

####6.1.1 Release####

Setup file structure:

~~~~~~
export MOFEM_SOURCE_CODE_DIR=/home/home
export MOFEM_INSTALL_DIR=/build
mkdir $MOFEM_INSTALL_DIR

cd $MOFEM_INSTALL_DIR
mkdir release
cd release
mkdir lib
mkdir usr_mods
~~~~~~

Cmake command:

~~~~~~
cd $MOFEM_INSTALL_DIR/release/lib

cmake \
-DCMAKE_BUILD_TYPE=Release \
-DCMAKE_CXX_FLAGS="-Wall" \
-DPETSC_DIR=/opt/petsc -DPETSC_ARCH=arch-linux2-c-opt \
-DMOAB_DIR=/opt/petsc/arch-linux2-c-opt \
-DADOL-C_DIR=/usr \
-DTETGEN_DIR=/opt/tetgen1.5.0 \
-DBUILD_SHARED_LIBS=yes \
-DCMAKE_INSTALL_PREFIX=$MOFEM_INSTALL_DIR/release/usr_mods \
/mnt/home/*mofem_source_code_directory*
~~~~~~

To build the MoFEM libraries run:

~~~~~~
make -j 4
~~~~~~
Install user’s modules:

~~~~~~
make install
~~~~~~

To build the user’s modules run:

~~~~~~
cd $MOFEM_INSTALL_DIR/release/users_modules

# Install "obsolete module" code needed by some examples. See list
# of available users modules on main page.
git clone https://bitbucket.org/likask/mofem_um_obsolete obsolete

cmake -DCMAKE_BUILD_TYPE=Release users_modules
make -j 4
~~~~~~

Note: building all of the users_modules would take considerable time and so
it’s better to only build that which you need to use.

Note 2: when changing cmake files or files which are copied to the build
directory, both MoFEM libraries and user's modules need to be re-built.

####6.1.2 Debug####

Setup file structure:

~~~~~~
cd $MOFEM_INSTALL_DIR
mkdir debug
cd debug
mkdir lib
mkdir usr_mods
~~~~~~

Cmake command:

~~~~~~
cd $MOFEM_INSTALL_DIR/debug/lib

cmake \
-DCMAKE_BUILD_TYPE=Debug\
-DCMAKE_CXX_FLAGS="-Wall" \
-DPETSC_DIR=/opt/petsc -DPETSC_ARCH=arch-linux2-c-debug \
-DMOAB_DIR=/opt/petsc/arch-linux2-c-debug \
-DADOL-C_DIR=/usr \
-DTETGEN_DIR=/opt/tetgen1.5.0 \
-DBUILD_SHARED_LIBS=yes \
-DCMAKE_INSTALL_PREFIX=$MOFEM_INSTALL_DIR/debug/usr_mods \
$MOFEM_SOURCE_CODE_DIR
~~~~~~

To build the MoFEM libraries run:

~~~~~~
make -j 4
~~~~~~
Install user’s modules:

~~~~~~
make install
~~~~~~

To build the user’s modules run:

~~~~~~
cd $MOFEM_INSTALL_DIR/debug/usr_mods

cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-Wall" users_modules
make -j 4
~~~~~~

Note: building all of the users_modules would take considerable time and so
it’s better to only build that which you need to use.

Note 2: when changing cmake files or files which are copied to the build
directory, both MoFEM libraries and user's modules need to be re-built.

##7. Copying Output

The output files generated are part of a Docker container but are required on
the host OS. To get the files there are several options:

1. Copy from your container to your local file system using scp or rsync.

2. Create a symbolic link of the executable in Docker to your host system. Note this will need to include any input files.

3. Change the location of the install directory to the host system
