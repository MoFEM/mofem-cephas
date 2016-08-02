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

##Download and Install Docker##

First install docker as per the instructions here: https://docs.docker.com/installation/#installation

##Setup Aliases##

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

# Set my hostname
DOCK_HOSTNAME="mofemdocker"

# Pull mofem container
alias dockpull='docker pull likask/ubuntu_mofem:v0.7'

# Start mofem docker
alias dockmofem='docker run -it --hostname $DOCK_HOSTNAME -v $HOME:/mnt/home -v mofem_build:/build likask/ubuntu_mofem:v0.7 /bin/bash'

# Remove last container
alias dockrmlast='docker rm `docker ps -l -q`'

# Start container with name work
alias dockstartwork='docker run -it --name work --hostname $DOCK_HOSTNAME -v $HOME:/mnt/home -v mofem_build:/build likask/ubuntu_mofem:v0.7 /bin/bash'
# ReAttach work container
alias dockattachwork='docker start -ai work'
# Remove container work
alias dockrmwork='docker rm work'

# Start container and remove it on exit
alias docktmpwork='docker run -it --rm --hostname $DOCK_HOSTNAME -v $HOME:/mnt/home -v mofem_build:/build likask/ubuntu_mofem:v0.7 /bin/bash'

# List all containers
alias dockpsall='docker ps -a'

#start sshd

docker-ip() {
  docker inspect --format '{{ .NetworkSettings.Ports.HostIp.22/tcp }}' "$@"
}


alias dockstartsshd='docker run -d -P -p 2222:22 --name sshd --hostname $DOCK_HOSTNAME -v $HOME:/mnt/home -v mofem_build:/build likask/ubuntu_mofem:v0.7 /usr/sbin/sshd -D'
alias dockstopsshd='docker rm -f sshd'
alias dockautorize='scp -P 2222 ~/.ssh/id_rsa.pub root@$0.0.0.0:.ssh/authorized_keys'
alias sshdock='ssh -X root@0.0.0.0 -p 2222'

~~~~~~
-----------------------------

##Download MoFEM Docker Image##

Download the ubuntu_mofem image which contains the pre-compiled libraries for MoFEM and also use this command to update to the most recent release image.

~~~~~~
docker pull likask/ubuntu_mofem:v0.7
~~~~~~

##Setup MoFEM in Docker##

Setting up MoFEM in Docker requires the creation of 2 containers from 2
different images and mounting the host OS in one of them. The diagram below
illustrates how this works with the relevant commands given thereafter.

![Docker_ImagesContainersVolumes](DockerImagesContainersVolumes.png)

###Create mofem_build Container###

The user generated data is kept on 2 volumes mounted on this container. This
generated data could be the compiled MoFEM code.

To create the mofem_build container from the Ubuntu image with the 2 volumes
run:

~~~~~~
docker volume create --name mofem_build
~~~~~~

Now you cans start to work in mofem container
~~~~~
docker run -it --hostname $DOCK_HOSTNAME -v $HOME:/mnt/home -v mofem_build:/build likask/ubuntu_mofem:v0.7
~~~~~

##Building MoFEM in Docker##

###Building MoFEM ###

From v0.2 the core MoFEM libraries are built separately from the user’s modules.
This adds further complications during the building process which will be
explained below. The reason for this is it allows the MoFEM core to be used as
a library for other projects.

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
git clone https://bitbucket.org/likask/mofem-cephas.git
~~~~~~

####Release####

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
-DCMAKE_Fortran_COMPILER=/usr/bin/gfortran \
-DCMAKE_CXX_FLAGS="-Wall" \
-DPETSC_DIR=/opt/petsc -DPETSC_ARCH=arch-linux2-c-opt \
-DMOAB_DIR=/opt/petsc/arch-linux2-c-opt \
-DADOL-C_DIR=/usr \
-DTETGEN_DIR=/opt/tetgen1.5.0 \
-DBUILD_SHARED_LIBS=yes \
-DCMAKE_INSTALL_PREFIX=$MOFEM_INSTALL_DIR/release/users_modules \
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

####Debug####

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
-DCMAKE_BUILD_TYPE=Debug \
-DCMAKE_Fortran_COMPILER=/usr/bin/gfortran \
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
