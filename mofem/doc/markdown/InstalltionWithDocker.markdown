##Installation with Docker##

Docker is an open platform that allows for the distribution and deployment of
applications across different systems. In the context of MoFEM it allows for
the distribution of pre-compiled libraries for the Docker platform that can
then run on any system. The Docker platform works by using facilities provided
by the Linux kernel to provide lightweight containers and thereby avoiding the
need to run costly virtual machines. It’s through the use of containers that
MoFEM is compiled and run.

In Mac OS X a lightweight Linux distribution is virtualized to run the Docker
containers in.

##Download and Install Docker##

First install docker as per the instructions here: https://docs.docker.com/installation/#installation

##Clone mofem repository

To build MoFEM the source code need to be downloaded. The best way to do it is
to clone repository
~~~~~~
cd $HOME
git clone https://bitbucket.org/likask/mofem-cephas.git
~~~~~~

##Build docker image

Next step of installation is to configure and compile MoFEM. First command creates
*mofem build image*. Second command creates *mofem build container* which
contains *mofem_build volume*. Volume in container will be shared between other
containers and docker *work* runs.
~~~~~~
cd mofem-cephas
docker build -t mofem_build:v0.1 --force-rm=true --file=Dockerfile-build $HOME/mofem-cephas
docker run --name mofem_build mofem_build:v0.1 /bin/bash
~~~~~~

If you do not exactly understand what is *docker image*, *docker container* and
*docker volume* do not worry. You do need to only know how to run and develop
code in docker what is explained in below. However if you like to fully explore
features avilable in docer and utilise its full potential pleas look into
documentation in https://docs.docker.com/engine/userguide/

##Running docker container

Installation is at that point finished. Now you can run docker container and
run some code.

You can start the *work container*, container mount *mofem build volume* from container which
has been in the previous step

    docker run \
    --rm=true -it \
    --volumes-from mofem_build  \
    -v $HOME/mofem-cephas/mofem:/mofem \
    -v $HOME:$HOME \
    -e HOSTHOME=$HOME \
    likask/ubuntu_mofem:latest /bin/bash

After execution of above command you are working inside docker, this is isolated system
hosted by your OS (MacOSX, Linux or Windows).

The *work container* mounts *mofem source directory* into *mofem* directory and
your home directory.

Note that:
- Changes in root direct make only effect for running this container.
- Changes in directory *mofem_build* are shared between other docker containers.
- Changes in home directory in container or host system are shared.

For example we can run make and run linear elasticity example,
~~~~~~
cd /mofem_build/um/basic_finite_elements/elasticity
make
mpirun -np 2 ./elasticity -my_file LShape.h5m -ksp_type gmres -pc_type lu -pc_factor_mat_solver_package mumps -ksp_monitor -my_order 2
mbconvert out.h5m $HOSTHOME/out.vtk
~~~~~~
In above example two processors are used to do calculations. Number of
processors which can be used effectively depends on your hardware and set up of
the docker container. For propose of this documentation program is executed in
*mofem_build* volume, thus resulting file *out.h5m* can be accessed from other
containers however is not visible from host file system. The last command which creates
VTK output file save results to HOME directory of your host system.

Note that working with docker you can work with several versions of MoFEM at once,
keep old versions locally or upload them int docker hub (https://hub.docker.com)

Any problems with this installation, please contact us by cmatgu@googlegroups.com
or on Slack https://mofem.slack.com/.
