Installation with Docker (Linux, macOS and some versions of Windows) {#install_docker}
=======================================================================

Docker is an open platform that allows for the distribution and deployment of
applications across different systems. In the context of MoFEM it allows for
the distribution of pre-compiled libraries for the Docker platform that can
then run on any system. The Docker platform works by using facilities provided
by the Linux kernel to provide lightweight containers and thereby avoiding the
need to run costly virtual machines. It’s through the use of containers that
MoFEM is compiled and run.

In macOS a lightweight Linux distribution is virtualized to run the Docker
containers in.
Entire installation procedure is also presented on [Youtube](https://www.youtube.com/watch?v=6opfKER7JHA) video.

[TOC]

# Download and Install Docker {#docker_install}

First install docker as per the instructions here: [https://docs.docker.com/installation/#installation](https://docs.docker.com/installation/#installation)

# Clone mofem repository {#docker_clone}
To build MoFEM the source code need to be downloaded. The best method to do it is
to clone repository
~~~~~~
cd $HOME
git clone --recurse-submodules https://bitbucket.org/likask/mofem-cephas.git
~~~~~~
You can clone specific branch, for example development branch with most up to
date bug fixes, new features, efficiency and functionality improvements
~~~~~~
git clone -b develop --recurse-submodules https://bitbucket.org/likask/mofem-cephas.git
~~~~~~

# Build docker image {#docker_image}

The propose of docker is to build MoFEM in the controlled programming
environment. Such an environment allows to reduce errors quickly, and make
subsequent updates of libraries.

MoFEM is building in two stages; the first build is the core library docker
image, and the second is built docker volume which users modules. This
environment is dedicated to the user who likes to add user module or modify
an existing one. Is less likely for the end-user, which want to run
simulations, in such case we can cate precompiled image which can be
downloaded from docker hub,

Next step of installation is to configure and compile MoFEM. First command creates
*mofem build image*; 
~~~~~~
docker build -t mofem_build --force-rm=true --file=$HOME/mofem-cephas/Dockerfile-build $HOME/mofem-cephas
~~~~~~
Second command creates *mofem build container* which
contains *mofem_build volume*. Volume in container will be shared between other
containers were MoFEM is compiled and run;
~~~~~~
docker run --name mofem_build mofem_build
~~~~~~
This command compiles users modules and runs tests. However, results of
compilation are not part of the container but are stored in the volume.
Several docker containers can share volume, by the option, *--volumes-from
mofem_build*, and use it as space where data can be easily exchanged.

If you do not exactly understand what is *docker image*, *docker container* and
*docker volume* do not worry. You do need to only know how to run and develop
code in docker, how to do it is explained in below. However if you like to fully explore
features available by running MoFEM in docker and utilize its full potential pleas look into
documentation in [Docker User Guide](https://docs.docker.com/engine/userguide/)

# Running docker container {#docker_run_container}

Installation is at that point done, now you can run docker container and
run some code.

To run code a *work container* need to be started, container mount *mofem build
volume* from container which has been created in the previous step
~~~~~~
docker run --rm=true -it --volumes-from mofem_build mofem_build /bin/bash
~~~~~~
However, if you need access or exchange data with home directory, or you like
to recompile MoFEM with source changes on your host hard drive, you can *run*
docker container as follows
~~~~~~
docker run --rm=true -it \
--volumes-from mofem_build \
-v $HOME/mofem-cephas/mofem:/mofem \
-v $HOME:$HOME \
-e HOSTHOME=$HOME mofem_build /bin/bash
~~~~~~
After execution of above command you are working inside docker, this is isolated
system hosted by your OS (MacOSX, Linux or Windows). You can run several
containers like this at once by executing above command in available terminal.

The *work container* mounts *mofem source directory* into *mofem* directory and
your home directory.

MoFEM docker container mount volumes as follows as follows:
- Changes in root direct make only effect for running this container.
- Changes in directory *mofem_build* are shared between other docker containers.
- Changes in home directory in container or host system are shared.

In container new we can compile linear elastic example and calculate simple problem
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
keep old versions locally or upload them int [Docker Hub](https://hub.docker.com/r/likask/ubuntu_mofem/).

# Adding additional user modules
Additional user modules can be added according to the user/developer
requirements. If the user/developer wants to add additional user module, e.g.
 "homogenisation". The steps are as follows:

Clone the homogenisation repository in the user_modules directory of your
source code on the host using:
~~~~~~
cd $HOME/mofem-cephas/mofem/users_modules/
git clone https://bitbucket.org/likask/mofem_um_homogenisation.git homogenisation
~~~~~~

The homogenisation user modules make use of "small_strain_plasticity" and
 "obsolete" user modules and therefore these should also be added to the
 user_modules directory using:

~~~~~~
git clone https://bitbucket.org/likask/mofem_um_small_strain_plasticity.git small_strain_plasticity
git clone https://bitbucket.org/likask/mofem_um_obsolete.git obsolete
~~~~~~

These additional users modules need to be configured and compile before use.
Use the following set of commands in the mofem_build/um/ directory on the
docker:

~~~~~~
cd ../mofem_build/um/

# Configuration:
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXE_LINKER_FLAGS="-L$MOFEM_INSTALL_DIR/local/lib" users_modules

# Build:
make -j4
~~~~~~


# What you will need on host system {#docker_prerequisites}

- Post processor to visualise results. We using [ParaView](http://www.paraview.org)
however you can find good alternatives like [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/).

- If you going to write your modules or modify existing MoFEM modules you will need
text editor. We recommend [Atom](https://atom.io).

- You will need some basic tools to make plots (f.e.
  [gnuplot](http://www.gnuplot.info)) or work with output files, tools like
  grep, [sed](https://en.wikipedia.org/wiki/Sed) or
  [awk](https://en.wikipedia.org/wiki/AWK). If you are working in Linux, simply
  install appropriate packages. If you are a macOS user,
  we recommend to install [HomeBrew](http://brew.sh), which installs missing
  packages into macOS system.

# Contact {#docker_contact}

Any problems with this installation, please contact us by [mofem-group@googlegroups.com](https://groups.google.com/forum/#!forum/mofem-group)
or on Slack [MoFEM Slack](https://mofem.slack.com/).
