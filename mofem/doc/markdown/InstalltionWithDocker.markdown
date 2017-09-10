##Installation with Docker## {#install_docker}

Docker is an open platform that allows for the distribution and deployment of
applications across different systems. In the context of MoFEM it allows for
the distribution of pre-compiled libraries for the Docker platform that can
then run on any system. The Docker platform works by using facilities provided
by the Linux kernel to provide lightweight containers and thereby avoiding the
need to run costly virtual machines. It’s through the use of containers that
MoFEM is compiled and run.

In Mac OS X a lightweight Linux distribution is virtualized to run the Docker
containers in.
Entire installation procedure is also presented on [Youtube](https://www.youtube.com/watch?v=6opfKER7JHA) video.

##Download and Install Docker##

First install docker as per the instructions here: [https://docs.docker.com/installation/#installation](https://docs.docker.com/installation/#installation)

##Clone mofem repository
To build MoFEM the source code need to be downloaded. The best method to do it is
to clone repository
~~~~~~
cd $HOME
git clone https://bitbucket.org/likask/mofem-cephas.git
~~~~~~

##Build docker image

Next step of installation is to configure and compile MoFEM. First command creates
*mofem build image*. Second command creates *mofem build container* which
contains *mofem_build volume*. Volume in container will be shared between other
containers were MoFEM is compiled and run;
~~~~~~
docker build -t mofem_build --force-rm=true --file=$HOME/mofem-cephas/Dockerfile-build $HOME/mofem-cephas
docker run --name mofem_build mofem_build
~~~~~~

If you do not exactly understand what is *docker image*, *docker container* and
*docker volume* do not worry. You do need to only know how to run and develop
code in docker, how to do it is explained in below. However if you like to fully explore
features available by running MoFEM in docker and utilize its full potential pleas look into
documentation in [Docker User Guide](https://docs.docker.com/engine/userguide/)

##Running docker container

Installation is at that point done, now you can run docker container and
run some code.

To run code a *work container* need to be started, container mount *mofem build
volume* from container which has been created in the previous step

~~~~~~
docker run --rm=true -it --volumes-from mofem_build  -v $HOME/mofem-cephas/mofem:/mofem -v $HOME:$HOME -e HOSTHOME=$HOME mofem_build /bin/bash
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

##What you will need on host system

- Post processor to visualise results. We using [ParaView](http://www.paraview.org)
however you can find good alternatives like [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/).

- If you going to write your modules or modify existing MoFEM modules you will need
text editor. We recommend [Atom](https://atom.io).

- You will need some basic tools make plots (f.e. [gnuplot](http://www.gnuplot.info)) or work with output files, tools like grep, [sed](https://en.wikipedia.org/wiki/Sed) or [awk](https://en.wikipedia.org/wiki/AWK). If you working in Linux simply install appropriate packages. If you are MacOS X user
we recommend to install [HomeBrew](http://brew.sh), which install missing packages into
MacOS X system.

##Contact

Any problems with this installation, please contact us by [mofem-group@googlegroups.com](https://groups.google.com/forum/#!forum/mofem-group)
or on Slack [MoFEM Slack](https://mofem.slack.com/).
