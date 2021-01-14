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

[TOC]

# Watch Docker & Spack installation tutorial

Flow link [https://youtu.be/I201ZaAnQG8](https://youtu.be/I201ZaAnQG8). In 
tutorial we briefly explain what is following bellow. 
# Download and Install Docker {#docker_install}

First, install docker as per the instructions here: 
[https://docs.docker.com/installation/#installation](https://docs.docker.com/installation/#installation)


# Running docker container {#docker_run_container}

Installation is at that point done; now you can run docker container and run
some code. You have to start by pulling MoFEM image from Docker Hub
~~~~~~
docker pull likask/mofem-spack-build
~~~~~~

To run code a *work container* need to be started, container mount *mofem build
volume* from container which has been created in the previous step
~~~~~~
docker run \
  -v $HOME:/host_home \
  --rm=true -it likask/mofem-spack-build 
~~~~~~

MoFEM docker container mount volumes as follows as follows:
- Changes in root direct make only effect for running this container.
- Changes in directory *mofem_install* are shared between other docker containers.
- Changes in home directory in container or host system are shared.

In container new, we can run a linear elastic example and calculate the
simple problem
~~~~~~
cd /mofem_install/um_view/elasticity/
mpirun -np 2 --allow-run-as-root \
  ./elasticity \
  -my_file LShape.h5m -my_order 2
mbconvert out_skin.h5m /host_home/out_skin.vtk
~~~~~~
In the above example, two processors are used to do calculations. A number of
processors which can be used effectively depends on your hardware and set up of
the docker container. For propose of this documentation, the program is executed in
*mofem_build* volume, thus resulting file *out.h5m* can be accessed from other
containers, however, is not visible from the host file system. The last command which creates VTK output file save results to the HOME directory of your host system.

# Making changes permanent

Changes which you do in Docker container are temporary. If you restart the computer,
or container content of the container is reset to initial state. If you like to 
keep changes in container or share volume of the container between another container
you can do as follows:
~~~~~~~
docker run --name mofem_volume likask/mofem-spack-build
~~~~~~~

Now we can run Docker container and attach the previously created volume to it
~~~~~~~
docker run \
  --name mofem_develop \
  -v $HOME:/host_home \
  --volumes-from mofem_volume \
  --rm=true -it likask/mofem-spack-build /bin/bash
~~~~~~~

You can run several containers and link them to the same volume. All changes in the volume are permanent unless you delete the volume. You can also create several
volumes and attach different names to it. 

# Developing with VSCode {#docker_vscode}

Run docker container in terminal, as above, and follow the tutorial in
[Attach to a running
container](https://code.visualstudio.com/docs/remote/attach-container) Source
code in VSCode, in Docker running docker container, is located in the root
directory in `/mofem_install/mofem-cephas`. Open this folder in VSCode to
edit source.

We already build created Release version of users modules, which you can modify,
compile and develop. However, for development purposes is useful to have
compiled version which debugging information. You can configure the version of users
modules with debugging information as follows,
~~~~~~
spack dev-build \
  --b build \
  --test root \
  --source-path $MOFEM_UM_SRC_DIR \
  mofem-users-modules@develop+docker \
  build_type=Debug   \
  ^/$MOFEM_CEPHAS_HASH
~~~~~~
Have above done, we can find directory `/mofem_install/um-build-Debug-fifofez`. 
Now, for example, we can look into the tutorial 
@ref basic_tutorials_poisson_homogeneous "SLC-1" oisson's equation, 
and follow instruction, cd directory, make code, partition mesh, and run the example.
~~~~~~
cd /mofem_install/um-build-Debug-fifofez/tutorials/scl-1
make
mofem_part \
  -my_file mesh2d.cub \
  -output_file mesh2d.h5m \
  -my_nparts 2 -dim 2 -adj_dim 1
mpirun -np 2 --allow-run-as-root \
  ./poisson_2d_homogeneous -file_name mesh2d.h5m -order 2
~~~~~~
Next, create VTK file in your home directory 
~~~~~~~
mbconvert out_result.h5m /host_home/out_result.vtk
~~~~~~~ 
and finally, open it in ParaView [ParaView](https://www.paraview.org/download/).

# Running Jupyter notebook {#docker_jupyter}


Pull MoFEM images and create mofem_volume
~~~~~~
docker run --name mofem_volume likask/mofem-spack-build
~~~~~~
and run container as follows
~~~~~~
docker run \
  -p 8888:8888  \
  --name mofem_jupyetr \
  -v $HOME:/host_home  \
  --volumes-from mofem_volume \
  --rm=true -ti \
  likask/mofem-spack-jupyter
~~~~~~
As results in terminal window you will see something similar to following:
~~~~~~
[I 20:33:35.653 NotebookApp] Jupyter Notebook 6.1.6 is running at:
[I 20:33:35.653 NotebookApp] http://localhost:8888/?token=bed64e1d532ab00f402d56432193a80f6aace612d4ffbec2
[I 20:33:35.653 NotebookApp]  or http://127.0.0.1:8888/?token=bed64e1d532ab00f402d56432193a80f6aace612d4ffbec2
[I 20:33:35.653 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
~~~~~~
Copy address into web-browser:
~~~~~~
http://127.0.0.1:8888/?token=bed64e1d532ab00f402d56432193a80f6aace612d4ffbec2
~~~~~~

You can start docker without starting volume, as a isolated system, as follows,
~~~~~~
docker run \
  -p 8888:8888  \
  --name mofem_jupyter \
  --rm=true -ti \
  likask/mofem-spack-jupyter
~~~~~~
# Clone MoFEM repository {#docker_clone}

To build MoFEM the source code need to be downloaded. The best method to do it is
to clone repository
~~~~~~
mkdir -p $HOME/mofem_install
cd $HOME/mofem_install
git clone --recurse-submodules https://bitbucket.org/likask/mofem-cephas.git
~~~~~~
You can clone specific branch, for example, development branch with most up to
date bug fixes, new features, efficiency and functionality improvements
~~~~~~
cd $HOME
mkdir -p mofem_install
git clone -b lukasz/develop --recurse-submodules https://bitbucket.org/likask/mofem-cephas.git
~~~~~~

Now you can use docker container with the source from your host system as follows
~~~~~~
docker run \
  --name mofem_develop \
  -v $HOME:/host_home \
  -v $HOME/mofem_install/mofem-cephas:/mofem_install/mofem-cephas \
  --volumes-from mofem_volume \
  --rm=true -it likask/mofem-spack-build /bin/bash
~~~~~~

# Build docker image {#docker_image}

Sometimes you like to build MoFEM environment, or build MoFEM from scratch, instead
downloading preinstalled code from Docker Hub. If you familiar with Docker, you
can modify `Dockerfile-spack-env` to build MoFEM dependent libraries. Or edit
`Dockerfile-spack-build` to build and test installation.

You can build MoFEM docker environment image as follows
~~~~~~
cd $HOME/mofem_install
docker build -t likask/mofem-spack-env -f Dockerfile-spack-env ,
~~~~~~
Once docker image is built, you build image can be created
~~~~~~
cd $HOME/mofem_install
docker build -t likask/mofem-spack-build -f Dockerfile-spack-build .
~~~~~~
and Jupyter notebook image
~~~~~~
cd $HOME/mofem_install
docker build -t likask/mofem-spack-jupyter -f Dockerfile-spack-jupyter .
~~~~~~

If you do not exactly understand what is *docker image*, *docker container* and
*docker volume* do not worry. You do need to only know how to run and develop
code in docker, how to do it is explained in below. If you like to fully explore
features available by running MoFEM in docker and utilize its full potential pleas look into documentation in [Docker User Guide](https://docs.docker.com/engine/userguide/)


# What you will need on host system {#docker_prerequisites}

- Postprocessor to visualise results. We using [ParaView](http://www.paraview.org)
however, you can find good alternatives like [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/).

- If you are going to write your modules or modify existing MoFEM modules, you will need
a text editor. We recommend [VSCode](https://code.visualstudio.com).
# Contact {#docker_contact}

Any problems with this installation, please contact us by [mofem-group@googlegroups.com](https://groups.google.com/forum/#!forum/mofem-group)
or on Slack [MoFEM Slack](https://mofem.slack.com/).
