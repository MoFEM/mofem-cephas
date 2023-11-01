Installation with Docker (Linux, macOS and some versions of Windows) {#install_docker_user}
=======================================================================

Docker is an open platform that allows for the distribution and deployment of
applications across different systems. In the context of MoFEM it allows for
the distribution of pre-compiled libraries for the Docker platform that can
then run on any system. The Docker platform works by using facilities provided
by the Linux kernel to provide lightweight containers and thereby avoiding the
need to run costly virtual machines. It’s through the use of containers that
MoFEM is compiled and run.

In macOS, a lightweight Linux distribution is virtualized to run the Docker
containers in.

[TOC]

## Download and Install Docker {#docker_install}

You can download and install Docker following instructions on the [Docker installation webpage](https://docs.docker.com/installation/#installation).

# How to get and run the container {#docker_getting_container}

> These instructions are for installing specific versions of MoFEM containers. To have a more advanced installation for developers see [Installation with Docker](#install_docker).

The containers available at the moment are: 

Purpose | Image | Tag 
---|---|---
Workshop2023 | likask/mofem-spack-jupyterhub | Workshop2023
Development | likask/mofem-spack-build | latest

## Using Docker user interface

This section describes how to set up and run the Docker container within the Docker user interface. For doing this in terminal, go to the [next section](#docker_terminal_installation).

After opening Docker you should be able to see and follow the following images, dependent on the operating system and version of Docker you have.

<!-- <img src="./../figures/docker_intro.png" alt="Docker - introductory screen" width="100%"/> -->
<img src="docker_intro.png" alt="Docker - introductory screen" width="100%"/>
<a id='figure_1'></a> 
    <center><b>Figure 1. Docker startup screen.</b></center>

At the top you can see a search bar and if you are connected to internet you can search for the container you want to run. Start typing `likask/mofem-spack-jupyterhub` and you should be able to see Workshop2023 option in the tag dropdown menu next to the Image. The steps are the same for the other containers listed above, searching for the relevant Image and tag names.

<!-- <img src="./../figures/docker_search.png" alt="Docker - searching for the required Docker image and Tag" width="100%"/> -->
<img src="docker_search.png" alt="Docker - searching for the required Docker image and Tag" width="100%"/>
<a id='figure_1'></a> 
    <center><b>Figure 2. Docker search for the relevant Image and tag.</b></center>

Press `Run` to pull and run the container at the same time. A popup window should appear once everything is downloaded, as in Figure 3. The download time varies with the size of the container and internet speed and may take some time.
> If the popup window doesn't appear, navigate to `Images` in the left menu and initialise the container by pressing triangle by the relevant image. 

Expand <span style="color:red"> `Optional settings`</span> and fill the fields as shown in the following Figure 3. These settings set the ports for ssh and browser to connect to the JupyterHub container.

<!-- <img src="./../figures/docker_container_run_settings.png" alt="Docker - container run settings for ports: 2222:22 and 8000:8000" width="100%"/> -->
<img src="docker_container_run_settings.png" alt="Docker - container run settings for ports: 2222:22 and 8000:8000" width="100%"/>
<a id='figure_1'></a> 
    <center><b>Figure 3. Docker container run settings - select Optional settings and fill in as above.</b></center>

After filling in the fields as above, press `Run`. The installed container can be found in the `Containers` section in the menu on the left. To open the container in your browser click on the second option in the `Port(s)` column of the required container or go to [http://localhost:8000](http://localhost:8000).

<!-- <img src="./../figures/docker_container_run.png" alt="Docker - containers stopping and running" width="100%"/> -->
<img src="docker_container_run.png" alt="Docker - containers stopping and running" width="100%"/>
<a id='figure_1'></a> 
    <center><b>Figure 4. Stopping and starting a docker container.</b></center>

If you want to stop the container, you can do so from the `Actions` column on the right. Square to stop, and triangle to start it again. If you used this method, you can skip the next section and continue to [Accessing the hub section](#docker_access_hub) to find login options and how to run notebooks.

## Using terminal {#docker_terminal_installation}

For this section, initialise Docker, open your terminal and use the following commands.

1. Pull image
~~~~
docker pull likask/mofem-spack-jupyterhub:Workshop2023
~~~~
replace `likask/mofem-spack-jupyterhub:Workshop2023` required by `image_name:tag_name` taken from [table listing images and tags](#docker_getting_container).

2. Run container
- If you would like to just try the container, removing the container after use, run docker as follows:
~~~~
docker run --rm --name workshop2023 -p 8000:8000 -p 2222:22 likask/mofem-spack-jupyterhub:Workshop2023
~~~~
- If you would like to switch it on for some time, we recommend running it as a daemon:
~~~~~
docker run -d --name workshop2023 -p 8000:8000 -p 2222:22 likask/mofem-spack-jupyterhub:Workshop2023
~~~~~
Once installed, you do not have to reinstall it. Instead, start it again by:
~~~~~
docker start workshop2023
~~~~~

#### ARM architecture case on Mac

If you have a Mac with an ARM chip, you have to switch platforms when you run the compiler,
~~~~~~
docker run -d --platform linux/amd64 --name workshop2023 -p 8000:8000 -p 2222:22 likask/mofem-spack-jupyterhub:Workshop2023
~~~~~~
That results in a suboptimal performance, however, it is a workable solution. 

> The base system is Ubuntu 20.04. To compile code for *arm* architecture, we would have to upgrade the system to Ubuntu 22.04, and then it would be possible to compile MoFEM ecosystem for M1 chip. That is tested and works. However, additionally, you would have to compile gMesh from scratch. Python pip installation for gMesh and *arm* architectures is not available. If you know how to do it, we will welcome PR from you to fix this problem.

# Accessing and running with JupyterHub {#docker_access_hub}

## Password and login {#docker_password_login}

If you run a container locally, open [http://localhost:8000](http://localhost:8000) in your browser

- The default login name is *mofem*
- On the first login, the first password you input will be you password from then onwards. 
Note this is the password to JupyterHub, not a password to the Linux environment.

You have admin rights, and you can add more users.

## Start running 

- **Important**. Before you start, execute *install.md* notebook. It will copy symbolic links to the executable binaries of MoFEM installation to your directory.

- Navigate to any of the Jupyter notebooks and more instructions should be included inside. 

#### Being a good citizen

This is a case when the container is running on a server, and you share resources with other users.

- If you run something with multiprocess which will run longer than 5-10 minutes, be nice, i.e. run the command as follows
~~~~
nice -n 10 mpirun -np 2 ./command_line
~~~~

# Development and debugging setup  

If you wish to develop and debug in MoFEM using the just created containers, read this section. Otherwise, you are done, and you can ignore the following instructions :) 


#### Password for debugging 
First of all, set up your password. This is different on from the one you use to access JupyterHub in the browser. To change it, run the following command from a terminal:

~~~~
docker exec -it workshop2023 /bin/bash
~~~~
where <span style="color:green">`workshop2023`</span> is the name of your container set in `Optional settings` and can be found in the `Name` column of the Containers section. To change the password use: 
~~~~
passwd mofem
~~~~
and afterwards choose a new password.

### Connecting to the container

We will only cover the process with using [Visual Studio Code](https://code.visualstudio.com). Start by making sure you have it installed. 

Next, go to `Extensions` within VS Code and install `Remote Explorer` extension. This will allow you to connect to servers or containers like the one we created. While here, also install `C/C++ Extension Pack` extension.

Copy the following code into ~/.ssh/config on your laptop or desktop:
~~~~
Host workshop2023
  HostName localhost
  ForwardX11 yes
  Compression yes
  User mofem
  Port 2222
~~~~
Here, we can see the port forwarding number 2222 we set up for our container earlier. The option to SSH to the container should now appear in the `Remote Explorer` extension on VS Code. Connect to the container in VS Code (use the password for debugging) and you should see the same things as you do when you open it through a browser ([http://localhost:8000](http://localhost:8000)). 

The next step can be done on either platform. Open `install_from_source` notebook and run all of the cells in it. This will create a separate version of MoFEM for the user `mofem`. New folders will appear in your starting directory.

- mofem_install - contains source code
- um_view_debug - symbolic links to the executable binaries of the debugging version

## VS Code debugging setup

To set up debugging, follow these steps:

- open `mofem_install/mofem-cephas/mofem/users_modules` folder through `File -> Open Folder...` 
- create `.vscode` folder
- create [`launch.json`](#launchjson) and [`tasks.json`](#tasksjson) files
- copy the code from the following sections into these files respectively
- replace the hash `5sehreo` with the hash in your folder
- adjust the files to fit your purpose as described bellow 

#### launch.json {#launchjson}
~~~
{
    "version": "0.2.0",
    "configurations": [
        {
            "type": "cppdbg",
            "request": "launch",
            "name": "nonlinear_elastic",
            "program": "/mofem_install/jupyter/mofem/mofem_install/mofem-cephas/mofem/users_modules/um-build-Debug-5sehreo/tutorials/vec-2/nonlinear_elastic",
            "args": [
                "-file_name",
                "beam_3D.cub",
                "-order",
                "2",
            ], // , "-log_quiet"
            "cwd": "/mofem_install/jupyter/mofem/mofem_install/mofem-cephas/mofem/users_modules/um-build-Debug-5sehreo/tutorials/vec-2/",
            "preLaunchTask": "build_vec_2",
            "environment": [],
            "externalConsole": false,
            "MIMode": "gdb",
            "setupCommands": [
                {
                    "description": "Enable pretty-printing for gdb",
                    "text": "-enable-pretty-printing",
                    "ignoreFailures": true
                },
                {
                    "description": "Set Disassembly Flavor to Intel",
                    "text": "-gdb-set disassembly-flavor intel",
                    "ignoreFailures": true
                }
            ]
            // "stopOnEntry": true
        },
    ]
}

~~~

`/mofem_install/jupyter/mofem/mofem_install/mofem-cephas/mofem/users_modules/um-build-Debug-5sehreo/tutorials/vec-2/` is the folder where the executable is and you can change it to what you want to debug. `nonlinear_elastic` is the executable name. I would suggest copying and then searching for these terms to change them in all of the locations. Lastly, check what mesh is located at the new location (wanted by `-file_name`).

#### tasks.json {#tasksjson}

This file defines which folder should be rebuild before running your code with debugging.

~~~
{
    "tasks": [
        {
            "label": "build_vec_2",
            "type": "shell",
            "command": "make -j6 -C /mofem_install/jupyter/mofem/mofem_install/mofem-cephas/mofem/users_modules/um-build-Debug-5sehreo/tutorials/vec-2",
        },
    ],
    "version": "2.0.0"
}
~~~

## Run and debug

- Go to the file you want to debug and place a breakpoint inside one of the functions.
- `Run -> Start Debugging` from the dropdown menu to start debugging

The program should now start and eventually stop at the breakpoint if everything was set correctly. For more information see [VS Code instructions for debugging](#https://code.visualstudio.com/Docs/editor/debugging).

The same procedure can be applied for any other users you create and everyone can debug separately. 



<!-- Use *workshop2023* when you log in through VS Code. Note that you are connecting to *jupyterhub cloud/docker container*. -->

## Video on JuputerHub, SSH and MoFEM

[![Watch the video](https://img.youtube.com/vi/xL3J8VHig68/hqdefault.jpg)](https://youtu.be/xL3J8VHig68)

Any problems with this installation, please contact us by [mofem-group@googlegroups.com](https://groups.google.com/forum/#!forum/mofem-group).