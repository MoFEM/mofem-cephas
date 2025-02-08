# Command: `make notebook`
  
The `make notebook` command is used to copy markdown files (excluding README.md) from specified directories within the source **tutorials** folder to the corresponding directories in the build **tutorials** folder. Once the files are copied, they are converted to Jupyter notebooks using the `jupytext` tool. This command is automatically executed when running `make`.

To specify the subdirectories for this functionality, include the `copy_and_convert_md_to_ipynb` function in their respective **CMakeLists.txt** file. For detailed instructions on how to include this functionality, refer to the **scl-8/README.md** file.

# How to run
To run the command, navigate to the build **tutorials** directory and execute the following command: `make notebook`.

# Installing jupytext 
If jupytext isn't found, `make notebook` command will not work. To install `jupytext` run, 
`pip install jupytext`.