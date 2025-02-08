
# Command: make notebook_scl-8
  
The `make notebook_scl-8` command is similar to the `make notebook` command, but it works locally for the scl-8 tutorial. When executed, this command will copy markdown files (excluding README.md) from the source directory **scl-8** to the build directory **scl-8** and convert them to Jupyter Notebooks. This command relies on the `make notebook_pluto` target, where ***pluto*** represents the markdown filename. For multiple markdown files, each file will have an associated `make notebook_filename_xyz` command, and `make notebook_scl-8` will be a collection of these commands.

## Modifying CMakeLists.txt
To enable the `make notebook_scl-8` functionality, add the following line to the **CMakeLists.txt** file:
```
copy_and_convert_md_to_ipynb(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
```

## How to execute
Navigate to the **scl-8** directory and run `make notebook_scl-8`.
