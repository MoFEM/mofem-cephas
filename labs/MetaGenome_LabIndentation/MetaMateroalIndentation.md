---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.5
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# Indentation Test
![Example3dAnim.gif](Example3dAnim.gif)
![meta-genome-logo.png](meta-genome-logo.png)


# MoFEM eccosystem
MoFEM has complex software eccosystem composed with many libraries. It scalable can be run on laptops with limited resurcers, up to large coumputer clusters.

<img src="https://bitbucket.org/likask/mofem-cephas/raw/2bba8d8fc3a7df1b9d976fbffee9a61e5cff64b5/media/Ecosystem/Ecosystem.002.png" alt="Ecosystem" width="800"/>

MoFEM is designed to run complex multi-physics problems. It supports research, it generated impact case for REF2021, and is used by indsutry. 

<td>
<img src="https://bitbucket.org/likask/mofem-cephas/raw/2bba8d8fc3a7df1b9d976fbffee9a61e5cff64b5/media/Ecosystem/Ecosystem.003.png" alt="Ecosystem" width="800"/>
</td>
</tr></table>


# Set problem parameters

Mesh flie geometry can be made in [gMesh](https://gmsh.info), [Salome](https://www.salome-platform.org)  or [Coreform](https://coreform.com). 

gMesh, Salome are free and open sourve, Coreform is free to use for research.

Note: MoFEM reads Cubit files, not Cube HDF5 files.

```python
from IPython.display import HTML

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from pyvirtualdisplay import Display
display = Display(backend="xvfb", visible=False, size=(800, 600))
display.start()

plt.rcParams['figure.figsize'] = [15, 10]

HTML("<style>.container { width:100% !important; }</style>")
```

```python
# mesh dimension
dimension=2

# material parameters 
young_modulus=1
poisson_ratio=0.3

# mesh file
mesh_file='anit-chiral2d.cub'
# output
time_depths_factor=1

# surface distance function
import importlib
import sdf
importlib.reload(sdf)
print('radius = %f' % sdf.r)
```

```python
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

import importlib
import sdf
importlib.reload(sdf)
print('radius = %f' % sdf.r)

mesh_file_vtk=mesh_file[:-3]+'vtk'
!/mofem_install/um_view/bin/mbconvert {mesh_file} {mesh_file[:-3]}vtk

mesh = pv.read(mesh_file_vtk)
mesh.point_data.keys()
mesh.point_data.remove('GLOBAL_ID')
mesh.point_data['xyz']=mesh.points

sphere = pv.Sphere(radius=sdf.r, center=(sdf.xc, sdf.yc, sdf.zc), phi_resolution=45, theta_resolution=45)

p = pv.Plotter(notebook=True, window_size=[500, 800])
mesh = mesh.shrink(0.9)
p.add_mesh(mesh, scalars='xyz')
p.add_mesh(sphere, opacity=0.50)
if dimension == 2:
    p.camera_position = "xy"
p.show(jupyter_backend='panel')

#print(mesh.point_data)
#mesh = mesh.shrink(0.95)
#mesh.plot(smooth_shading=True, scalars='xyz', notebook=True, jupyter_backend='ipygany')  
```

```python
wd=!pwd
wd=wd[0]

if dimension == 3:
    code='/mofem_install/um_view/tutorials/adv-1/contact_3d'
elif dimension == 2:
    code='/mofem_install/um_view/tutorials/adv-1/contact_2d'
else:
    raise Exception("Only works for 2d or 3d")
    
print(wd)
print(code)

#number of processors
nb_proc=4

# time stepping
time=30
time_step=0.25

# computational parameters
cn=1e-2/young_modulus
spring_stiffness=0
alpha_dumping=0.0001 # should be small value, used for dynamic relaxation
!/mofem_install/um_view/bin/mofem_part -file_name {mesh_file} -my_nparts {nb_proc} -dim {dimension} -adj_dim 0 -output_file one.h5m
```

```python
# time stepping
time=10
time_step=0.1

# computational parameters
cn=1e-2/young_modulus
spring_stiffness=0
alpha_dumping=0.0001 # should be small value, used for dynamic relaxation

log_file='log'

!rm out*h5m
!export OMPI_MCA_btl_vader_single_copy_mechanism=none && \
nice -n 5 mpirun --allow-run-as-root -np {nb_proc} {code} -file_name one.h5m \
-ts_type beuler \
-ts_max_time {time} \
-ts_dt {time_step} \
-log_no_color \
-snes_linesearch_type basic \
-cn {cn} \
-spring_stiffness {spring_stiffness} \
-alpha_dumping {alpha_dumping} \
-young_modulus 1 \
-poisson_ratio {poisson_ratio} -order 3 2>&1 | tee {log_file}
```

```python
# Plot convergence
!sed 's/\[0\] <inform> \[petsc\]//g' {log_file} > log_snes
!grep Function log_snes | sed 's/\[//g' | sed 's/,//g' > snes
newton_data=pd.read_fwf('snes', header=None)
newton_data=newton_data.rename(columns={0: "it", 4: "res", 5: "equlibrium", 6: 'constrain'})
newton_data=newton_data.drop([1, 2, 3], axis=1)
print(newton_data)

# newton_data

plt.plot(newton_data['res'].to_numpy(),'r^-')
plt.title('Neton method convergence')
plt.ylabel('absolute residial')
plt.xlabel('accumulated iterations')
plt.yscale('log')
plt.grid(True)
plt.show()
```

```python
# converting analysis files to format readable by post processors

!rm -f *.vtk

import re, os
list_of_files=!ls -c1 out_*.h5m
def extract_numner(s):
    return int(re.findall(r'\d+',s)[0])

size=len(list_of_files)
mod=int(size / 10)

list_to_process=[]
for f in list_of_files:
    n=extract_numner(f)
    if n == 0 or n == 1:
        list_to_process.append(f)
    elif n % mod == 0:
        list_to_process.append(f)
        
sorted_list_of_files = sorted(list_to_process, key=extract_numner)

out_to_vtk = !ls -c1 out_*h5m
last_file=out_to_vtk[0]
print(last_file)
!/mofem_install/um_view/bin/mbconvert {last_file} {last_file[:-3]}vtk
# for i in sorted_list_of_files:
#     !mbconvert {i} {i[:-3]}vtk

scale_displacements=1

list_of_files=!ls -c1 out*.vtk
list_of_files=sorted(list_of_files, key=extract_numner)
print(list_of_files)

mesh = pv.read(list_of_files[len(list_of_files)-1])
print(mesh)

mesh=mesh.warp_by_vector('U',factor=scale_displacements)

#geom = pv.Arrow()  # This could be any dataset
#glyphs = mesh.glyph(orient="t", scale="t", factor=0.003, geom=geom)

sphere = pv.Sphere(radius=sdf.r, center=(sdf.xc, sdf.yc, sdf.zc), phi_resolution=45, theta_resolution=45)

field='U'
my_cmap = plt.colormaps["rainbow"]
p = pv.Plotter(notebook=True, window_size=[500, 800])
p.add_mesh(mesh, scalars=field, show_edges=True)
p.add_mesh(sphere, opacity=0.25)
if dimension == 2:
    p.camera_position = "xy"
#p.add_mesh(glyphs, opacity=0.25)
p.show(jupyter_backend='panel')

```

```python
# Plot load disp-path
!sed 's/\[0\] <inform> \[indent\]//g' {log_file} > log_path
!grep Force log_path > reaction
!grep Ux log_path > ux
!grep Uy log_path > uy
!grep Uz log_path > uz

data_reaction=pd.read_csv('reaction',sep='\s+',header=None)
data_reaction=data_reaction.rename(columns={2: "time", 3: "tx", 4: "ty", 5: "tz"})
data_reaction=data_reaction.drop([0, 1], axis=1)
#print(data_reaction)

def get_u_data(str):
    data=pd.read_csv(str,sep='\s+',header=None)
    data=data.rename(columns={2: "time", 4: "min", 6: "max"})
    data=data.drop([0, 1, 3, 5], axis=1)
    return data

data_ux=get_u_data('ux')
data_uy=get_u_data('uy')
if dimension == 3:
    data_uz=get_u_data('uz')
#print(data_ux)

fig, (ax1,ax2) = plt.subplots(2,1)

yt=young_modulus
tf=time_depths_factor

ax1.plot(tf*data_reaction['time'].to_numpy(), yt*data_reaction['tx'].to_numpy(), 'ro-', label='tx')
ax1.plot(tf*data_reaction['time'].to_numpy(), yt*data_reaction['ty'].to_numpy(), 'g^-', label='ty')
if dimension == 3:
    ax1.plot(tf*data_reaction['time'].to_numpy(), yt*data_reaction['tz'].to_numpy(), 'b+-', label='tz')
legend = ax1.legend(loc='upper left', shadow=True)
ax1.set(xlabel='indentation depth', ylabel='force',
        title='Load path indentation depth vs force')
ax1.grid(True)

ax2.plot(data_uy['min'].to_numpy(), yt*data_reaction['tx'].to_numpy(), 'ro-', label='min(Uy)-tx')
ax2.plot(data_uy['min'].to_numpy(), yt*data_reaction['ty'].to_numpy(), 'g^-', label='min(Uy)-ty')
if dimension == 3:
    ax2.plot(data_uy['min'].to_numpy(), yt*data_reaction['tz'].to_numpy(), 'b+-', label='min(Uy)-tz')
ax2.plot(data_uy['max'].to_numpy(), yt*data_reaction['tx'].to_numpy(), 'ro-', label='max(Uy)-tx')
ax2.plot(data_uy['max'].to_numpy(), yt*data_reaction['ty'].to_numpy(), 'g^-', label='max(Uy)-ty')
if dimension == 3:
    ax2.plot(data_uy['max'].to_numpy(), yt*data_reaction['tz'].to_numpy(), 'b+-', label='max(Uy)-tz')
legend = ax2.legend(loc='upper left', shadow=True)
ax2.set(xlabel='displacement', ylabel='force',
        title='Load displacement path min[uy] vs max[uy]')
ax2.grid(True)
```
