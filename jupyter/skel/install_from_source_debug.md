---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.16.0
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

This notebook is for installing a Debug version of MoFEM.

> Debug version is slower and mainly used for code developent

Clone MoFEM core and users modules.

```bash
mkdir -p mofem_install
cd mofem_install
git clone -b Workshop2023 --recurse-submodules https://bitbucket.org/likask/mofem-cephas.git
cd $HOME/mofem_install/mofem-cephas/mofem/users_modules
git checkout Workshop2023
```

Print User Id

```bash
id $UID
```

Install MoFEM core - Debug version

```bash
export TARGET=x86_64
spack --config-scope /mofem_install/spack_config_dir dev-build -j4 \
  --source-path $HOME/mofem_install/mofem-cephas \
  --keep-prefix \
  mofem-cephas@lukasz~copy_user_modules \
  target=$TARGET build_type=Debug install_id=$UID ^petsc+X ^boost+python+numpy
```

Check existing installations of `mofem-cephas`

```bash
export TARGET=x86_64 
spack --config-scope /mofem_install/spack_config_dir find -lv mofem-cephas install_id=$UID
```

Install user modules including tutorials which you can learn from.

```bash
spack --config-scope /mofem_install/spack_config_dir dev-build -j4 \
  --source-path $HOME/mofem_install/mofem-cephas/mofem/users_modules \
  mofem-users-modules@lukasz build_type=Debug install_id=$UID  \
  ^mofem-cephas@lukasz+adol-c~copy_user_modules~docker~ipo+med+mgis~shared+slepc+tetgen build_system=cmake build_type=Debug dev_path=/mofem_install/jupyter/$USER/mofem_install/mofem-cephas install_id=$UID \
  ^petsc+X ^boost+python+numpy
```

Check existing installations of `mofem-users-modules`

```bash
export TARGET=x86_64 
spack --config-scope /mofem_install/spack_config_dir find -lv mofem-users-modules install_id=$UID
```

Create symbolic links to the Debug version of the code

```bash
rm -rf um_view_debug
spack view symlink -i um_view_debug mofem-users-modules build_type=Debug install_id=$UID
```
