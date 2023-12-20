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

This notebook is for installing a Release version of what you might have changed during your exercises with debugging.

> Do not run this notebook if you have not run install_from_source_Debug

If you have run this notebook before and just want to recompile the latest changes, run the last cell only.

Print User Id

```bash
id $UID
```

Install MoFEM core - Release version

```bash
export TARGET=x86_64
spack --config-scope /mofem_install/spack_config_dir dev-build -j4 \
  --source-path $HOME/mofem_install/mofem-cephas \
  --keep-prefix \
  mofem-cephas@lukasz~copy_user_modules \
  target=$TARGET build_type=Release install_id=$UID
```

Check existing installations of mofem-cephas

```bash
export TARGET=x86_64 
spack --config-scope /mofem_install/spack_config_dir find -lv mofem-cephas install_id=$UID
```

Install updated user modules which you have ammended.

```bash
spack --config-scope /mofem_install/spack_config_dir dev-build -j4 \
  --source-path $HOME/mofem_install/mofem-cephas/mofem/users_modules \
  mofem-users-modules@lukasz build_type=Release install_id=$UID  \
  ^mofem-cephas@lukasz+adol-c~copy_user_modules~docker~ipo+med+mgis~shared+slepc+tetgen build_system=cmake build_type=Release dev_path=/mofem_install/jupyter/mofem/mofem_install/mofem-cephas install_id=$UID
```

Check existing installations of mofem-users-modules

```bash
export TARGET=x86_64 
spack --config-scope /mofem_install/spack_config_dir find -lv mofem-users-modules install_id=$UID
```

Create symbolic links to the Release version of the code

```bash
rm -rf um_view_release
spack view symlink -i um_view_release mofem-users-modules build_type=Release install_id=$UID
```

To run the notebooks provided with the amended version of the code, replace all of the `um_view` in the path definitions in the notebooks you want to run with `um_view_release`.


# Cell to run to update changes when debugging
(assuming you run all of the cells once before)

```bash
cd mofem_install/mofem-cephas/mofem/users_modules/um-build-Release-* && make -j4 install
```
