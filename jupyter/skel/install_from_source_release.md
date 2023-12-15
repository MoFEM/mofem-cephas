---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.15.1
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

This notebook assumes that you have already run install_from_source_debug.md

Print User Id

```bash
id $UID
```

Install users modules
> Release version is faster but you cannot run debugger with it.

```bash
spack --config-scope /mofem_install/spack_config_dir dev-build -j 4 \
  --test root  \
  --source-path $HOME/mofem_install/mofem-cephas/mofem/users_modules \
  mofem-users-modules@lukasz build_type=Release install_id=$UID  \
  ^mofem-cephas@lukasz+adol-c~copy_user_modules~docker~ipo+med~shared+slepc+tetgen build_system=cmake build_type=Release dev_path=/mofem_install/mofem-cephas install_id=0 \
  ^boost+python ^petsc+X
```

Create symbolic links to generic workshop installation, and user version

```bash
spack find -lv mofem-users-modules build_type=Release install_id=$UID
```

```bash
rm -rf um_view_release
spack view symlink -i um_view_release mofem-users-modules build_type=Release install_id=$UID
```

```python

```
