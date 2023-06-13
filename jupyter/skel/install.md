---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.1
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

```python
!whoami
!$MOFEM_ENV_FILE
!rm -rf um_view
!spack view symlink -i um_view mofem-users-modules
#!rm -f adv-0 && ln -s um_view/tutorials/adv-0 .
#!rm -f vec-1 && ln -s um_view/tutorials/vec-1 .
#!rm -f scl-8 && ln -s um_view/tutorials/scl-8 .
#!rm -f nonlinear_elasticity && ln -s um_view/nonlinear_elasticity .
#!rm -f unsaturated_2Dflow && um_view/softmech/unsaturated_2Dflow .
```

```python
|
```
