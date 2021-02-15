---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.10.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

```python
!whoami
!$MOFEM_ENV_FILE
!rm -rf um_view
!spack view symlink -i um_view mofem-softmech
```
