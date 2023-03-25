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

```bash
rm -f out*
rm -f one.h5m
rm -f log*
rm -f snes reaction ux uy uz
rm -f *ipynb
```

```bash
apt-get update
apt-get install -y zip unzip
```

```bash
zip all *  -x "__pycache__*"
```
