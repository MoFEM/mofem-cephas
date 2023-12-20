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

```python
import os
import zipfile
zipf = zipfile.ZipFile('vtks_and_pngs.zip', 'w', zipfile.ZIP_DEFLATED)
path = "."
for root, dirs, files in os.walk(path):
    for file in files:
        if file.endswith('.vtk') or file.endswith('.png'):        
            zipf.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), os.path.join(path, '..')))
zipf.close()
```
