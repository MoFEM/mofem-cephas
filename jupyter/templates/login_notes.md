### Password and login

- On the first login, select the password
- The shell password on the SSH login is different than the JupyterHub password
- Ask admin for SSH password
- Change password in the shell using [passwd](https://www.computerhope.com/unix/upasswor.htm) command
- Set up your shell with the command (e.g. bash):
```
$ usermod --shell /bin/bash your_login_name
```

### SSH config if you have account on cactus

Copy the following code into .ssh/config on your laptop or desktop.
```
Host forward_greenstick
  HostName cactus.eng.gla.ac.uk
  ForwardX11 yes
  Compression yes
  User YOUR_LOGIN_ON_HUB
  Port 2222
```

Use *forward_greenstick* when you login in VSCode. Note tha you are connecting to *cactus* from which is tunnel to *geenstick*.

