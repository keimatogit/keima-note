# SwitchyOmegaでサーバーでjupyterを使う

PROFILES -> New Profile ... name: SOCKS_1080
Proxy servers
Protocol: SOCKS5
Server: localhost
Port: 1080
Bypass Listはそのままでいい

PROFILES -> auto switch
下２つを登録
----
Condition Type: Host wildcard
Condition Details: idnfront1
profile: SOCKS_1080
----
Condition Type: Host wildcard
Condition Details: idnsgiuv
profile: SOCKS_1080

```
ssh -f -N -D 1080 idnfront1
ssh  idnfront1
qsub -I
jupyter notebook --ip=* --no-browser
```