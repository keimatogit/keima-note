# ssh接続でホスト名を登録

[【ssh接続】ホスト名を設定してIPアドレスを省略する設定方法](https://tofusystem.work/programming-hack/easy-ssh-login-setting/)

`~/.ssh/config`を作成して設定を書き込む。複数可。

```
# ~/.ssh/config

Host  host_name1（好きな名前）
  Hostname  IP_address1
  User  user_name
  IdentityFile  private_key_file_path

Host  host_name2
  Hostname  IP_address2
  User  user_name
  IdentityFile  private_key_file_path2
```

接続するとき
```
ssh host_name1
```
