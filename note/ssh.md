# sshの設定

```
ssh user_name@host_name
```

初めて接続するときはknown_hostに登録するか聞かれるので、登録しておこう（`~/.ssh/known_host`に書き込まれる）。

```
Host	idnfront
	HostName	133.1.178.1
	User	user_name
	IdentityFile	~/.ssh/id_rsa

Host *
	ServerAliveInterval 300
```

鍵の登録（パスワードが不要になる）
ssh-keygen # id_rsa（秘密鍵）とid_rsa.pub（公開鍵）を作成。パスフレーズは何も入力しなくてOK。
ssh username@idnfront1 "cat >> ~/.ssh/authorized_keys" < ~/.ssh/id_rsa.pub
~/.ssh/authorized_keys (サーバ側のファイル．ログインされることを許す公開鍵一覧)




```shell
cd ~/.ssh
ssh-keygen # key作成（パスフレーズを聞かれるけど何も入力しなくてOK）


ssh username@idnfront1 "cat >> ~/.ssh/authorized_keys" < ~/.ssh/id_rsa.pub
~/.ssh/id_rsa (秘密鍵．ssh-keygenにより生成)
~/.ssh/id_rsa.pub (公開鍵．ssh-keygenにより生成)
~/.ssh/authorized_keys (サーバ側のファイル．ログインされることを許す公開鍵一覧)
```
