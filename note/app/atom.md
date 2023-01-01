# ATOM

[ATOM](https://atom.io/)からインストール

## 入れておきたいパッケージ

インストール方法：
Atom -> Preferences -> Install -> 検索してインストール

- platformio-ide-terminal: ターミナルを使う
- remote-ftp: サーバーのファイルを扱う
- language-r: Rスクリプトのハイライトを行う
- markdown-preview-enhanced: Markdownプレビュー（最初から入っているmarkdown-previewはdisableしておく）
- markdown-toc: Markdownに目次を入れる


### platformio-ide-terminal用の設定

Atom -> Keymapでkeymap.cson開き、
```
'atom-workspace atom-text-editor:not([mini])': 'cmd-enter': 'platformio-ide-terminal:insert-selected-text'
```
を書き込むと、cmd+enterでカーソル一の行や選択した部分をターミナルで実行できる。shifやcontrolが良ければ`cmd-enter`の代わりに`shift-enter`や`ctrl-enter`にする。

windowsでubuntuを使いたいのにコマンドプロンプトやpowershellが起動してしまう場合、`wsl`と叩くとでデフォルトのディストリビューションが起動する（ubuntuにしておこう）。platformio-ide-terminalのsettingのCore -> Auto Run Commandにwslを入れて、起動毎に自動実行するようにすると良い。


### remote-ftpの使い方

1. 接続用のフォルダを用意する（`~/remote/サーバー名`など）
2. File -> Add Project Folderで接続用のProjectフォルダを選択
3. Packages -> Remote FTP -> Create FTP/SFTP config file
4. Projectフォルダ直下に作成された.ftpconfigを編集
5. Packages -> Remote FTP -> Toggle -> Remoteタブを選択 -> Connect

.ftpconfigの例。あらかた自動で書いてくれるので、host名やuser名やremote（サーバー内の場所）を書き込もう。
```
{
    "protocol": "sftp",
    "host": "idnfront1",
    "port": 22,
    "user": "user_name",
    "pass": "",
    "promptForPass": false,
    "remote": "/user/user_name",
    "local": "",
    "agent": "",
    "privatekey": "",
    "passphrase": "",
    "hosthash": "",
    "ignorehost": true,
    "connTimeout": 10000,
    "keepalive": 10000,
    "keyboardInteractive": false,
    "keyboardInteractiveForPass": false,
    "remoteCommand": "",
    "remoteShell": "",
    "watch": [],
    "watchTimeout": 500
}
```

他のサーバーや、サーバー内の別フォルダに接続して使いたい場合は、それぞれの接続用フォルダを用意して各フォルダ直下に.ftpconfigを作成するのが良い。

サーバーでの更新内容が反映されていないときは、Remoteタブで更新したいフォルダの上で右クリック -> Refreshで更新される

## その他やっておきたい設定

パッケージの無効化（Atom -> Preferences -> Packages -> 検索してDisable）
- markdown-preview： markdown-preview-enhancedを使う場合は無効化しておこう
- spell-check: 赤線がたくさん出て嫌なので無効化しておこう
