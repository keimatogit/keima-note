---
html:
  embed_local_images: false
  embed_svg: true
  offline: true
  toc: true
toc:
  depth_from: 2
  depth_to: 3
  ordered: true
export_on_save:
  html: true
---

# Visual Studio Code

[Visual Studio Code](https://code.visualstudio.com/)からインストール

## カラーテーマの設定

Code -> Preferences -> Color Theme
お気に入りはAtom One Dark Theme（Mahmoud Ali）


## リモート接続

Extension: Remote Development (Microsoft)を追加。左側のメニューにRemotte Explorerが表示される。上部のプルダウンメニューでSSH Targetを選ぶと、~/.ssh/configで設定した接続先が表示される。接続ができたら「Open Folder」から接続先のフォルダを開いて操作できるようになる。

ssh先で使いたいExtensionは、接続先に入れないといけないみたい。

ubuntuにつなぐ場合も同じパッケージでOK（プルダウンでWSL Targetsを選ぶ）:[Linux 用 Windows サブシステムで Visual Studio Code の使用を開始する](https://docs.microsoft.com/ja-jp/windows/wsl/tutorials/wsl-vscode)


## cmd+enterでターミナルで実行できるようにする

[VS Code execute current line or selection to in the integrated console](https://stackoverflow.com/questions/45667252/vs-code-execute-current-line-or-selection-to-in-the-integrated-console)

AtomやRstudioと同じように、cmd+enterで現在の行を実行して行送りor選択部分の実行を行えるようにする。pythonやRのスクリプトでも、ターミナルの方でpythonやRを起動しておけば実行できるので便利だと思う。

1) Extension: macrosRe (l7ssha)を追加。（Extension: macros (geddski)は想定通りに動かない時がある）

2) settings.jsonに追加（command + shift + P -> JSONで検索 -> “Open Settings (JSON)”）（一番外側の{}内に入れる）
```
"macros":{
        "execCurLn": [
            "expandLineSelection",
            "workbench.action.terminal.runSelectedText",
            "cancelSelection",
        ]
},
```

3) keybindings.jsonに追加（command + shift + P -> JSONで検索 -> “Open Keyboard Shortcuts (JSON)”）（一番外側の[]内に入れる）
```
    {
        "key": "cmd+enter",
          "command": "workbench.action.terminal.runSelectedText",
            "when": "editorTextFocus && editorHasSelection"
    },
    {
        "key": "cmd+enter",
          "command": "macros.execCurLn",
            "when": "editorTextFocus && !editorHasSelection"
    },
```

## jupyter Notebook

[Jupyter Notebooks in VS Code](https://code.visualstudio.com/docs/datascience/jupyter-notebooks)

- VSCodeでマイクロソフトのJupyterパッケージをインストールする。
- Command + p -> Python: Select Interpreterを選び、使いたい環境を選ぶ。

### Python 

Extension: Python

Command + Shift + P -> "select interpreter" -> 使いたいpythonを選ぶ

デフォルトだとshift + enterで実行（行送りなし）。command + enterで行送りもつける方法は以下。
[Running code in VS Code with Ctrl + Enter: move to next line](https://actuarialdatascience.com/shortcut_vscode.html)
[VS Code execute current line or selection to in the integrated console](https://stackoverflow.com/questions/45667252/vs-code-execute-current-line-or-selection-to-in-the-integrated-console)


settings.jsonに追加（command + shift + P -> JSON -> “Open Settings (JSON)”）（一番外側の{}内に入れる）
```
"macros": {
    "pythonExecSelectionAndCursorDown": [
        "python.execSelectionInTerminal",
        "cursorDown",
    ]
} ,
```

keybindings.jsonに追加（command + shift + P -> JSONで検索 -> “Open Keyboard Shortcuts (JSON)”）（一番外側の[]内に入れる）
```
    {
      "key": "cmd+enter",
      "command": "macros.pythonExecSelectionAndCursorDown",
      "when": "editorTextFocus && editorLangId == 'python'"
    }
```





### R

[R in Visual Studio Code](https://code.visualstudio.com/docs/languages/r)
(https://i-doctor.sakura.ne.jp/font/?p=44103#toc_id_4)
[vscodeでRをいい感じで使えるようにする](https://qiita.com/aishidajt9/items/319659c9344bb7ae4302)


1) Rのパッケージlanguageserverとjsonliteを入れておく
```
conda install -c r r-jsonlite
conda install -c conda-forge r-languageserver
```

2) R extension for VScode をインストール

3) 使用するRのパスを設定：
Code -> prefaces -> settings -> Extensions -> R
R > Rpath: MacとR > Rterm: Mac
`/Users/user1/miniconda3/envs/r-env/bin/R`など。


- デフォルトでcommand + enterで実行&行送りができる。選択してcommand + enterで選択部分を全部実行できる。
- 構文チェックが鬱陶しいときは、R>Lsp:Diagnosticのチェックをはずす。


### Jupyter

[Jupyter Notebooks in VS Code](https://code.visualstudio.com/docs/datascience/jupyter-notebooks)


## 色々な設定

Code -> Preferences -> Settings

設定はユーザー単位とプロジェクト単位で指定できます。
ユーザー単位の設定ファイル：/Users/[USER_NAME]/Library/Application\ Support/Code/User/settings.json
プロジェクト設定ファイル：[PROJECT_FOLDER]/.vscode/settings.json

マイ設定

Editor: Font Family
`'Ricty Diminished', Menlo, Monaco, 'Courier New', monospace`
[Ricty Diminished](https://github.com/edihbrandon/RictyDiminished)
(Code -> Download ZIP -> 解凍 -> ダブルクリックで開く -> フォントをインストール)

Editor: Font Size
`14`

Terminal › Integrated: Font Size
`14`

## Markdown関連

[Visual Studio Code で Markdown 編集環境を整える](https://qiita.com/kumapo0313/items/a59df3d74a7eaaaf3137)

### Markdown Preview Enhanced (Yiyi Wang)

[Markdown Preview Enhanced](https://shd101wyy.github.io/markdown-preview-enhanced/#/)

#### プレビューの表示

右クリック -> Markdown Preview Enhanced: Open Preview to the Side、または`Cmd(Ctl) + K` -> `V`または右上にアイコンが出ているかも。

プレビューのスタイルを変更するには、プレビュー画面で右クリック -> Preview Theme -> 好きなテーマを選ぶ

#### TOCの作成

cmd-shift-p -> choose Markdown Preview Enhanced: Create Toc

- `orderedList` Use orderedList or not.
- `depthFrom`, `depthTo` [1~6] inclusive.
- `ignoreLink` If set to true, then TOC entry will not be hyperlinks.

TOCに含めたくない場合はheadingの後ろに` {ignore=true}`を追記する

#### html出力

プレビュー画面で右クリック -> HTML -> HTML(offline)

html出力時の設定はフロントマターで。

デフォルト（toc: trueでサイドバーtocを表示）

```
---
html:
  embed_local_images: false
  embed_svg: true
  offline: false
  toc: undefined

print_background: false
---
```

tocの設定も書ける。

```
---
toc:
  depth_from: 2
  depth_to: 3
  ordered: true
---
```

以下のフロントマターを追加すると、保存のたびにhtml出力される。cssの変更などは拾ってくれていない気がするので、ヘンだったらプレビュー画面で更新してからから右クリックでHTML出力しよう。

```
---
export_on_save:
  html: true
---
```

フロントマターの例

```
---
html:
  embed_local_images: false
  embed_svg: true
  offline: true
  toc: true
toc:
  depth_from: 2
  depth_to: 3
  ordered: true
export_on_save:
  html: true
---
```

#### html 出力時にサイドバー TOC を有効化（セキュリティリスクがあるらしい）

settings.jsonに追加（command + shift + P -> JSONで検索 -> "Open User Settings (JSON)"）（一番外側の{}内に入れる）

```
"markdown-preview-enhanced.enableScriptExecution": true
```


三本線のメニューボタン位置を左下から左上にするには、style.less（ctrl+shift+p -> Markdown Preview Enhanced: Customize Css）に以下を追記。

```
// サイドバーTOCメニューボタンを左上にする
#sidebar-toc-btn {
  bottom: unset;
  top: 8px;
}
// TOCメニューボタンが左上に来た分、TOCの上部を空ける
.md-sidebar-toc.md-sidebar-toc {
  margin-top: 60px !important;
}

```

本体のTOCを消すには、TOC部分を`<div class="display-none"></div>`で囲んで、style.lessに`.display-none {display:none;}`を追加しておく。


#### カスタムCSS

cmd(ctr) + shift + P -> Markdown Preview Enhanced: Customize.cssでstyle.cssが開くので編集する。
style.lessの例

```
/* Please visit the URL below for more information: */
/*   https://shd101wyy.github.io/markdown-preview-enhanced/#/customize-css */
@font-face {
  font-family: "游ゴシック体",
    "YuGothic",
    "游ゴシック",
    "Yu Gothic",
    "ヒラギノ角ゴ Pro W3",
    "Hiragino Kaku Gothic Pro",
    sans-serif;
}

.markdown-preview.markdown-preview {
  // modify your style here
  // eg: background-color: blue;
  font-family: "游ゴシック体",
    "YuGothic",
    "游ゴシック",
    "Yu Gothic",
    "ヒラギノ角ゴ Pro W3",
    "Hiragino Kaku Gothic Pro",
    sans-serif;
  
  pre[class*="language-"]>code {
    font-size: 1em !important;
  }

  h1 {
    counter-reset: chapter;
    font-size: 1.8em;
    border-bottom: solid 1px;
    margin-bottom: 40px;
  }
  h2 {
    counter-reset: sub-chapter;
    font-size: 1.4em;
    margin-top: 60px;
    margin-bottom: 20px;
    padding: 0.4em 0.6em;
    background: #3182bd;
    color: white;
    line-height: 2;
  }
  h3 {
    counter-reset: section;
    font-size: 1.2em;
    margin-top: 40px;
    border-left: solid 10px #3182bd;
    padding-left: 10px;
    line-height: 1.8;
  }

  h2::before {
    counter-increment: chapter;
    content: counter(chapter) ". ";
  }
  h3::before {
    counter-increment: sub-chapter;
    content: counter(chapter) "." counter(sub-chapter) ". ";
  }
}

// メニューの三本線アイコンを左下から左上に
#sidebar-toc-btn {
  bottom: unset;
  top: 8px;
}

// 三本線アイコンの分上を空ける
.md-sidebar-toc.md-sidebar-toc {
  // sidebar TOC style
  padding-top: 60px !important;
}
```

#### ローカルCSS

フロントマターでidやクラスを設定できる。@importでcssやlessを読み込める。

```
---
id: "my-id"
class: "my-class1 my-class2"
---

@import "my-style.less"
```

my-style.lessの例

```
#my-id {
  background-color: #222;
  color: #fff;

  h1,
  h2,
  h3,
  h4,
  h5,
  h6 {
    color: #fff;
  }
}
```