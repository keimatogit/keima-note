# 外部キー制約、check制約

カラムに格納する値に制約をつける。
外部キー制約は、関連したテーブル間をで整合性を保証するもので、子テーブルに設定する。
check制約はTrue/Falseとなる条件式で制約を指定する。
外部キー制約をつけたカラムは自動的にインデックス（キー）が作成される。

テーブルの制約(CONSTRAINT)設定状態は、`SHOW CREATE TABLE tabelname;`で確認できる（`CONSTRAINT〜`のところ）。


<!-- @import "[TOC]" {cmd="toc" depthFrom=2 depthTo=6 orderedList=false} -->
<!-- code_chunk_output -->

- [外部キー制約](#外部キー制約)
  - [テーブル作成時に外部キー制約を設定](#テーブル作成時に外部キー制約を設定)
  - [外部キー制約の追加](#外部キー制約の追加)
  - [外部キー制約の削除](#外部キー制約の削除)
  - [ON DELETE、ON UPDATE（親テーブルでの削除時と更新時）の設定](#on-delete-on-update親テーブルでの削除時と更新時の設定)
  - [外部キー制約の挙動](#外部キー制約の挙動)
- [CHECK制約](#check制約)
  - [テーブル作成時にCHECK制約を設定](#テーブル作成時にcheck制約を設定)
  - [CHECK制約を追加](#check制約を追加)
  - [CHECK制約を削除](#check制約を削除)
  - [各フィールド情報でCHECK制約をつける（非推奨）](#各フィールド情報でcheck制約をつける非推奨)

<!-- /code_chunk_output -->


## 外部キー制約

### テーブル作成時に外部キー制約を設定

```sql
CREATE TABEL tabelname(
  columnname INT, ...,

  # 詳しい設定
  [CONSTRAINT [constraint_name]] -- constraint_nameを明記したい時
    FOREIGN KEY [index_name] (columnname) -- index_nameを明記したい時
      REFERENCE parent_tabelname(parent_column_name)
      [ON DELETE {RESTRICT, CASCADE, SET NULL}] -- 親テーブルでの値の削除時の挙動を指定したい時
      [ON UPDATE {RESTRICT, CASCADE, SET NULL}], -- 親テーブルでの値の更新時の挙動を指定したい時

  # 基本はこれだけで良い
  FOREIGN KEY (columnname2)
    REFERENCE parent_tabelname(parent_column_name2)
)


```

### 外部キー制約の追加
詳しく指定したい場合は上記参照

```sql
ALTER TABLE table_name ADD
  FOREIGN KEY (columnname)
  REFERENCES parent_tabelname(parent_columnname);
```

### 外部キー制約の削除

constraint_nameは、`SHOW CREATE TABLE table_name;`で確認
```sql
ALTER TABLE table_name DROP FOREIGN KEY constraint_name;
```

インデックス(key)は残る。インデックスの削除はこちら。
インデックス名は、`SHOW CREATE TABLE table_name;`で確認（`FOREIGN KEY ()`内）。
```sql
ALTER TABLE table_name DROP INDEX index_name;
```

### ON DELETE、ON UPDATE（親テーブルでの削除時と更新時）の設定
- `RESTRICT`、`CASCADE`、`SET NULL`の３種類。デフォルトは`RESTRICT`。
- `RESTRICT`: 子テーブルにある値を親テーブルで削除や更新できない（エラーが出る）
- `CASCADE`: 子テーブルにある値を親テーブルで削除したら、子テーブルの該当レコードも削除される。親テーブルで値を更新したら、子テーブルの値も更新される。（子テーブルから先に値を更新することはできない）
- `SET NULL`: 子テーブルにある値を親テーブルで削除や更新すると、子テーブルの値はNULLになる。

### 外部キー制約の挙動
- 親テーブルを先に作成していないと、外部キー制約つきテーブルはcreateできない。
- 子テーブルに制約外となる値がある場合、外部キー制約を追加できない（ON DELETE/UPDATE設定によらず）。
- pythonのpandas.dataframe.to_sqlでは、１行でもFK制約を満たさない行があれば全部インポートされなかった。
- sqlのload dataでは、制約を満たさない行は無視され、他の行はインポートされた。外部キー制約つきカラムが空白の行はNULL（FK設定カラムのデフォルト）としてインポートされた。


## CHECK制約

### テーブル作成時にCHECK制約を設定

```sql
CREATE TABEL tabelname(
  columnname INT, ...

  [CONSTRAINT [constraint_name]] -- constraint_nameを明記したい時
    CHECK (columnname > 100)
)
```

### CHECK制約を追加

```sql
ALTER TABLE table_name ADD
  CHECK (columnname > 100);
```

### CHECK制約を削除

constraint_nameは、`SHOW CREATE TABLE table_name;`で確認
```sql
ALTER TABLE table_name DROP CONSTRAINT constraint_name;
```

### 各フィールド情報でCHECK制約をつける（非推奨）

各フィールド情報にCHECK制約をつけることもできる
```sql
CREATE TABEL tabelname(
  columnname INT CHECK (columnname > 100), ...
)
```

なので、各フィールド情報の編集でCHECK制約を追加することもできる
```sql
ALTER TABLE tabelname MODIFY colmnname2 INT CHECK (columnname > 200);
```
このようにALTER TABLE tabelname MODIFYでフィールドにCHECK制約を追加した場合、既存レコードの値は制約外のものもそのまま残り、新たに格納されるレコードにのみCHECK制約が適応される（CONSTRAINTでは制約外のレコードがあると制約の追加自体できない）。

なんだかややこしいので、この方法は使わず、テーブル全体のCONSTRAINTでつけることにしようと思う。

CONSTRAINTでついているかフィールドについているかは、`SHOW CREATE TABLE table_name;`で分かる。
