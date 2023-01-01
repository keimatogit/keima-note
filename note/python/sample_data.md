# seabornサンプルデータの使い方

[学習用データセット – seaborn【Python】](https://biotech-lab.org/articles/1408)
[seabornのデータセット（GitHub）](https://github.com/mwaskom/seaborn-data)

途中で[SSL: CERTIFICATE_VERIFY_FAILED]と出た場合は、以下を実行してから再度行います（セキュリティ的にはよくないらしいので気をつけましょう）。

```python
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
```

----

取得できるデータセット名の一覧

```python
import seaborn as sns
print(sns.get_dataset_names())
```

こんな感じで表示されます
```
['anagrams', 'anscombe', 'attention', 'brain_networks', 'car_crashes', 'diamonds', 'dots', 'exercise', 'flights', 'fmri', 'gammas', 'geyser', 'iris', 'mpg', 'penguins', 'planets', 'tips', 'titanic']
```

取得したいデータセット名を指定し、pandas.DataFrameとして読み込みます。

```python
iris = sns.load_dataset('iris')
```
seabornを呼び出さずに、pandasでurlを指定して取得することもできます。
```python
import pandas as pd
iris = pd.read_csv('https://raw.githubusercontent.com/mwaskom/seaborn-data/master/iris.csv')
```
