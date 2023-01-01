# CAZy

[CAZy](http://www.cazy.org/)

hmmer.mdも参照。

```
cd ~/db/cazy/221209

# dbCAN(https://bcb.unl.edu/dbCAN2/)からのダウンロード
sbatch -c 22 -p c6420 -J wget_dbCAN --wrap="wget https://bcb.unl.edu/dbCAN2/download/CAZyDB.08062022.fa"

# CAZyホームページ。これ何?
wget http://www.cazy.org/IMG/cazy_data/cazy_data.zip
```


# run_dbcan

[run_dbcan3](https://github.com/linnabrown/run_dbcan)