# GTDB-Tk

[GTDB-Tk](https://ecogenomics.github.io/GTDBTk/index.html)
[Ecogenomics/GTDBTk](https://github.com/Ecogenomics/GTDBTk)

```
>conda create -n gtdb-tk
>conda install -c bioconda gtdbtk=2.1.0=pyhdfd78af_5 -y

    GTDB-Tk v2.1.0 requires ~63G of external data which needs to be downloaded
    and extracted. This can be done automatically, or manually.

    Automatic:

        1. Run the command "download-db.sh" to automatically download and extract to:
            /user/ngsdata/kmatsumoto/miniconda3/envs/gtdb-tk/share/gtdbtk-2.1.0/db/

    Manual:

        1. Manually download the latest reference data:
            wget https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/gtdbtk_r207_v2_data.tar.gz

        2. Extract the archive to a target directory:
            tar -xvzf gtdbtk_r207_v2_data.tar.gz -c "/path/to/target/db" --strip 1 > /dev/null
            rm gtdbtk_r207_v2_data.tar.gz

        3. Set the GTDBTK_DATA_PATH environment variable by running:
            conda env config vars set GTDBTK_DATA_PATH="/path/to/target/db"
```

download-db.shでエラー。

```
(gtdb-tk) bash-4.2$ download-db.sh
[INFO] - Downloading the GTDB-Tk database to: /user/ngsdata/kmatsumoto/miniconda3/envs/gtdb-tk/share/gtdbtk-2.1.0/db
--2022-07-13 16:18:10--  https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/gtdbtk_r207_v2_data.tar.gz
Resolving data.gtdb.ecogenomic.org (data.gtdb.ecogenomic.org)... 203.101.230.53
Connecting to data.gtdb.ecogenomic.org (data.gtdb.ecogenomic.org)|203.101.230.53|:443... connected.
ERROR: cannot verify data.gtdb.ecogenomic.org's certificate, issued by ‘/C=US/O=Let's Encrypt/CN=R3’:
  Issued certificate has expired.
To connect to data.gtdb.ecogenomic.org insecurely, use `--no-check-certificate'.
```


/user/ngsdata/kmatsumoto/miniconda3/envs/gtdb-tk/bin/download-db.shの中身のwgetの部分(L28)に`--no-check-certificate`を追加して実行。
