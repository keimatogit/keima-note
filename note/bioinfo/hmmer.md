# HMMER

```
conda create -n hmmer -c bioconda hmmer
```

## hmmer to dbCAN

[dbCAN](https://bcb.unl.edu/dbCAN2/)

dbCAN HMMdb v11.0
[README](https://bcb.unl.edu/dbCAN2/download/Databases/V11/readme.txt)

> ** if you want to run dbCAN CAZyme annotation on your local linux computer, do the following:
** 1. download dbCAN-fam-HMMs.txt, hmmscan-parser.sh 
** 2. download HMMER 3.0 package [hmmer.org] and install it properly
** 3. format HMM db: hmmpress dbCAN-fam-HMMs.txt
** 4. run: hmmscan --domtblout yourfile.out.dm dbCAN-fam-HMMs.txt yourfile > yourfile.out
** 5. run: sh hmmscan-parser.sh yourfile.out.dm > yourfile.out.dm.ps (if alignment > 80aa, use E-value < 1e-5, otherwise use E-value < 1e-3; covered fraction of HMM > 0.3)
** 6. run: cat yourfile.out.dm.ps | awk '$5<1e-15&&$10>0.35' > yourfile.out.dm.ps.stringent (this allows you to get the same result as what is produced in our dbCAN2 webpage)
Cols in yourfile.out.dm.ps:
1. Family HMM
2. HMM length
3. Query ID
4. Query length
5. E-value (how similar to the family HMM)
6. HMM start
7. HMM end
8. Query start
9. Query end
10. Coverage
** About what E-value and Coverage cutoff thresholds you should use (in order to further parse yourfile.out.dm.ps file), we have done some evaluation analyses using arabidopsis, rice, Aspergillus nidulans FGSC A4, Saccharomyces cerevisiae S288c and Escherichia coli K-12 MG1655, Clostridium thermocellum ATCC 27405 and Anaerocellum thermophilum DSM 6725. Our suggestion is that for plants, use E-value < 1e-23 and coverage > 0.2; for bacteria, use E-value < 1e-18 and coverage > 0.35; and for fungi, use E-value < 1e-17 and coverage > 0.45.
** We have also performed evaluation for the five CAZyme classes separately, which suggests that the best threshold varies for different CAZyme classes (please see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4132414/ for details). Basically to annotate GH proteins, one should use a very relax coverage cutoff or the sensitivity will be low (Supplementary Tables S4 and S9); (ii) to annotate CE families a very stringent E-value cutoff and coverage cutoff should be used; otherwise the precision will be very low due to a very high false positive rate (Supplementary Tables S5 and S10)
** On our dbCAN2 website, we use E-value < 1e-15 and coverage > 0.35, which is more stringent than the default ones in hmmscan-parser.sh


```
cd ~/db/dbCAN
sbatch -J wget -p c6420 -c 22 --wrap="wget https://bcb.unl.edu/dbCAN2/download/dbCAN-HMMdb-V11.txt"
wget https://bcb.unl.edu/dbCAN2/download/Databases/V11/hmmscan-parser.sh
```

hmmscan用レファレンス作成？

```
srun -c 22 -p c6420 hmmpress dbCAN-HMMdb-V11.txt

# ↓4ファイルできる
# Pressed and indexed 699 HMMs (699 names and 9 accessions).
# Models pressed into binary file:   dbCAN-HMMdb-V11.txt.h3m
# SSI index for binary model file:   dbCAN-HMMdb-V11.txt.h3i
# Profiles (MSV part) pressed into:  dbCAN-HMMdb-V11.txt.h3f
# Profiles (remainder) pressed into: dbCAN-HMMdb-V11.txt.h3p
```

hmmscan

`hmmscan [-options] <hmmdb> <seqfile>`
`--domtblout`: save parseable table of per-domain hits to file
Evalueのデフォは10

```
hmmscan \
  -E 1e-10 \
  --cpu 22 \
  --domtblout output.out.dm \
  -o output.txt \
  ~/db/dbCAN/dbCAN-fam-HMMs.txt \
  input.fasta
```

hmmscan-parser.sh
> if alignment > 80aa, use E-value < 1e-5, otherwise use E-value < 1e-3; covered fraction of HMM > 0.3

```
sh hmmscan-parser.sh \
  output.out.dm > output.out.dm.ps 
```

dbCAN2 webpageと同じ結果にするには（E-value < 1e-15 and coverage > 0.35）
上のhmmscan-parser.shの設定はwebpageより甘いので。

```
cat output.out.dm.ps | awk '$5<1e-15&&$10>0.35' > output.out.dm.ps.stringent
```