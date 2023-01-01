# illumina BaseSpace

わりと面倒そうでまだ使ってみてない

[BaseSpace login](https://basespace.illumina.com/settings/account)
[BaseSpace help](https://help.basespace.illumina.com/)


[BaseSpace - Biosample Workflow Files](https://help.basespace.illumina.com/files-used-by-basespace/biosample-workflow-files)




アップロードするfastqファイルのフォーマットにきびしい。
[Import Data Into Projects](https://help.basespace.illumina.com/manage-data/import-data)

> FASTQ files need to adhere to Illumina standards, as specified below:
- Data for a single sample can constitute multiple files. The total number of files per sample and their combined size are limited to 16 and 25 GB respectively.
The uploader will only support gzipped FASTQ files generated on Illumina instruments.
- The name of the FASTQ files must conform the following convention: SampleName_SampleNumber_Lane_Read_FlowCellIndex.fastq.gz (i.e. SampleName_S1_L001_R1_001.fastq.gz / SampleName_S1_L001_R2_001.fastq.gz)
- The read descriptor in the FASTQ files must conform to the following convention: @Instrument:RunID:FlowCellID:Lane:Tile:X:Y ReadNum:FilterFlag:0:SampleNumber:
- Read 1 descriptor would look like this: @M00900:62:000000000-A2CYG:1:1101:18016:2491 1:N:0:13
- Read 2 would have a 2 in the ReadNum field, like this: @M00900:62:000000000-A2CYG:1:1101:18016:2491 2:N:0:13

- 形式はfastq.gz
- ファイル名は「SampleName_SampleNumber_Lane_Read_FlowCellIndex.fastq.gz」(ex. SampleName_S1_L001_R1_001.fastq.gz, SampleName_S1_L001_R2_001.fastq.gz)
- リード名は「@Instrument:RunID:FlowCellID:Lane:Tile:X:Y ReadNum:FilterFlag:0:SampleNumber:」

OK(ENA:PRJEB38984): @HWI-D00483:156:C9JY7ANXX:5:1101:1497:2187 1:N:0:GAATTCGTCCTATCCT
NG(d178/F3319): @V35_D2B_211109_0030461L1C001R0020000002/1
