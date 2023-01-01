# 例外処理（try, tryCatch）


途中でエラーが発生してもプログラムを止めずに続けたり、エラーや警告に対応できる。

try
tryCatch


oCb <- try(CellBlock(oWs, intXlsxTop + (!is.na(strCompVarLabel2)), intXlsxLeftTemp, intXlsxTopTemp - intXlsxTop - (!is.na(strCompVarLabel2)) + blnOutWhole, intvcColVarN[i] * intFormatSize, FALSE), silent = TRUE)
if (class(oCb) == "try-error"){
	oCb <- CellBlock(oWs, intXlsxTop + (!is.na(strCompVarLabel2)), intXlsxLeftTemp, intXlsxTopTemp - intXlsxTop - (!is.na(strCompVarLabel2)) + blnOutWhole, intvcColVarN[i] * intFormatSize, FALSE)
}

tryCatch({data <- filter(data_original, eval(parse(text = args_target)))},
         error = function(e) {
           cat("データ抽出に失敗しました。'args_target'を見直してください：\n\n")
           cat(args_target, "\n\n")
           cat("参考サイト https://dplyr.tidyverse.org/reference/filter.html \n\n")
           stop(e)
           },
         message = function(w){
           cat("'args_target'でのデータ抽出時に警告がありました：/n", w)
           }
         )
