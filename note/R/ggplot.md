# ggplot


テーマのベース

```
# 白背景 + 枠線
theme(text = element_text(size = 18),
	legend.title = element_blank(),
	panel.background = element_rect(fill = "white", colour = "black"))

# 白背景 + 軸線
theme(text = element_text(size = 18),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = NULL),
         axis.line = element_line(lineend = "square"))
```

軸の余白部分

```
scale_x_continuous(expand = c(expand, expand))
```