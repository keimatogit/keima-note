### サブクエリの利用〜SELECT結果を条件につかう
```sql
SELECT * FROM table_name
where (column1, column2) in (
		SELECT column1, column2
		from table_name
		group by group_column
		having count(*) > 1000
	)
```

### サブクエリの利用〜SELECT結果からSELECT
```sql
SELECT column1, count(*) FROM
		SELECT column1, group_column, count(*)
		from table_name
		group by group_column
		having count(*) > 1000
	) AS count_table
	GROUP BY column_1
```
