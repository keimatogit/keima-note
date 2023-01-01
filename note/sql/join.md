# テーブルの結合(JOIN)

[JOIN 構文](https://dev.mysql.com/doc/refman/5.6/ja/join.html)

USAGE:
```sql
SELECT column1, column2, ...
  FROM table1
  LEFT JOIN table2 -- INNER JOIN, RIGHT JOINも可
  ON conditional_expr
  ```
テーブル間でフィールド名に重複がある場合は、`table1.column1` のようにテーブル名を明示して指定する。



## サンプルデータでの例

MySQLのサンプルデータ「world」を使用します。データベース「world」には、「city」「country」「countrylanguage」の３つのテーブルがあります。

### サンプルデータ「world」

city
```
MariaDB [world]> select * from  city limit 5;
+----+----------------+-------------+---------------+------------+
| ID | Name           | CountryCode | District      | Population |
+----+----------------+-------------+---------------+------------+
|  1 | Kabul          | AFG         | Kabol         |    1780000 |
|  2 | Qandahar       | AFG         | Qandahar      |     237500 |
|  3 | Herat          | AFG         | Herat         |     186800 |
|  4 | Mazar-e-Sharif | AFG         | Balkh         |     127800 |
|  5 | Amsterdam      | NLD         | Noord-Holland |     731200 |
+----+----------------+-------------+---------------+------------+
```
country
```
MariaDB [world]> select Code, Name, Continent, Region, SurfaceArea limit 5;
+------+-------------+---------------+---------------------------+-------------+
| Code | Name        | Continent     | Region                    | SurfaceArea |
+------+-------------+---------------+---------------------------+-------------+
| ABW  | Aruba       | North America | Caribbean                 |      193.00 |
| AFG  | Afghanistan | Asia          | Southern and Central Asia |   652090.00 |
| AGO  | Angola      | Africa        | Central Africa            |  1246700.00 |
| AIA  | Anguilla    | North America | Caribbean                 |       96.00 |
| ALB  | Albania     | Europe        | Southern Europe           |    28748.00 |
+------+-------------+---------------+---------------------------+-------------+
```

countrylanguage
```
MariaDB [world]> select * from  countrylanguage limit 5;
+-------------+------------+------------+------------+
| CountryCode | Language   | IsOfficial | Percentage |
+-------------+------------+------------+------------+
| ABW         | Dutch      | T          |        5.3 |
| ABW         | English    | F          |        9.5 |
| ABW         | Papiamento | F          |       76.7 |
| ABW         | Spanish    | F          |        7.4 |
| AFG         | Balochi    | F          |        0.9 |
+-------------+------------+------------+------------+
```

###

例１）「country」に「countrylanguage」を結合
```
SELECT * FROM country
  LEFT JOIN countrylanguage
  ON countrylanguage.CountryCode = country.Code;
  ```

例２）「country」に「city」「countrylanguage」を結合
```
SELECT * FROM country
  LEFT JOIN (city, countrylanguage)
  ON (city.CountryCODE=country.Code AND countrylanguage.CountryCODE=country.Code);
  ```

例３）「country」に「countrylanguage」を結合し、Englishのpercentageが50以上の国を抽出してpercet降順にsort
  ```
  SELECT name, countrylanguage.Language, countrylanguage.Percentage
    FROM country
    LEFT JOIN countrylanguage
    ON countrylanguage.CountryCode = country.Code
    WHERE countrylanguage.Language = "English" AND countrylanguage.Percentage >= 50
    ORDER BY countrylanguage.Percentage DESC;

    +----------------------+----------+------------+
    | name                 | Language | Percentage |
    +----------------------+----------+------------+
    | Bermuda              | English  |      100.0 |
    | Ireland              | English  |       98.4 |
    | United Kingdom       | English  |       97.3 |
    | Trinidad and Tobago  | English  |       93.5 |
    | Gibraltar            | English  |       88.9 |
    | New Zealand          | English  |       87.0 |
    | United States        | English  |       86.2 |
    | Virgin Islands, U.S. | English  |       81.7 |
    | Australia            | English  |       81.2 |
    | Canada               | English  |       60.4 |
    | Belize               | English  |       50.8 |
    +----------------------+----------+------------+
  ```

例４）「city」に「country」を結合し、Continentごとのcity.Populationの平均などを計算して降順にsort
```
SELECT country.Continent, avg(city.Population), max(city.Population), min(city.Population)
  FROM city
  LEFT JOIN country
  ON country.Code = city.CountryCode
  GROUP BY country.Continent
  ORDER BY avg(city.Population) DESC;

  +---------------+----------------------+----------------------+----------------------+
  | Continent     | avg(city.Population) | max(city.Population) | min(city.Population) |
  +---------------+----------------------+----------------------+----------------------+
  | Asia          |          395019.3109 |             10500000 |                21484 |
  | Africa        |          371143.6585 |              6789479 |                 1500 |
  | South America |          366037.9979 |              9968485 |                 1636 |
  | North America |          289587.5749 |              8591309 |                  595 |
  | Europe        |          287684.6766 |              8389200 |                  455 |
  | Oceania       |          252475.4364 |              3276207 |                   42 |
  +---------------+----------------------+----------------------+----------------------+
  ```
