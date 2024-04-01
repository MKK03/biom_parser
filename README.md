# Biom_parser
## Requirments
* Python 3.7+ (Tested for Windows 10, Python 3.10)
* [BIOM package](https://biom-format.org/)
## Description
### English:
– It is a simple and fast data processing program based on standard Python multithreading tools. It accepts a taxonal unit as input in a file format.biom (for example: sk__, k__, s__), and, bypassing all files in the specified directory, outputs the final statistics calculating the average number of small rRNA subunits (SSU rRNA). The implementation includes the use of standard Python multithreading modules: threading, queue. There is also a standard API (Applied Program Interface) for the Python language for processing format files.a biom called "biom". The implementation feature allows you to quickly get the desired result and easily understand the source code of the parser. In addition, multithreading allows you to process input data very quickly, and varying the argument of the built-in sleep function allows you to process large (located close to the root of the hierarchical tree of taxa) taxonal units.

The BIOM file format is intended as a general use format for representing biological experiment data using tables of observation features (metadata). The biome format is intended for general use in wide areas of comparative mathematics. For example, in gene marker surveys, the main use of this format is the presentation of OTU (Operational Taxonomic Unit) tables: the observations in this case are OTU, and the matrix contains numbers corresponding to the number of cases each OTU observed in each sample.

This data can be used for further research, and the data processing program itself allows not only to automate, but also to significantly speed up the process of obtaining data/ statistics. Besides , the format .biom is widespread in the scientific community, so such a program is useful not only for international scientific cooperation, but also simply for analyzing data from any research.
### Russian:
– это простая и быстрая программа обработки данных, основанная на стандартных инструментах многопоточности языка Python. Она принимает на вход таксональную единицу в формате файла .biom (например: sk__, k__, s__), и, обходя все файлы в указанной директории, выводит итоговую статистику, подсчитывающую среднее число малых субъединиц рРНК (SSU rRNA). Реализация включает в себя использование стандартных модулей многопоточности Python: threading, queue. Также стандартную API (Applied Program Interface) для языка Python для обработки файлов формата .biom под названием «biom». Особенность имплементации позволяет быстро получать необходимый результат и легко разбираться в исходном коде парсера. Кроме того, многопоточность позволяет очень быстро обрабатывать входные данные, а варьирование аргумента встроенной функции sleep позволяет обрабатывать большие (расположенные близко к корню иерархического дерева таксонов) таксональные единицы.

Формат файла BIOM предназначен как формат общего использования для представления данных биологического эксперимента с помощью таблиц особенностей (метаданных) наблюдения. Формат биома предназначен для общего использования в широких областях сравнительной -омики. Например, в обследованиях маркера генов основным использованием этого формата является представление таблиц OTU (Operational Taxonomic Unit): наблюдения в этом случае являются OTU, а матрица содержит количество, соответствующие количеству случаев каждый OTU, наблюдается в каждом образце.

Эти данные могут использоваться для дальнейших исследований, а сама программа обработки данных позволяет, не только автоматизировать, но и в значительной степени ускорить процесс получения данных/статистик. К тому же формат .biom распространен в научной среде, поэтому такая программа полезна не только для международного научного сотрудничества, но и просто для анализа данных любых исследований.
