# Biom_parser
– это простая и быстрая программа обработки данных в файлах формата .biom. Реализация включает в себя использование стандартных модулей многопоточности Python: threading, queue. Использует стандартную API для языка Python для обработки файлов формата .biom под названием «biom».


Формат файла BIOM предназначен как формат общего использования для представления данных биологического эксперимента с помощью таблиц особенностей (метаданных) наблюдения. Формат биома предназначен для общего использования в широких областях сравнительной -омики. Например, в обследованиях маркера генов основным использованием этого формата является представление таблиц OTU (Operational Taxonomic Unit): наблюдения в этом случае являются OTU, а матрица содержит количество, соответствующие количеству случаев каждый OTU, наблюдается в каждом образце.

# Biom_parser
– это простая и быстрая программа обработки данных, основанная на стандартных инструментах многопоточности языка Python. Она принимает на вход таксональную единицу в формате файла .biom (например: sk__, k__, s__), и, обходя все файлы в указанной директории, выводит итоговую статистику, подсчитывающую среднее число малых субъединиц рРНК (SSU rRNA). Реализация включает в себя использование стандартных модулей многопоточности Python: threading, queue. Также стандартную API (Applied Program Interface) для языка Python для обработки файлов формата .biom под названием «biom». Особенность имплементации позволяет быстро получать необходимый результат и легко разбираться в исходном коде парсера. Кроме того, многопоточность позволяет очень быстро обрабатывать входные данные, а варьирование аргумента встроенной функции sleep позволяет обрабатывать большие (расположенные близко к корню иерархического дерева таксонов) таксональные единицы.

Эти данные могут использоваться для дальнейших исследований, а сама программа обработки данных позволяет, не только автоматизировать, но и в значительной степени ускорить процесс получения данных/статистик. К тому же формат .biom распространен в научной среде, поэтому такая программа полезна не только для международного научного сотрудничества, но и просто для анализа данных любых исследований.
