import csv
from biom import load_table
from collections import Counter
import threading
import queue
from time  import sleep
import os


# Определим таксон по умолчанию
TARGET = "s__"

# Количество потоков по умолчанию
Q_THREADS = 4

# Инициализируем количество упоминаний заданного таксона и суммарную SSU rRNA и 
# составим словарь-счетчик всех упоминаний различных  таксонов во всех файлах и их SSU rRNA

SSU_rRNA_by_taxon = Counter()
Number_by_taxon = Counter()

# Обьявляем очередь для передачи данных из потоков в Общие счетчики-словари
queue_collect_SSU = queue.Queue()

def find_taxon(taxon,list_str):
    k = False
    for i in list_str:
        if taxon in i:
            k = True
            break
    return k


# Эта функция выгружает данные из файла biom в обьект класса biom. 
# Аргументы:
# 1. Массив полных имен файлов 
# 2. Имя таксона 
# Затем методами класса biom проводит фильтрацию по всем подходящим таксонам 
# создает словарь и выгружает словарь этого файла в очередь
# Функция запускается в отдельных потоках
 
def parse_biom(arr, taxon):

    Dict_SSU_temp = {}
    Name_of_taxon_temp = ''
    SSU_rRNA_by_taxon_temp = 0

    for i in arr:
        try:
            table = load_table( i )
        except TypeError:
            continue

        
# Выбор фильтра в зависимости от типа "taxonomy"

        for i in table.ids(axis='observation'):

            if isinstance(table.metadata(id = i, axis = 'observation').get('taxonomy'),list):
                f_filter = lambda values, id_, md: find_taxon(taxon,md['taxonomy'])

            else:
                f_filter = lambda values, id_, md: taxon in md['taxonomy']

            break

        filtered_table = table.filter(f_filter, axis='observation', inplace=True)

        for i in filtered_table.ids(axis='observation'):

            Name_of_taxon_temp_buff = filtered_table.metadata(id = i, axis = 'observation').get('taxonomy')
            
            if isinstance(Name_of_taxon_temp_buff,list):
                Name_of_taxon_temp = ";".join(Name_of_taxon_temp_buff)
            else:  
                Name_of_taxon_temp = Name_of_taxon_temp_buff.replace(",",";")        
                
            SSU_rRNA_by_taxon_temp = filtered_table.data(id = i, axis='observation')[0]
            Dict_SSU_temp[Name_of_taxon_temp] = SSU_rRNA_by_taxon_temp

  
# Выгрузка происходит в  очередь

        queue_collect_SSU.put(Dict_SSU_temp)

    return

# Эта функция выгружает данные из файла tsv в массив. 
# Аргументы:
# 1. Массив полных имен файлов 
# 2. Имя таксона 
# Затем проводит фильтрацию массива по всем подходящим таксонам 
# создает словарь и выгружает словарь этого файла в очередь 
# Функция запускается в отдельных потоках

def parse_tsv(arr,taxon):

    Dict_SSU_temp = {}
    Name_of_taxon_temp = ''
    SSU_rRNA_by_taxon_temp = 0

    for i in arr:

        with open(i, 'r') as f:
            try:
                tsv_reader = csv.reader(f, delimiter='\t')
            except TypeError:
                continue
            headers = next(tsv_reader)
            if len(headers) < 4:
                 headers = next(tsv_reader)
            try:
                index_taxonomy = headers.index('taxonomy')
                index_SSU_rRNA = headers.index('SSU_rRNA')
            except:
                print(i)
                continue

            for row in tsv_reader:
                if len(row) == 4:
                    if taxon in row[index_taxonomy]:
                        Name_of_taxon_temp = row[index_taxonomy].replace(",",";")
                        SSU_rRNA_by_taxon_temp = float(row[index_SSU_rRNA])
                        Dict_SSU_temp[Name_of_taxon_temp] = SSU_rRNA_by_taxon_temp

            queue_collect_SSU.put(Dict_SSU_temp)

    return

# Функция Collector  забирает информацию из очереди и добавляет ее в соответствующие счетчики

def collector():

    while True:

        dict_temp = queue_collect_SSU.get()
        SSU_rRNA_by_taxon.update(dict_temp)

        for i in dict_temp : dict_temp[i] = 1

        Number_by_taxon.update(dict_temp)


# Функция array_of_files(str) получает в качестве аргумента расширение имен файлов, а возвращает
#  список имен файлов с этим расширением  из папки ~/data 

#Путь к основной директории
MAIN_DIRECTORY = ""

def array_of_files(ext_str):
    dir_name = MAIN_DIRECTORY + ext_str + "\\"
    #dir_name = os.getcwd() + "/data"
    arr_temp = list(filter(lambda  val: ("."+ext_str) in val, os.listdir(dir_name)))

    for i in range(len(arr_temp)):
        arr_temp[i] = dir_name + "\\" + arr_temp[i]

    return arr_temp

if __name__ == "__main__":

# Ввод имени таксона. Если имя не вводится то используется таксон по умолчанию

    taxon = input("Type the name of the taxon or press Enter (Default Value %s )\n\n" % TARGET)

    if taxon == "" : taxon = TARGET


# Получение списка имен файлов biom из папки

    biom_arr = array_of_files("biom")   

    Number_of_Files = len(biom_arr)


#Разделим исходный массив на части для каждого потока

    L_Q_Threads = Number_of_Files//Q_THREADS
    
    biom_arr_by_threads = [[]*Q_THREADS for i in range(Q_THREADS)]

    for i in range(Q_THREADS-1):
        biom_arr_by_threads[i] = biom_arr[L_Q_Threads*i:L_Q_Threads*(i+1)]

    biom_arr_by_threads[Q_THREADS - 1] = biom_arr[L_Q_Threads*(Q_THREADS-1) : Number_of_Files]


# Формируем массив потоков
    threads = [threading.Thread(target = parse_biom, args = (biom_arr_by_threads[i],taxon,), daemon=True) for i in range(Q_THREADS)]

# Формируем поток для collector

    thread_collector = threading.Thread(target = collector, daemon=True)

# Старт потоков
    thread_collector.start()

    for thread in threads:
        thread.start()
    
    for thread in threads:
        thread.join()


# Получение списка имен файлов tsv из папки

    tsv_arr = array_of_files("tsv")   

    Number_of_Files = len(tsv_arr)


#Разделим исходный массив на части для каждого потока


    L_Q_Threads = Number_of_Files//Q_THREADS
    
    tsv_arr_by_threads = [[]*Q_THREADS for i in range(Q_THREADS)]

    for i in range(Q_THREADS - 1):
        tsv_arr_by_threads[i] = tsv_arr[L_Q_Threads*i:L_Q_Threads*(i+1)]

    tsv_arr_by_threads[Q_THREADS - 1] = tsv_arr[L_Q_Threads*(Q_THREADS-1) : Number_of_Files]


# Формируем массив потоков
    threads = [threading.Thread(target = parse_tsv, args = (tsv_arr_by_threads[i],taxon,), daemon=True) for i in range(Q_THREADS)]

# Старт потоков
    for thread in threads:
        thread.start()
    
    for thread in threads:
        thread.join()

# Вывод результата обработки в файл
    #это необходимо для больших таксональных единиц, например sk__
    sleep(20)
    #Путь к .csv файлу, куда будут записаны данные
    with open(MAIN_DIRECTORY +'\\taxon.csv', 'w', newline='') as csvfile:

        writer = csv.writer(csvfile, quoting=csv.QUOTE_NONNUMERIC)
        writer.writerow(["Taxon", "Average_SSU_rRNA", "Total_number"])
        arr_temp = [[],[],[]]

        for i in SSU_rRNA_by_taxon:
            arr_temp[0] = i
            arr_temp[1] = SSU_rRNA_by_taxon[i]/Number_by_taxon[i]
            arr_temp[2] = Number_by_taxon[i]
            writer.writerow(arr_temp)
