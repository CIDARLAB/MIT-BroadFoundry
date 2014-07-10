import re
from itertools import product
from xlrd import *

global cache
cache = {}

def maine() :
    cell(workbook, worksheet, cellname)
    w_cell(workbook, worksheet, cellname, cellvalue)

def cell(workbook, worksheet, cellname) :

    k = workbook+worksheet+cellname

    if k in cache :

        value = cache[k]

    else :
        
        book = open_workbook(workbook)
        
        sheet = book.sheet_by_name(worksheet)

        row = int(re.sub("\D", "", cellname))-1
        col = re.sub("\d", "", cellname)
    
        col_list = generate_colnames()

        ncol = col_list.index(col)

        value = sheet.cell(row, ncol).value

        cache[k] = value

    return value

def w_cell(workbook, worksheet, cellname, cellvalue) :

    book = workbook
    sheet = worksheet

    row = int(re.sub("\D", "", cellname))-1
    col = re.sub("\d", "", cellname)
    
    col_list = generate_colnames()

    ncol = col_list.index(col)

    sheet.write(row, ncol, cellvalue)

def generate_colnames():
    
    col1 = product('ABCDEFGHIJKLMNOPQRSTUVWXYZ',repeat=1)
    col2 = product('ABCDEFGHIJKLMNOPQRSTUVWXYZ',repeat=2)
    cs = []

    for x in col1 :
        c = ''.join(x)
        cs.append(c)
        
    for y in col2 :
        c = ''.join(y)
        cs.append(c)

    return cs

if __name__ == "__maine__" :
    maine()
