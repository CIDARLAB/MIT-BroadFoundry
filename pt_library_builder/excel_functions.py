import re
from itertools import product
from xlrd import *

def maine() :
    cell(workbook, worksheet, cellname)

def cell(workbook, worksheet, cellname) :

    book = open_workbook(workbook)
    sheet = book.sheet_by_name(worksheet)

    row = int(re.sub("\D", "", cellname))-1
    col = re.sub("\d", "", cellname)
    
    col_list = generate_colnames()

    ncol = col_list.index(col)

    value = sheet.cell(row, ncol).value

    return value

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
