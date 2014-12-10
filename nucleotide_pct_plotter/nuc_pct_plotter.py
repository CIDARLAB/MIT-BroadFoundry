from openpyxl import *
from openpyxl.cell import *
import matplotlib.pyplot as plt

def maine():
    nuc_pct_plotter()

def nuc_pct_plotter():
    
    wb = load_workbook('\Users\Eric\Dropbox\MIT-Broad\Data\Parts Analysis\Nucleotide Percentages.xlsx')
    
    ws = wb.get_sheet_by_name('PT Stepcharts')

    array = ws.rows

    cache = {}

    ATCG = ['A', 'T', 'C', 'G']
    
    for row in array :
        
        for cell in row :
            xy = coordinate_from_string(cell.coordinate)
            col = column_index_from_string(xy[0])
            row = xy[1]
    
            if col is 1 :
                k = cell.value
                dD = {}
                if k not in cache :
                    cache[k] = dD
                    cache[k]['x_values'] = []

                    for n in ATCG :
                        cache[k][n] = []
                        
            if col is 2 :
                x = cell.value
                cache[k]['x_values'].append(x)

            if 3 <= col <= 6 :
                y = cell.value
                cache[k][ATCG[col-3]].append(y)

    count = 1
    fig = plt.figure(1)
    
    while count <= 36 :
        plt.subplot(9,4,count)
        
        x_val = cache['>PRO%(num)s'%{'num':count}]['x_values']

        plt.plot(x_val, cache['>PRO%(num)s'%{'num':count}]['A'], 'g',
                 x_val, cache['>PRO%(num)s'%{'num':count}]['T'], 'r',
                 x_val, cache['>PRO%(num)s'%{'num':count}]['C'], 'b',
                 x_val, cache['>PRO%(num)s'%{'num':count}]['G'], 'k')

        plt.title('Promoter %(num)s'%{'num':count}, loc='left')
        plt.ylabel('%')
        plt.ylim([0,100])


        count += 1
    
    fig.legend((l1,l2,l3,l4), ('A','T','C','G'), 'upper left')
        
    plt.show()
    
if __name__ == "__maine__" :
    maine()
