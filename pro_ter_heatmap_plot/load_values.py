from openpyxl import *
from openpyxl.cell import *

def maine():

    directory   = '\Users\Eric\Dropbox (MIT)\Data\Flow'   
    wb_filename = 'P-T Combo All Values'
    sheetname   = 'Log Vals'
    
    load_value_array(directory, wb_filename, sheetname)
    
def load_value_array(directory, wb_filename, sheetname):

    wb = load_workbook('%(dc)s\%(f)s.xlsx'
                           % {'dc':directory,
                              'f':wb_filename},
                       data_only=True)
    
    ws = wb.get_sheet_by_name(sheetname)

    cell_array = ws.rows

    z_values = []
    x_labels = []
    y_labels = []

    value_array = []

    for row in cell_array :
        gfp_values = []

        for cell in row :

            xy = coordinate_from_string(cell.coordinate)
            col = column_index_from_string(xy[0])
            row = xy[1]


            if row is 1 :
                if col >= 2 :

                    x_labels.append(cell.value)

            if row >= 2 :
                if col is 1 :
                    y_labels.append(cell.value)

            if row > 2:              
                if col > 2:
                    if 0 <= cell.value <= 10:
                        z_values.append(cell.value)
                        gfp_values.append(cell.value)
                    else :
                        z_values.append(0)
                        gfp_values.append(0)

        if row > 2:
                value_array.append(gfp_values)

    return [value_array, x_labels, y_labels]
