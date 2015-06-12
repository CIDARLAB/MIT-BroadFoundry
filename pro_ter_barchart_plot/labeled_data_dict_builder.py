from openpyxl import *
from openpyxl.cell import *

def maine():

    directory = '''\Users\Eric\Dropbox (MIT)\Data\Flow'''
    
    filename = 'P-T Combo All Values'

    work_sheet = 'Log Sorted Vals'

    x_names_row = '1'
    y_names_col = 'A'
    first_data_cell = ['C', '3']

    labeled_data_dict_builder(directory,
                              filename,
                              work_sheet,
                              y_names_col,
                              x_names_row,
                              first_data_cell)


def labeled_data_dict_builder(directory,
                              filename,
                              work_sheet,
                              y_names_col,
                              x_names_row,
                              first_data_cell):
    
    workbookname = '%(d)s\%(file)s.xlsx' % {'d': directory,
                                           'file': filename}
    
    input_file = load_workbook(workbookname,data_only=True)

    sheet = input_file.get_sheet_by_name(work_sheet)
    
    cell_array = sheet.rows

    dataD = {}

    namesD = {'1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6,
              'A': 1, 'B': 2, 'C': 3, 'D': 4, 'E': 5, 'F': 6}

    first_data_col = first_data_cell[0]
    first_data_row = first_data_cell[1]

    k = namesD[first_data_col]
    m = namesD[first_data_row]

    x_labelL = []
    y_labelL = []

    # Make two lists of the x and the y labels for dictionary keys
    for row in cell_array :
        for cell in row :
            xy = coordinate_from_string(cell.coordinate)
            col = column_index_from_string(xy[0])
            row = xy[1]

            if row is namesD[x_names_row] :
                if col >= k :
                    x_labelL.append(cell.value)

            if col is namesD[y_names_col] :
                if row >= m :
                    y_labelL.append(cell.value)

    dataD['xlabels'] = x_labelL
    dataD['ylabels'] = y_labelL

    # Construct data using the labels as keys
    for row in cell_array :
        for cell in row :
            xy = coordinate_from_string(cell.coordinate)
            col = column_index_from_string(xy[0])
            row = xy[1]

            if row >= namesD[first_data_row] :

                y = y_labelL[row-m]

                if y not in dataD :

                    dataD[y] = {}

                if col >= namesD[first_data_col] :

                    x = x_labelL[col-k]

                    dataD[y][x] = cell.value

    return dataD
 
    
if __name__ == "__maine__" :
    maine()
