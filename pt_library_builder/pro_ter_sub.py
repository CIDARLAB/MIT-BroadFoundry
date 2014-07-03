from random import random
from xlrd import *
from excel_functions import cell
from common_functions import clear

def maine() :
    pro_ter_sub(int(raw_input("Design Number: ")))

def pro_ter_sub(design_number) :

    # Define location of data
    iB = 'Library Design.xlsx'
    lS = 'Levels'
    dS = 'Design %(n)s' % {'n' : design_number}

    # Define number of pathways and genes
    np = 192
    ng = 6

    # Clear output file
    clear('Design_%(n)s.txt' % {'n': design_number})

    # Define columns for each gene
    ng = ['B', 'C', 'D', 'E', 'F', 'G']
    

    # Define iterating list for picking PT combinations
    pt_l = ['A', 'B', 'C', 'D', 'E', 'F']
    
    # Populate a dictionary with all promoter terminator combinations
    # and associated expression from GFP cytometry 
    pt_D = {'L': {'A': [cell(iB,lS,'C3'),cell(iB,lS,'B3')],
                  'B': [cell(iB,lS,'C4'),cell(iB,lS,'B4')],
                  'C': [cell(iB,lS,'C5'),cell(iB,lS,'B5')],
                  'D': [cell(iB,lS,'C6'),cell(iB,lS,'B6')],
                  'E': [cell(iB,lS,'C7'),cell(iB,lS,'B7')],
                  'F': [cell(iB,lS,'C8'),cell(iB,lS,'B8')]},
            'M': {'A': [cell(iB,lS,'G3'),cell(iB,lS,'F3')],
                  'B': [cell(iB,lS,'G4'),cell(iB,lS,'F4')],
                  'C': [cell(iB,lS,'G5'),cell(iB,lS,'F5')],
                  'D': [cell(iB,lS,'G6'),cell(iB,lS,'F6')],
                  'E': [cell(iB,lS,'G7'),cell(iB,lS,'F7')],
                  'F': [cell(iB,lS,'G8'),cell(iB,lS,'F8')]},
            'H': {'A': [cell(iB,lS,'K3'),cell(iB,lS,'J3')],
                  'B': [cell(iB,lS,'K4'),cell(iB,lS,'J4')],
                  'C': [cell(iB,lS,'K5'),cell(iB,lS,'J5')],
                  'D': [cell(iB,lS,'K6'),cell(iB,lS,'J6')],
                  'E': [cell(iB,lS,'K7'),cell(iB,lS,'J7')],
                  'F': [cell(iB,lS,'K8'),cell(iB,lS,'J8')]},
            'V': {'A': [cell(iB,lS,'O3'),cell(iB,lS,'N3')],
                  'B': [cell(iB,lS,'O4'),cell(iB,lS,'N4')],
                  'C': [cell(iB,lS,'O5'),cell(iB,lS,'N5')],
                  'D': [cell(iB,lS,'O6'),cell(iB,lS,'N6')],
                  'E': [cell(iB,lS,'O7'),cell(iB,lS,'N7')],
                  'F': [cell(iB,lS,'O8'),cell(iB,lS,'N8')]}}

    # Define counter over pathways
    pathway_number = 0

    while pathway_number < np :

        # Define list to store already used combinations
        already_used = []

        # Define counter over gene number
        gene_number = 0
        
        while gene_number < 6 :

            # Pull desired strength from JMP output in design worksheet
            strength = cell(iB,dS,'%(gn)s%(pn)s'
                            % {'gn': ng[gene_number], 'pn': pathway_number+2})

            if 0 < strength  <= 1000 :
                level = 'L'
                
            if 1000 < strength <= 1500 :
                level = 'M'
                
            if 1500 < strength <= 5000 :
                level = 'H'
                
            if 5000 < strength <=10000 :
                level = 'V'
                

            # Generate random number to choose a combination
            choose = random()*6

            if 0 < choose <= 1 :
                opt = 'A'

            if 1 < choose <= 2 :
                opt = 'B'

            if 2 < choose <= 3 :
                opt = 'C'

            if 3 < choose <= 4 :
                opt = 'D'

            if 4 < choose <= 5 :
                opt = 'E'

            if 5 < choose <= 6 :
                opt = 'F'

            if level+opt not in already_used :

                combination = pt_D[level][opt]

                with open('Design_%(n)s.txt' % {'n': design_number}, 'a') as output_file :

                    output_file.write("%(p)s %(g)s %(c)s %(s)s\n"
                                      % {'p': pathway_number+1, 'g': gene_number+1,
                                         'c': combination[0], 's': int(combination[1])}) 
                
                already_used.append(level+opt)

                gene_number = gene_number + 1

        pathway_number = pathway_number + 1

            
if __name__ == "__maine__" :
    maine()
