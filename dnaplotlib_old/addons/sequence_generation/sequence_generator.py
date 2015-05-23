from plot_SBOL_designs import load_part_information as load_parts
from plot_SBOL_designs import load_dna_designs as load_designs

def maine():
    library_variant_sequence_generator()

def library_variant_sequence_generator():

    parts = load_parts('part_information_yeast.csv')

    designs = load_designs('dna_designs_dsm_yeast.csv', parts)

    design_list = sorted(designs.keys())

    for x, n in enumerate(design_list):
        seq = ''

        for y, q in enumerate(designs[design_list[x]]):
            next_seq = designs[design_list[x]][y]['opts']['sequence']
            seq = seq+next_seq

        with open('%(design)s.txt' % {'design': design_list[x]}, 'w') as f:
            f.write(seq)
    
if __name__ == "__maine__" :
    maine()
