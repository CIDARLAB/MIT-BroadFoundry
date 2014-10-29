
import sys
sys.path.append('../')
import gene_cluster_library as gct

a = gct.GeneClusterLibrary()
a.load('test_find_functions.txt')

print 'V1 extract_seq_range:', a.extract_seq_range('V1', 3, 2, 3)
print '-------------------'
print 'V2: extract_seq_range', a.extract_seq_range('V2', 3, 2, 3)
print '-------------------'

print 'V1: find_seq_idx_range', a.find_seq_idx_range('V1', 3)
print '-------------------'
print 'V2: find_seq_idx_range', a.find_seq_idx_range('V2', 3)
print '-------------------'

print 'V1: find_next_part_idx (part_type=None)', a.find_next_part_idx('V1', 1, part_type=None)
print '-------------------'
print 'V2: find_next_part_idx (part_type=None)', a.find_next_part_idx('V2', 1, part_type=None)
print '-------------------'

print 'V1: find_next_part_idx (part_type=CDS)', a.find_next_part_idx('V1', 1, part_type='CDS')
print '-------------------'
print 'V2: find_next_part_idx (part_type=CDS)', a.find_next_part_idx('V2', 4, part_type='CDS')
print '-------------------'

print 'V1 extract_seq_range_around_part:', a.extract_seq_range_around_part('V1', 3, 0, 0)
print '-------------------'
print 'V2: extract_seq_range_around_part', a.extract_seq_range_around_part('V2', 3, 0, 0)
print '-------------------'

nifs = gct.GeneClusterLibrary()
nifs.load('../nif_cluster_data/nif_stata_library.txt')

my_col_mapping = {'P1':2, 'P2':4, 'P3':6,
                  'Rm1':1, 'Rm2':2, 'Rs1':3, 'Rs2':4, 'Ru1':5, 'Ru2':6, 'Rv1':7, 
                  'Rv2':8, 'Rw1':9, 'Rw2':10, 'Rz1':11, 'Rz2':12,
                  'T1':14,
                  'nifM':2, 'nifS':4, 'nifU':6, 'nifV':8, 'nifW':10, 'nifZ':12}
nifs.save_PigeonCAD('./pigeon_output/pigeon_nif_', col_mapping=my_col_mapping)

