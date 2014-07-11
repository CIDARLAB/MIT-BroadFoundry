#!/usr/bin/env python
"""
Scar Choice Tests
=================

    Script to test some of the scar choice functions.
"""
#    Scar Choice
#    Copyright (C) 2014 by
#    Thomas E. Gorochowski <tom@chofski.co.uk>
#    All rights reserved.
#    OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

import scar_choice as sc

print sc.find_scars(4, seed_set=[], max_homology=2, num_to_find=3, random_search=True)
print '---------'
print sc.find_scars(4, seed_set=['AATA'], max_homology=2, num_to_find=3, random_search=True)
print '---------'
print sc.find_scars(4, seed_set=['AATA'], max_homology=2, num_to_find=3, random_search=True)
print '---------'
print sc.find_scars(4, seed_set=['AATA', 'AGTA'], max_homology=2, num_to_find=3, random_search=False)
print '---------'
full_scar_set, new_scars, removed_from_seed = sc.find_scars(4, seed_set=[], max_homology=2, num_to_find=None, random_search=False, allowed_set=[])
print 'Should be 32 scars in total:', len(full_scar_set)
print 'Full set =', len(full_scar_set), full_scar_set
print 'Added set =', len(new_scars), new_scars
print 'Removed =', len(removed_from_seed), removed_from_seed
