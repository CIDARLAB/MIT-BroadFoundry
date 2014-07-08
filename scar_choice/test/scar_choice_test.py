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

print sc.find_scars(4, seed_set=[], max_homology=3, num_to_find=3, random_search=True)
print '---------'
print sc.find_scars(4, seed_set=['AATA'], max_homology=3, num_to_find=3, random_search=True)
print '---------'
print sc.find_scars(4, seed_set=['AATA'], max_homology=3, num_to_find=3, random_search=True)
print '---------'
print sc.find_scars(4, seed_set=['AATA'], max_homology=3, num_to_find=3, random_search=False)
print '---------'
print sc.find_scars(4, seed_set=['AATA'], max_homology=3, num_to_find=None, random_search=False)
