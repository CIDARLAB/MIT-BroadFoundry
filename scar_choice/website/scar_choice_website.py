#!/usr/bin/env python
"""
Scar Choice Website
===================

    Flask-based website to enable easy access to the scar choice scripts.
"""
#    Scar Choice Website
#    Copyright (C) 2014 by
#    Thomas E. Gorochowski <tom@chofski.co.uk>
#    All rights reserved.
#    OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'

from flask import Flask
app = Flask(__name__)

@app.route('/')
def sc_main():
    return 'Hello World!'

if __name__ == '__main__':
    app.run()
