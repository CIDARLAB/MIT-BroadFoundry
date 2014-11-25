#DNAplotlib

![DNAplotlib](http://www.chofski.co.uk/img/dnaplotlib/dnaplotlib.jpg)

DNAplotlib is a computational toolkit that enables highly customizable visualization of individual genetic constructs and libraries of design variants. Publication quality vector-based output is produced and all aspects of the rendering process can be easily customized or replaced by the user. DNAplotlib is capable of SBOL Visual compliant diagrams in addition to a format able to better illustrate the precise location and length of each genetic part. This alternative visualization method enables direct comparison with nucleotide-level information such as RNA-seq read depth. While it is envisaged that access will be predominantly via the programming interface, several easy to use text-based input formats can be processed by a command-line scripts to facilitate broader usage. A web front-end is also available.

##Installation
The DNAplotlib library is contained within the `dnaplotlib.py` file and requires Python 2.6 and matplotlib 1.2 or newer. To install add the location of this file to your `PYTHONPATH` and you are good to:

``import dnaplotlib``

##Getting Started
We provide an extensive gallery of use cases for DNAplotlib in the `gallery` directory. These cover accessing the library directly with Python, to using the quick construct creator (`quick.py`) and the CSV file input script (`plot_SBOL_designs.py`). If you would like to set up your own web-based server we provide the entire web front-end in the `web` directory.

##Usage and Licensing
DNAplotlib is cross-platform and open-source software developed using Python and released under the OSI recognized NPOSL-3.0 license.
