import pymongo
from pymongo import MongoClient
client = MongoClient('localhost', 27017)
db = client.act

import argparse

import sys

def parse_args(args):
	
	
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
			description=__doc__)

	parser.add_argument("--input", "-i", dest="input", type=argparse.FileType('r'),
			help="Input file to read data from.")
	parser.add_argument("--delimeter", "-d", dest="delimeter", default = "\t", help="Text\
			delimeter in input file that separates columns.")
	
	args = parser.parse_args(args)
	
	return args

#def insert(


def main(args):
	options = parse_args(args)
	
	actzymes = options.input

	labels = actzymes.readline().strip("\n").split(options.delimeter)

	for line in actzymes:
		if not line:
			continue
		line = line.strip("\n").split(options.delimeter)
		db.parts.insert(dict(zip(labels, line)))	
	
if __name__ == "__main__":
	try:
		sys.exit(main(sys.argv[1:]))
	except EnvironmentError as (errno,strerr):
		sys.stderr.write("ERROR: " + strerr + "\n")
		sys.exit(errno)