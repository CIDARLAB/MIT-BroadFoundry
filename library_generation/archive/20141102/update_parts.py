import pymongo
from pymongo import MongoClient
from bson.objectid import ObjectId

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
	
def update(table):
	"""Function takes a list of lists. This is a matrix of data (a table). The first row
	must be column labels. The first column must be ObjectID for the documents to be
	edited/updated.
	"""
	
	labels = table[0]
	
	for row in table[1:]:
   		db.parts.update({labels[0]: ObjectId(row[0])}, {"$set": dict(zip(labels[1:], row[1:]))})
	


def main(args):
	options = parse_args(args)
		
	table = [[ele for ele in row.split(options.delimeter)] for row in options.input.read()]
	update(table)

	
if __name__ == "__main__":
	try:
		sys.exit(main(sys.argv[1:]))
	except EnvironmentError as (errno,strerr):
		sys.stderr.write("ERROR: " + strerr + "\n")
		sys.exit(errno)