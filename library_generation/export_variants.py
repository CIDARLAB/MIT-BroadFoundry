import pymongo
from pymongo import MongoClient

client = MongoClient('localhost', 27017)
db = client.act

import argparse

import sys

def parse_args(args):
	
	
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
			description=__doc__)

	parser.add_argument("--output", "-o", dest="output", type=argparse.FileType('w'),
			help="Output file to write data to.")
	parser.add_argument("--atts", "-a", dest="atts", default="",
		help="Attributes to select documents. att:value, comma separated")
	parser.add_argument("--data", dest="data", default="parts",
		help="Attributes to return for selected documents (comma separated). Default is parts (ObjectId is always returned as first column")
	parser.add_argument("--delimeter", "-d", dest="delimeter", default = "\t", help="Text\
			delimeter in output file that separates columns.")
	
	args = parser.parse_args(args)
	
	return args
	
def write_tab(file, labels, data, delimeter):
	file.write("ObjectID{}{}\n".format(delimeter, delimeter.join(labels)))
	for row in data:
		file.write("{}{}{}\n".format(row["_id"], delimeter, delimeter.join([row.setdefault(att, "") for att in labels])))
	file.close()


def main(args):
	options = parse_args(args)
	
	atts = dict([att.split(":") for att in atts.split(",")])
	
	data = []
	for doc in db.variants.find(atts):
    	data.append(doc)
    write_tab(options.output, options.data.split(","), data, options.delimeter)
    
	
if __name__ == "__main__":
	try:
		sys.exit(main(sys.argv[1:]))
	except EnvironmentError as (errno,strerr):
		sys.stderr.write("ERROR: " + strerr + "\n")
		sys.exit(errno)