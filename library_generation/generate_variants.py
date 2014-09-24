import pymongo
from pymongo import MongoClient
from bson.objectid import ObjectId
import sys
import argparse
client = MongoClient('localhost', 27017)
db = client.act

def combine(s1, s2):
    if len(s1) == 0:
        return s2
    if len(s2) == 0:
        return s1
    combs = []
    for e1 in s1:
        for e2 in s2:
            combs.append(e1 + e2)
    return combs
    
def part_lookup(spec):
    spec = spec.split(";")
    atts = dict([att.split(":") for att in spec])
    return [[ObjectId(part['_id'])] for part in db.parts.find(atts)]
    
def spec_parse(spec):
	"""function accepts a string representing a 'spec'. This included (), *, +, and data
	specifying documents in the database. Example:
		A+B = AB
		A*B = AB, BA
		A+B*C = (A+B)*C = ABC, CAB
		(A*B)+(C*D) = ABCD, BACD, ABDC, BADC
	Each variable (A, B, C, D) can be a part:
		substrate:L-TRYPTOPHAN;product:5-HYDROXY-L-TRYPTOPHAN+substrate:5-HYDROXY-L-TRYPTOPHAN;product:SEROTONIN
	this concatenates all combinations of parts that satisfy that spec.
	"""
		
    if not spec.count("(") == spec.count(")"):
        print "Missing Parentheses"
        return
    
    if spec.endswith(")"):
        record = 1
        rec = []
        opened = 0
        closed = 0
        for ind, char in enumerate(spec[::-1]):
            if record:
                if char == "(":
                    opened += 1
                if char == ")":
                    closed += 1
                rec.append(char)
                if opened == closed:
                    if ind == len(spec) - 1:
                        return spec_parse("".join(rec[1:-1])[::-1])
                    record = 0
                    continue
            else:
                rec = "".join(rec[1:-1])[::-1]
                if char == "*":
                    return combine(spec_parse(spec[:len(spec)-ind-1]), spec_parse(rec)) + combine(spec_parse(rec), spec_parse(spec[:len(spec)-ind-1]))
                if char == "+":
                    return combine(spec_parse(spec[:len(spec)-ind-1]), spec_parse(rec))
                else:
                    return 'uh oh'
        
    plus = spec.rfind("+")
    star = spec.rfind("*")
    if star > plus:
        l, r = spec.rsplit("*", 1)
        return combine(spec_parse(l), spec_parse(r)) + combine(spec_parse(r), spec_parse(l))
    if plus > star:
        l, r = spec.rsplit("+", 1)
        return combine(spec_parse(l), spec_parse(r))
    if star == plus == -1:
        return part_lookup(spec)
        
def parse_args(args):
	
	
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
			description=__doc__)

	parser.add_argument("--input", "-i", dest="input", type=argparse.FileType('r'),
			help="Input file to write data to.")
	parser.add_argument("--delimeter", "-d", dest="delimeter", default = "\t", help="Text\
			delimeter in output file that separates columns.")
	
	args = parser.parse_args(args)
	
	return args
        
def insert_variants(table):
	
	labels = table[0]
	for row in table[1:]:
	    doc = dict(zip(labels,row))
	    for path in spec_parse(doc['spec']):
	        db.variants.insert(dict(zip(labels, row) + ['parts', path]))
	        
def main(args):
	
	options = parse_args(args)
    table = [[ele for ele in row.split(options.delimeter)] for row in options.input.read()]
	insert_variants(table)
	
        
if __name__ == "__main__":
	try:
		sys.exit(main(sys.argv[1:]))
	except EnvironmentError as (errno,strerr):
		sys.stderr.write("ERROR: " + strerr + "\n")
		sys.exit(errno)