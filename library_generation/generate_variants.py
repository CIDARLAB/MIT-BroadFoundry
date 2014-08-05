import pymongo
from pymongo import MongoClient
from bson.objectid import ObjectId
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
        
design = open("../act/test_construction.tab")
labels = design.readline().strip().split("\t")
for line in design:
    line = line.strip().split("\t")
    for path in spec_parse(line[4]):
        db.variants.insert({'name': line[0], 'substrate':line[1], 'product':line[2], 'type':line[3], 'spec':line[4], 'parts': path})