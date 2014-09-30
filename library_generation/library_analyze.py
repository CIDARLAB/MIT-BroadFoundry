
import dbmap
from sqlalchemy.orm import sessionmaker

def parse_line(parts):
	"""Function accepts a line of text, with key=value pairs separated by tabs. Returns
	a dictionary with each key: value pair.
	"""
	line_parts = dict()
	for pair in parts:
		pair = pair.split("=", 1)
		line_parts[pair[0]] = pair[1]
	return line_parts
	
		
def import_dicts(self, data):
	"""function accepts an open data file of tab delimeted, xml-like data. Function
	returns a list of dictionaries, one dictionary for each Part/Variant. Rxns
	are loaded with key: 'meta' and value: list of dictionaries, one dictionary for
	each reaction.
	"""
	dicts = []
	for line in data:
		line = line.strip()
		if line.startswith(">"):
			part_dict = parse_line(line[1:].split("\t"))
		else:
			dicts[-1].setdefault('meta', [])
			meta = parse_line(line.strip("\t").split("\t"))
			dicts[-1]['meta'].append(meta)
	return dicts

def combine(s1, s2):
	"""Function accepts two lists of objects that accept concatenation, and returns all
	combinations of list one with list two. If either list is empty, then the other list
	is returned. If both lists are empty, an empty list is returned.
	"""
	if len(s1) == 0:
		return s2
	if len(s2) == 0:
		return s1
	combs = []
	for e1 in s1:
		for e2 in s2:
			combs.append(e1 + e2)
	return combs

class connection():
	
	def __init__(self, sqldb, mongodb):
		engine = dbmap.connect(sqldb, mongodb)
		Session = sessionmaker(bind=engine)
		self.session = Session()
		
	def dict_sql_match(self, match_type, match_dict, match_columns = []):
		"""Function accepts a match type (object defined in the dbmap orm) and a match
		dict. Function returns the first object of the type match_type that matches all
		of the columns values (specified as key value pairs) in the match_dict.
		If you want to match a subset of the columns in the match_dict, you can include
		them in the list match_columns. If match_columns includes a column name that
		the match_dict doesn't have a key-value pair for, it is ignored.
		"""
		if not match_columns:
			match_columns = match_dict.keys()
		
		q = self.session.query(match_type)
		for column in set(match_columns) & set(match_dict.keys()):
			q = q.filter(column == match_dict[column])
		return q.first()
		
	def dict_sql_match_one(self, match_type, match_dict, match_columns):
		"""Function accepts a match type (object defined in the dbmap orm) and a match
		dict. Function returns the first object of the type match_type that matches one
		of the columns values (specified as key value pairs) in the match_dict). The 
		match_columns list is the order that columns should be checked for a match. 
		A column and its corresponding value from match_dict is checked in the database.
		If there is a match, it is returned. If there is no match, the next column in the
		match_columns list is queried. If none of them have matches, then a None object
		is returned.
		"""
		match_columns = [col for col in match_columns if col in match_dict.keys()]
		if not match_columns:
			return None
		for column in match_columns:
			m = self.session.query(match_type).filter(column == match_dict[column]).first()
			if m:
				return m
			else:
				continue
		return m

	def add_rxns(self, rxn_dicts):
		rxns = []
		for rxn_dict in rxn_dicts:
			rxn_dict.setdefault('cosubstrate', None)
			rxn_dict.setdefault('coproduct', None)
			
			rxn = self.dict_sql_match(dbmap.Rxn, rxn_dict, ['substrate', 'product', 'cosubstrate', 'coproduct'])
			if not rxn:
				rxn = dbmap.Rxn(rxn_dict)
			else:
				rxn.reset(rxn_dict)
			self.session.add(rxn)
			self.session.flush()
			rxns.append(rxn)
		return rxns
		
	def import_parts(self, data):
		"""Function accepts a list of parts represented as dictionaries, and imports
		them into the db. 
		"""
		for part_dict in data:
			rxns = part_dict.pop('meta')
			part = self.dict_sql_match_one(dbmap.Part, part_dict, ['did', 'dna', 'name'])
			new = False
			if not part:
				part = dbmap.Part(part_dict)
				new = True
			else:
				part.reset(part_dict)
			if new:
				self.session.add(part)
			part.rxns.extend(add_rxns(rxns))
			self.session.flush()
				
	def import_designs(self, data):
		"""Function accepts a list of design specs represented as dictionaries, and
		imports them into the db.
		"""
		for design_dict in data:
			rxns = design_dict.pop('meta')
			design = self.dict_sql_match_one(dbmap.Design, design_dict, ['spec'])
			new = False
			if not deisgn:
				design = dbmap.Design(design_dict)
				new = True
			else:
				design.reset(design_dict)
			if new:
				self.session.add(design)
			for path in self.spec_parse(design.spec):
					
		
		table = [[ele for ele in row.split("\t")] for row in data.read()]
		labels = table[0]
		for row in table[1:]:
			doc = dict(zip(labels,row))
			for path in self.spec_parse(doc['spec']):
				db.designs.insert(dict(zip(labels, row) + ['parts', path]))
	
	def part_lookup(self, spec):
		spec = spec.split(";")
		atts = dict([att.split(":") for att in spec])
		return self.dict_sql_match(dbmap.Part, atts)
	
	def spec_parse(self, spec):
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
						return combine(self.spec_parse(spec[:len(spec)-ind-1]), self.spec_parse(rec)) +\
								combine(self.spec_parse(rec), self.spec_parse(spec[:len(spec)-ind-1]))
					if char == "+":
						return combine(self.spec_parse(spec[:len(spec)-ind-1]), self.spec_parse(rec))
					else:
						return 'uh oh'
		
		plus = spec.rfind("+")
		star = spec.rfind("*")
		if star > plus:
			l, r = spec.rsplit("*", 1)
			return combine(self.spec_parse(l), self.spec_parse(r)) +\
					combine(self.spec_parse(r), self.spec_parse(l))
		if plus > star:
			l, r = spec.rsplit("+", 1)
			return combine(self.spec_parse(l), self.spec_parse(r))
		if star == plus == -1:
			return self.part_lookup(spec)
				
	
				
	def load_db(self):
		self.parts = self.session.query(dbmap.Part).all()
		self.rxns = self.session.query(dbmap.Rxn).all()
		self.designs = self.session.query(dbmap.Design).all()
		self.tubes = self.session.query(dbmap.Tube).all()
		self.constructs = self.session.query(dbmap.Construct).all()