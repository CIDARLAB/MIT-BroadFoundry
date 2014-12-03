
import dbmap
from sqlalchemy.orm import sessionmaker
import sys

import dnaplotlib as dpl

def convert_float(num):
	"""Function accepts text, converts to float if possible.
	"""
	try:
		return float(num)
	except ValueError:
		return num

def analyze_cell(c):
	"""Function accepts a cell of data as text. If it's a number, it gets converted to a
	float, if it is semicolon delimited values, they are converted to a list of floats 
	(if the values are numbers).
	"""
	if ';' in c:
		return [convert_float(ci) for ci in c.split(";")]
	else:
		return convert_float(c)
		

def parse_data(f):
	"""Function accepts an open file object or read in text from a tab delimited text
	file. Returns a list of dictionaries, one dictionary for each row, with the keys
	being the column headings, and the values being the data for the relevant row.
	"""
	labels = f.readline().strip().split("\t")
	rows = []
	for line in f:
		rows.append(dict([(label, analyze_cell(value)) for label, value in zip(labels, line.strip().split("\t")) if value]))
	return rows

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
	
def reverse(combs):
	"""Function accepts a list of designs (which is itself a list of tuples, where each
	tuple is a (part, direction), and direction is -1 or 1) and gives the reverse
	complement of each sequence. This means the direction in each of the designparts is
	reversed, and each design list itself is also reversed.
	"""
	for comb in combs:
		for part in comb:
			part[1] *= -1
	return [part for comb in combs for part in comb[::-1]]		

class connection():
	
	def __init__(self, sqldb, mongodb):
		self.session, self.mdb = dbmap.connect(sqldb, mongodb)
		
	def dict_match(self, match_type, match_dict, match_columns = []):
		"""Function accepts a match type (object defined in the dbmap orm) and a match
		dict. Function returns the all objects of the type match_type that matches all
		of the columns values (specified as key value pairs) in the match_dict.
		If you want to match a subset of the columns in the match_dict, you can include
		them in the list match_columns. If match_columns includes a column name that
		the match_dict doesn't have a key-value pair for, it is ignored.
		"""
		if not match_columns:
			mongo_columns = [col for col in match_dict.keys() if col not in match_type.defined]
			sql_columns = [col for col in match_dict.keys() if col in match_type.defined]
		else:
			mongo_columns = dict([col for col in match_columns if col in match_dict.keys() and col not in match_type.defined])
			sql_columns = [col for col in match_columns if col in match_dict.keys() and col in match_type.defined]
			
		ids = dbmap.get_mdb_ids(match_type.__tablename__, dict([(col, match_dict[col]) for col in mongo_columns]))
				
		q = self.session.query(match_type)
		col = getattr(match_type, match_type.defined[0])
		q.filter(col.in_(ids))
		for column in sql_columns:
			col = getattr(match_type, column)
			q = q.filter(col == match_dict[column])
		return q.all()
		
	def dict_match_one(self, match_type, match_dict, match_columns):
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
			col = getattr(match_type, column)
			m = self.session.query(match_type).filter(col == match_dict[column]).first()
			if m:
				return m
			else:
				continue
		return m
		
	def import_parts(self, data):
		"""Function accepts a list of parts represented as dictionaries, and imports
		them into the db. 
		"""
		for part_dict in data:
			part = dbmap.Part(part_dict)
				
	def import_designs(self, data):
		"""Function accepts a list of design specs represented as dictionaries, and
		imports them into the db.
		"""
		for design_dict in data:
			num = 1
			for variantlist in self.spec_parse(design_dict['spec']):
				design = dbmap.Design(design_dict)
				last = 0
				for part in variantlist:
					designpart = dbmap.DesignPart({'start':last+1, 'end':last+part[0].length})
					designpart.part = part[0]
					designpart.direction = part[1]
					designpart.design = design
					last = last + part[0].length
				
	
	def part_lookup(self, spec):
		spec = spec.split(";")
		atts = dict([att.split(":") for att in spec])
		parts = self.dict_match(dbmap.Part, atts)
		return [[[part, 1]] for part in parts]
	
	def spec_parse(self, spec):
		"""function accepts a string representing a 'spec'. This includes (), *, +, -, ~, and data
		specifying documents in the database. Example:
			A+B = AB
			-A+B = reverse complement of A, plus B
			-(A+B) = reverse complment of AB
			~A+B = AB and reverse complement of A, plus B
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
	
		if spec.startswith("(") or spec.startswith("-("):
			if spec.startswith("-("):
				spec = spec[1:]
				direction = -1
			elif spec.startswith("~("):
				spec = spec[1:]
				direction = 0
			else:
				direction = 1
			recording = True
			recorded = []
			opened = 0
			closed = 0
			for ind, char in enumerate(spec):
				if recording:
					if char == "(":
						opened += 1
					if char == ")":
						closed += 1
					recorded.append(char)
					if opened == closed:
						if ind == len(spec) - 1:
							combs = spec_parse("".join(rec[1:-1]))
							if direction == -1:
								return reverse(combs)
							else:
								return combs
						continue
				else:
					rec = "".join(rec[1:-1])
					if char == "*":
						combs = combine(self.spec_parse(spec[:-(ind+1)]), self.spec_parse(rec)) +\
								combine(self.spec_parse(rec), self.spec_parse(spec[:-(ind+1)]))
						if direction == -1:
							return reverse(combs)
						elif direction == 0:
							return combs + reverse(combs)
						else:
							return combs
					if char == "+":
						combs = combine(self.spec_parse(spec[:-(ind+1)]), self.spec_parse(rec))
						if direction == -1:
							return reverse(combs)
						elif direction == 0:
							return combs + reverse(combs)
						else:
							return combs
					else:
						sys.stderr('You broke it. Smooth. Something is not quite right with the parantheses')
						return
		
		plus = spec.find("+")
		star = spec.find("*")
		if star > plus:
			l, r = spec.split("*", 1)
			return combine(self.spec_parse(l), self.spec_parse(r)) +\
					combine(self.spec_parse(r), self.spec_parse(l))
		if plus > star:
			l, r = spec.split("+", 1)
			print 'going...'
			return combine(self.spec_parse(l), self.spec_parse(r))
		if star == plus == -1:
			if spec.startswith("-"):
				combs = self.part_lookup(spec[1:])
				return reverse(combs)
			elif spec.startswith("~"):
				combs = self.part_lookup(spec[1:])
				return combs + reverse(combs)
			else:
				return self.part_lookup(spec)
				
	def load_db(self):
		self.parts = self.session.query(dbmap.Part).all()
		self.designs = self.session.query(dbmap.Design).all()
		self.tubes = self.session.query(dbmap.Tube).all()
		self.constructs = self.session.query(dbmap.Construct).all()
