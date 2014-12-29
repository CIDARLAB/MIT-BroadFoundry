
import dbmap
from sqlalchemy.orm import sessionmaker
import sys

import dnaplotlib as dpl

import sqlalchemy

import re
import itertools


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
	labels = f.readline().strip("\n").split("\t")
	rows = []
	for line in f:
		rows.append(dict([(label, analyze_cell(value)) for label, value in zip(labels, line.strip("\n").split("\t")) if value]))
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
	tuple is (part, direction), and direction is -1 or 1) and gives the reverse
	complement of each sequence. This means the direction in each of the designparts is
	reversed, and each design list itself is also reversed.
	"""
	new_combs = []
	for comb in combs:
		new_comb = []
		for part in comb:
			new_comb.append([part[0], -1*part[1]])
		new_combs.append(new_comb[::-1])
	return new_combs
	

	
def multi_split(str, chars):
	parts = []
	for c in str:
		if c in chars:
			parts.append([c])
			parts.append([])
		else:
			if len(parts) == 0:
				parts.append([])
			parts[-1].append(c)
	return ["".join(part) for part in parts]
					

class connection():
	
	def __init__(self, sqldb = ":memory:", mongodb = "temp"):
		self.session, self.mdb = dbmap.connect(sqldb, mongodb)
	
	def disconnect(self):
		dbmap.close_all()
		
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
		#List of regulators with the source filled out. The target needs to be updated
		# after the first iteration of the for loop in case the target has not been
		# added to the database yet.
		regulators = []
		for part_dict in data:
			neg_reg = part_dict.pop('neg_reg', None)
			pos_reg = part_dict.pop('pos_reg', None)
			try:
				part = dbmap.Part(part_dict)
			except sqlalchemy.exc.IntegrityError:
				self.session.rollback()
				print "WAS NOT UNIQUE:\n", part_dict
			color = None if 'color' not in part_dict else part_dict['color']
			if neg_reg:
				regulators.append((part, -1, neg_reg, color))
			if pos_reg:
				regulators.append((part, 1, pos_reg, color))
		for regulator in regulators:
			for regulated in self.part_lookup(regulator[2]):
				if regulator[3]:
					reg = dbmap.Regulator({'type':regulator[1], 'color':regulator[3]})
				else:
					reg = dbmap.Regulator({'type':regulator[1]})
				reg.source = regulator[0]
				reg.target = regulated
		self.load_db()
				
	def import_designs(self, data):
		"""Function accepts a list of design specs represented as dictionaries, and
		imports them into the db.
		"""
		self.load_db()
		for design_dict in data:
			num = 1
			variants = self.spec_parse(design_dict['spec'])
			for variantlist in variants:
				if num % 1000 == 0:
					print num, ' of ', len(variants)
				num += 1
				design = dbmap.Design(design_dict)
				last = 0
				for part in variantlist:
					designpart = dbmap.DesignPart({'start':last+1, 'end':last+part[0].length})
					designpart.part = part[0]
					designpart.direction = part[1]
					designpart.design = design
					last = last + part[0].length
		self.load_db()
				
	
	def part_lookup(self, spec):
		spec = spec.split(";")
		atts = dict([att.split(":") for att in spec])
		return self.dict_match(dbmap.Part, atts)
		
	def parens(self, spec):
		if len([c for c in spec if c in "()+*-~"]) == 0:
			return [spec]
		new = []
		prev_counter = 0
		counter = 0
		save = True
		for char in multi_split(spec, ")(+*-~"):
			if not char:
				continue
			if char == "(":
				counter += 1
				if prev_counter == 0:
					save = False
					new.append([])
			if char == ")":
				counter -= 1
				if prev_counter == 1:
					save = False
				
			if counter == 0 and prev_counter == 0:
				new.append(char)
			else:
				if save:
					new[-1].append(char)
			
			if not save:
				save = True
			prev_counter = counter
		return ["".join(n) for n in new]
	
	def combine(self, group):
		combs = [list(itertools.chain(*c)) for c in itertools.product(*group)]
		return combs

	def star(self, group):
		#All I have to say is... Sorry. But a list comprehension was easier here.
		perms = [sub for group in list(itertools.permutations([list(g)[0] \
				for k,g in itertools.groupby(group, lambda x: x == '*') if not k])) \
				for sub in self.combine(group)]
		return perms

	def plus(self, group):
		combs = self.combine([list(g)[0] for k,g in itertools.groupby(group,
								lambda x: x == '+') if not k])
		return combs

	def spec_parse(self, spec):
		separated = self.parens(spec)
		new_sep = []
		dir = 1
		for comp in separated:
			if comp == '-':
				dir = -1
				continue
			if comp == '~':
				dir = 0
				continue
			if comp in '*+':
				new_sep.append(comp)
			elif '(' in comp or '*' in comp or '+' in comp:
				if dir == 1:
					new_sep.append(self.spec_parse(comp))
				if dir == -1:
					new_sep.append(reverse(self.spec_parse(comp)))
				if dir == 0:
					new_sep.append(self.spec_parse(comp) + reverse(self.spec_parse(comp)))
				dir = 1
			else:
				if dir == 1:
					new_sep.append([[[p, 1]] for p in self.part_lookup(comp)])
				if dir == -1:
					new_sep.append(reverse([[[p, 1]] for p in self.part_lookup(comp)]))
				if dir == 0:
					stuff = [[[p, 1]] for p in self.part_lookup(comp)]
					new_sep.append(stuff + reverse(stuff))
				dir = 1
					
		separated = new_sep
		
		if '*' in separated:
			return self.star(separated)
		elif '+' in separated:
			return self.plus(separated)
		else:
			return separated[0]
				
	def load_db(self):
		self.parts = self.session.query(dbmap.Part).all()
		self.designs = self.session.query(dbmap.Design).all()
		self.tubes = self.session.query(dbmap.Tube).all()
		self.constructs = self.session.query(dbmap.Construct).all()
