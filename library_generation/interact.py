
import dbmap
from sqlalchemy.orm import sessionmaker
import sys

import dnaplotlib as dpl
import plot_SBOL_designs as plsbol

def parse_data(f):
	"""Function accepts an open file object or read in text from a tab delimited text
	file. Returns a list of dictionaries, one dictionary for each row, with the keys
	being the column headings, and the values being the data for the relevant row.
	"""
	labels = f.readline().strip().split("\t")
	rows = []
	for line in f:
		rows.append(dict([(label, value) for label, value in zip(labels, line.strip().split("\t"))]))
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
		engine, self.mdb = dbmap.connect(sqldb, mongodb)
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
			col = getattr(match_type, column)
			q = q.filter(col == match_dict[column])
		return q.all()
		
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
			print part_dict
			part = self.dict_sql_match_one(dbmap.Part, part_dict, ['pid', 'dna', 'name'])
			new = False
			if not part:
				part = dbmap.Part(part_dict)
				new = True
			else:
				part.reset(part_dict)
			if new:
				self.session.add(part)
			self.session.flush()
			part.insert_mon()
				
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
				self.session.add(design)
				self.session.flush()
				design.insert_mon()
				
	
	def part_lookup(self, spec):
		spec = spec.split(";")
		atts = dict([att.split(":") for att in spec])
		parts = self.dict_sql_match(dbmap.Part, atts)
		return [[[part, 1]] for part in parts]
	
	def spec_parse(self, spec):
		"""function accepts a string representing a 'spec'. This includes (), *, +, -, and data
		specifying documents in the database. Example:
			A+B = AB
			-A+B = reverse complement of A, plus B
			-(A+B) = reverse complment of AB
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
						else:
							return combs
					if char == "+":
						combs = combine(self.spec_parse(spec[:-(ind+1)]), self.spec_parse(rec))
						if direction == -1:
							return reverse(combs)
						else:
							return combs
					else:
						sys.stderr('You broke it. Smooth. Something is not quite right with the parantheses')
						return
		
		plus = spec.find("+")
		star = spec.find("*")
		if star > plus:
			l, r = spec.split("*", 1)
			combs = combine(self.spec_parse(l), self.spec_parse(r)) +\
					combine(self.spec_parse(r), self.spec_parse(l))
		if plus > star:
			l, r = spec.split("+", 1)
			return combine(self.spec_parse(l), self.spec_parse(r))
		if star == plus == -1:
			if spec.startswith("-"):
				combs = self.part_lookup(spec[1:])
				return reverse(combs)
			else:
				return self.part_lookup(spec)
				
	def load_db(self):
		self.parts = self.session.query(dbmap.Part).all()
		self.designs = self.session.query(dbmap.Design).all()
		self.tubes = self.session.query(dbmap.Tube).all()
		self.constructs = self.session.query(dbmap.Construct).all()
		
	def plot_designs(self, plot_parameters={'linewidth':1.5,'scale':1.0,'fig_x':8.0,
			'fig_y':2.0,'show_title':'N','backbone_pad_left':3.0,'backbone_pad_right':3.0}):
		dr = dpl.DNARenderer()
		reg_renderers = dr.std_reg_renderers()
		part_renderers = dr.SBOL_part_renderers()
		plot_parameters = {'linewidth':1.5, 'scale':1.0, 'fig_x':8.0, 'fig_y':2.0 'show_title':'N', 
				'backbone_pad_left':3.0, 'backbone_pad_right':3.0}
		
		dna_designs = self.designs
		
		if 'axis_y' not in plot_params.keys():
			plot_params['axis_y'] = 55
		left_pad = 0.0
		right_pad = 0.0
		scale = 1.0
		linewidth = 1.0
		if 'backbone_pad_left' in plot_params.keys():
			left_pad = plot_params['backbone_pad_left']
		if 'backbone_pad_right' in plot_params.keys():
			right_pad = plot_params['backbone_pad_right']
		if 'scale' in plot_params.keys():
			scale = plot_params['scale']
		if 'linewidth' in plot_params.keys():
			linewidth = plot_params['linewidth']
		dr = dpl.DNARenderer(scale=scale, linewidth=linewidth,
							 backbone_pad_left=left_pad, 
							 backbone_pad_right=right_pad)

		# We default to the standard regulation renderers
		reg_renderers = dr.std_reg_renderers()

		# We default to the SBOL part renderers
		part_renderers = dr.SBOL_part_renderers()

		# Create the figure
		fig = plt.figure(figsize=(plot_params['fig_x'],plot_params['fig_y']))
		# Cycle through the designs an plot on individual axes

		design_list = sorted(dna_designs.keys())
		if(regs_info != None):
			regs_list   = sorted(regs_info.keys())
	
		num_of_designs = len(design_list)
		#print len(design_list),len(regs_list)
		ax_list = []
		max_dna_len = 0.0
		for i in range(num_of_designs):
			# Create axis for the design and plot

			regs = None
			if(regs_info != None):
				regs   =  regs_info[i]
			design =  dna_designs[design_list[i]]

			ax = fig.add_subplot(num_of_designs,1,i+1)
			if 'show_title' in plot_params.keys() and plot_params['show_title'] == 'Y':
				ax.set_title(design_list[i], fontsize=8)
			start, end = dr.renderDNA(ax, design, part_renderers, regs, reg_renderers)

			dna_len = end-start
			if max_dna_len < dna_len:
				max_dna_len = dna_len
			ax_list.append(ax)
		for ax in ax_list:
			ax.set_xticks([])
			ax.set_yticks([])
			# Set bounds
			ax.set_xlim([(-0.01*max_dna_len)-left_pad,
						max_dna_len+(0.01*max_dna_len)+right_pad])

			ax.set_ylim([-plot_params['axis_y'],plot_params['axis_y']])

			ax.set_aspect('equal')
			ax.set_axis_off()
		# Save the figure
		plt.tight_layout()
		fig.savefig(out_filename, transparent=True)
		# Clear the plotting cache
		plt.close('all')
		