



import sqlalchemy
import sqlalchemy.orm
import sqlalchemy.ext.declarative
import os
import re

import pymongo

#In order to parallelize with pp, you cannot have from ... import clauses, so I import sqlalchemy
# and then reassign different functions within it:

declarative_base = sqlalchemy.ext.declarative.declarative_base
create_engine = sqlalchemy.create_engine
ForeignKey = sqlalchemy.ForeignKey
Table = sqlalchemy.Table
Column = sqlalchemy.Column
Integer = sqlalchemy.Integer
String = sqlalchemy.String
MetaData = sqlalchemy.MetaData
Float = sqlalchemy.Float
sessionmaker = sqlalchemy.orm.sessionmaker
relationship = sqlalchemy.orm.relationship
backref = sqlalchemy.orm.backref
mapper = sqlalchemy.orm.mapper

MongoClient = pymongo.MongoClient

#Open declarative base builder to build object-db mappers
Base = declarative_base()

pr = Table('parts_rxns', Base.metadata,
	Column('part_id', Integer, ForeignKey('parts.pid'), primary_key=True),
	Column('rxn_id', Integer, ForeignKey('rxns.rid'), primary_key=True)
	)
	
vr = Table('designs_rxns', Base.metadata,
	Column('design_id', Integer, ForeignKey('designs.did'), primary_key=True),
	Column('rxn_id', Integer, ForeignKey('rxns.rid'), primary_key=True)
	)
	
ct = Table('constructs_tubes', Base.metadata,
	Column('construct_id', Integer, ForeignKey('constructs.cid'), primary_key=True),
	Column('tube_id', Integer, ForeignKey('tubes.tid'), primary_key=True)
	)
	
def connect(sqldb, mongodb, verb=False):
	engine = create_engine('sqlite:///'+sqldb, echo=verb)
	metadata = Base.metadata
	metadata.create_all(engine)
	global mdb
	client = MongoClient('localhost', 27017)
	mdb = client[mongodb]
	return engine
	
def syn_mongo(collection, atts):
	return list(mdb[collection].find(atts))

class Part(Base):

	__tablename__ = 'parts'
	
	pid = Column(Integer, primary_key=True)
	ec = Column(String)
	host = Column(String)
	protein = Column(String)
	dna = Column(String, unique=True, nullable=True)
	type = Column(String)
	accession = Column(String)
	enzyme = Column(String)
	name = Column(String, unique=True, nullable=True)
	version = Column(String)
	
	#list of designs that contain this part
	
	#list of reactions that this part performs
	rxns = relationship('Rxn', secondary=pr, backref=backref('parts'))
	
	#tubes variable created by one-to-many relationship specified in Tube()
	#tubes = relationship(Tube, backref=backref('part'))
	
	def __init__(self, dict):
		self.reset(dict)
	
	def reset(self, dict):
		for k, v in dict.items():
			setattr(self, k, v)
		
	@sqlalchemy.orm.reconstructor
	def get_mongo(self):
		self.mongo = sync_mongo('parts', {'pid': self.pid})

class Design(Base):
	
	__tablename__ = 'designs'
	
	did = Column(Integer, primary_key=True)
	type = Column(String)
	spec = Column(String)
	name = Column(String, unique=True, nullable=True)
	
	#parts variable created by many-to-many relationship specified in Part()
	#parts = list of parts in the design
	
	#list of reactions that this design performs
	rxns = relationship('Rxn', secondary=vr, backref=backref('designs'))
	
	#constructs variable created by one-to-many relationship specified in Construct()
	#constructs = list of constructs that are specified by this design
	
	def __init__(self, dict):
		self.reset(dict)
	
	def reset(self, dict):
		for k, v in dict.items():
			setattr(self, k, v)
	
	@sqlalchemy.orm.reconstructor
	def get_mongo(self):
		self.mongo = sync_mongo('designs', {'did': self.vid})
		
class DesignPart(Base):
	
	__tablename__ = 'designparts'
	
	pid = Column(Integer, ForeignKey('parts.pid'), primary_key = True)
	did = Column(Integer, ForeignKey('designs.did'), primary_key = True)
	
	start = Column(Integer)
	end = Column(Integer)
	direction = Column(Integer)
	
	parts = relationship('Part', backref=backref('designs'))
	
	designs = relationship('Design', backref=backref('parts'))
	
	def __init__(self, dict):
		self.reset(dict)
	
	def reset(self, dict):
		for k, v in dict.items():
			setattr(self, k, v)

class Rxn(Base):

	__tablename__ = 'rxns'
	
	rid = Column(Integer, primary_key=True)
	substrate = Column(String, nullable=False)
	product = Column(String, nullable=False)
	cosubstrate = Column(String)
	coproduct = Column(String)
	name = Column(String, unique=True, nullable=True)
	
	#parts variable created by many-to-many relationship specified in Part()
	#parts = list of parts that perform this reaction
	
	#designs variable created by many-to-many relationship specified in Design()
	#designs = list of designs that perform this reaction
	
	def __init__(self, dict):
		self.reset(dict)
	
	def reset(self, dict):
		for k, v in dict.items():
			setattr(self, k, v)

class Tube(Base):

	__tablename__ = 'tubes'
	
	tid = Column(Integer, primary_key=True)
	pid = Column(Integer, ForeignKey('parts.pid'))
	dna = Column(String, unique=True)
	author = Column(String)
	date = Column(String)
	plate = Column(String)
	well = Column(String)
	location = Column(String)
	name = Column(String, unique=True, nullable=True)
	
	#constructs variable created by many-to-many relationship specified in construct()
	#constructs = list of constructs that contain this tube
	
	#part variable created by one-to-many relationship specified in Tube()
	part = relationship('Part', backref=backref('tubes'))
	
	def __init__(self, dna, dict = {}):
		self.dna = dna
		
		self.reset(dict)
	
	def reset(self, dict):
		for k, v in dict.items():
			setattr(self, k, v)
	
	@sqlalchemy.orm.reconstructor
	def get_mongo(self):
		self.mongo = sync_mongo('tubes', {'tid': self.tid})

class Construct(Base):

	__tablename__ = 'constructs'
	
	cid = Column(Integer, primary_key=True)
	did = Column(Integer, ForeignKey('designs.did'))
	dna = Column(String)
	author = Column(String)
	date = Column(String)
	plate = Column(String)
	well = Column(String)
	location = Column(String)
	name = Column(String, unique=True, nullable=True)
	
	#design template for construct
	template = relationship('Design', backref=backref('constructs'))
	
	#tubes (physical version of parts) included in this construct
	tubes = relationship('Tube', secondary=ct, backref=backref('constructs'))
	
	def __init__(self, dna, dict = {}):
		self.dna = dna
		
		self.reset(dict)
	
	def reset(self, dict):
		for k, v in dict.items():
			setattr(self, k, v)

	@sqlalchemy.orm.reconstructor
	def get_mongo(self):
		self.mongo = sync_mongo('constructs', {'cid': self.cid})	
	