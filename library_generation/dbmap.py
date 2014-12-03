

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
NullPool = sqlalchemy.pool.NullPool

MongoClient = pymongo.MongoClient

#Open declarative base builder to build object-db mappers
Base = declarative_base()
	
ct = Table('constructs_tubes', Base.metadata,
	Column('construct_id', Integer, ForeignKey('constructs.cid'), primary_key=True),
	Column('tube_id', Integer, ForeignKey('tubes.tid'), primary_key=True)
	)
	
def connect(sqldb, mongodb, verb=False):
	engine = create_engine('sqlite:///'+sqldb)
	metadata = Base.metadata
	metadata.create_all(engine)
	Session = sessionmaker(bind=engine)
	global session
	session = Session()
	
	global mdb
	client = MongoClient('localhost', 27017)
	mdb = client[mongodb]
	if mongodb=='temp':
		mdb.connection.drop_database('temp')
		mdb = client[mongodb]
	return session, mdb
	
def close_all():
	session.close()
	mdb.close()
	
def sync_mongo(collection, atts):
	return list(mdb[collection].find(atts))[0]
	
def update_mongo(collection, id, atts):
	mdb[collection].update(id, atts)
	return list(mdb[collection].find(id))[0]
	
def get_mdb_ids(collection, atts, rfields={'_id':1}):
	return list(mdb[collection].find(atts, rfields))
	
class Part(Base):

	__tablename__ = 'parts'
	
	pid = Column(Integer, primary_key=True)
	ec = Column(String)
	host = Column(String)
	protein = Column(String)
	dna = Column(String, unique=True, nullable=True)
	length = Column(Integer)
	type = Column(String)
	accession = Column(String)
	enzyme = Column(String)
	name = Column(String, unique=True, nullable=True)
	version = Column(String)
	
	#list of designs that contain this part
	#designparts = relationship with designs that contain this part.
	
	#tubes variable created by one-to-many relationship specified in Tube()
	#tubes = relationship(Tube, backref=backref('part'))
	
	defined = ['pid', 'ec', 'host', 'protein', 'dna', 'type', 'accession', 'enzyme',
			'name', 'version', 'length']
	
	def __init__(self, dict):
		self.mon = {}
		self.reset(dict)
		session.add(self)
		session.flush()
		self.insert_mon()
	
	def reset(self, dict):
		for k, v in dict.items():
			if k in self.defined:
				setattr(self, k, v)
			else:
				self.mon[k] = v
		if self.dna:
			self.length = len(self.dna)
		
	@sqlalchemy.orm.reconstructor
	def get_mon(self):
		self.mon = sync_mongo('parts', {'pid': self.pid})
		
	def insert_mon(self):
		self.mon['pid'] = self.pid
		mdb['parts'].insert(self.mon)
		
	def push_mon(self):
		self.mon = update_mongo('parts', {'pid': self.pid}, self.mon)


class Design(Base):
	
	__tablename__ = 'designs'
	
	did = Column(Integer, primary_key=True)
	type = Column(String)
	spec = Column(String)
	name = Column(String)
	
	#DesignPart variable created by many-to-many relationship specified in Part()
	#designparts = list of links to parts in the design
	
	#constructs variable created by one-to-many relationship specified in Construct()
	#constructs = list of constructs that are specified by this design
	
	defined = ['did', 'type', 'spec', 'name']
	
	def __init__(self, dict):
		self.mon = {}
		self.reset(dict)
		session.add(self)
		session.flush()
		self.insert_mon()
	
	def reset(self, dict):
		
		for k, v in dict.items():
			if k in self.defined:
				setattr(self, k, v)
			else:
				self.mon[k] = v
	
	@sqlalchemy.orm.reconstructor
	def get_mon(self):
		self.mon = sync_mongo('designs', {'did': self.did})
		
	def insert_mon(self):
		self.mon['did'] = self.did
		mdb['designs'].insert(self.mon)
		
	def push_mon(self):
		self.mon = update_mongo('designs', {'did': self.did}, self.mon)
	
	def plot(self):
		return self.name, [designpart_plot.plot() for designpart_plot in self.designparts]
		
class DesignPart(Base):
	
	__tablename__ = 'designparts'
	
	dpid = Column(Integer, primary_key=True)
	pid = Column(Integer, ForeignKey('parts.pid'))
	did = Column(Integer, ForeignKey('designs.did'))
	
	start = Column(Integer)
	end = Column(Integer)
	direction = Column(Integer)
	
	part = relationship('Part', backref=backref('designparts'))
	
	design = relationship('Design', backref=backref('designparts'))
	
	defined = ['dpid', 'pid', 'did', 'start', 'end', 'direction']
	
	def __init__(self, dict):
		self.mon = {}
		self.reset(dict)
		session.add(self)
		session.flush()
		self.insert_mon()
	
	def reset(self, dict):
		for k, v in dict.items():
			if k in self.defined:
				setattr(self, k, v)
			else:
				self.mon[k] = v
			
	@sqlalchemy.orm.reconstructor
	def get_mon(self):
		self.mon = sync_mongo('designparts', {'dpid': self.dpid})
		
	def insert_mon(self):
		self.mon['dpid'] = self.dpid
		mdb['designparts'].insert(self.mon)
		
	def push_mon(self):
		self.mon = update_mongo('designparts', {'dpid': self.dpid}, self.mon)
			
	def plot(self):
		fwd = True if self.direction == 1 else False
		return {'name': self.part.name, 'start': self.start, 'end': self.end, 'fwd': fwd,
				'type': self.part.type, 'opts': self.part.mon}

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
	
	defined = ['tid', 'pid', 'dna', 'author', 'date', 'plate', 'well', 'location', 'name']
	
	def __init__(self, dna, dict = {}):
		self.dna = dna
		self.reset(dict)
	
	def reset(self, dict):
		for k, v in dict.items():
			if k in self.defined:
				setattr(self, k, v)
			else:
				self.mon[k] = v
	
	@sqlalchemy.orm.reconstructor
	def get_mon(self):
		self.mon = sync_mon('tubes', {'tid': self.tid})
		
	def insert_mon(self):
		self.mon['tid'] = self.tid
		mdb['tubes'].insert(self.mon)
		
	def push_mon(self):
		self.mon = update_mongo('tubes', {'tid': self.tid}, self.mon)

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
			if k in self.defined:
				setattr(self, k, v)
			else:
				self.mon[k] = v
	
	@sqlalchemy.orm.reconstructor
	def get_mon(self):
		self.mon = sync_mon('constructs', {'cid': self.cid})
		
	def insert_mon(self):
		self.mon['cid'] = self.cid
		mdb['constructs'].insert(self.mon)
		
	def push_mon(self):
		self.mon = update_mongo('constructs', {'cid': self.cid}, self.mon)
	