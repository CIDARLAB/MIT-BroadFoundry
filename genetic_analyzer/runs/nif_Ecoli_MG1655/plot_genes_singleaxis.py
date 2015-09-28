#!/usr/bin/env python

import csv
import random
import dnaplotlib as dpl
import matplotlib.pyplot as plt

parts_per_axis = 50;

def rand_rgb_string():
	r = 1.0*random.randint(1,100)/100
	g = 1.0*random.randint(1,100)/100
	b = 1.0*random.randint(1,100)/100
	rgb = str(r) + ";" + str(g) + ";" + str(b)
	return rgb

def rand_rgb():
	r = 1.0*random.randint(1,100)/100
	g = 1.0*random.randint(1,100)/100
	b = 1.0*random.randint(1,100)/100
	rgb = (r,g,b)
	return rgb

#part_name,type,x_extent,y_extent,start_pad,end_pad,color,hatch,arrowhead_height,arrowhead_length,linestyle,linewidth

parts_filename = 'Stetzurigenome.csv'

parts = [];
parts_reader = csv.reader(open(parts_filename, 'rU'), delimiter=',')

a = 0
b = 45
c = 0

for row in parts_reader:
	c = c + 1
	if(c > a and c < b):
		part = {}
		part['name']     = row[0]
		part['bp_start'] = int(row[1])
		part['bp_end']   = int(row[2])
		part['strand']   = row[3]
		part['color']    = rand_rgb()
		parts.append(part)
	

dna_designs = []

counter = 0;

for part in parts:

	part_name = part['name']
	part_type = "CDS"
	x_extent = int((part['bp_end']-part['bp_start'])/40)
	y_extent = 3
	start_pad = ""
	end_pad = ""
	color = part['color']
	hatch = ""
	arrowhead_height = 2
	arrowhead_length = 8
	linestyle = "-"
	linewidth = 1

	fwd = False
	if(part['strand'] == "+"):
		fwd = True

	if(x_extent < arrowhead_length):
		arrowhead_length = x_extent

	#fields = [part_name,part_type,x_extent,y_extent,start_pad,end_pad,color,hatch,arrowhead_height,arrowhead_length,linestyle,linewidth]
	#myString = ','.join(map(str, fields)) 
	#print myString

	#label_rotation
	#label_y_offset
	#label_size

	dpl_part  = {'type':part_type, 'name':part_name, 'fwd':fwd, 'opts':{'color':color, 'x_extent':x_extent, 'y_extent':y_extent, 'arrowhead_length':arrowhead_length, 'arrowhead_height':arrowhead_height, 'start_pad':0, 'end_pad':0, 'label':part_name, 'label_size':9, 'label_y_offset':-10, 'label_rotation':90}}

	if(counter > 0):
		prev_part = parts[counter-1]
		x_extent = int((part['bp_start']-prev_part['bp_end'])/40)
		dpl_spacer = {'type':'EmptySpace', 'name':'S'+str(counter), 'fwd':True, 'opts':{'x_extent':x_extent}}
		dna_designs.append(dpl_spacer)

	dna_designs.append(dpl_part)
	counter = counter + 1



fig = plt.figure(figsize=(15,2))	

ax  = fig.add_subplot(1,1,1, axisbg='white')

dr = dpl.DNARenderer()
part_renderers = dr.SBOL_part_renderers()

start, end = dr.renderDNA(ax, dna_designs, part_renderers)

ax.set_xlim([start, end])
ax.set_ylim([-15,15])
#ax.set_aspect('equal')
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')


plt.subplots_adjust(hspace=0.01, wspace=0.01, left=0.05, right=0.95, top=0.99, bottom=0.01)

fig.savefig('output_gene_single.pdf', dpi=300)

# Clear the plotting cache
plt.close('all')
