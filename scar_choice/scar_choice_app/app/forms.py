from flask.ext.wtf import Form
from wtforms import TextAreaField, SubmitField, RadioField, SelectField, IntegerField

class ScarChoice(Form):
	length = SelectField(u'Scar Length: ', choices=[(2, '2 bp'), (3, '3 bp'), (4, '4 bp'), (5, '5 bp')])
	req_num = IntegerField(u'Quantity: ')
	homology = SelectField(u'Max Homology: ', choices=[(2, '2 bp'), (3, '3 bp'), (4, '4 bp')])
	search_type = SelectField(u'Search Type: ', choices=[(0, 'Random'), (1, 'Enumerate')])
	exclude_scars = TextAreaField(u'Scars to Exclude: ')
	submit = SubmitField(u'Generate Scars')
