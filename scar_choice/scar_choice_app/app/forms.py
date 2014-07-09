from flask.ext.wtf import Form
from wtforms import TextAreaField, SubmitField, RadioField, SelectField, IntegerField

class ScarChoice(Form):
	length = SelectField(u'Scar Length: ', choices=[(2, '2 bp'), (3, '3 bp'), (4, '4 bp'), (5, '5 bp')], default=4)
	req_num = IntegerField(u'Quantity: ', default='1')
	homology = SelectField(u'Max Homology: ', choices=[(1, '1 bp'), (2, '2 bp'), (3, '3 bp'), (4, '4 bp')], default=2)
	search_type = SelectField(u'Search Type: ', choices=[(0, 'Random'), (1, 'Enumerate')], default=0)
	exclude_scars = TextAreaField(u'Existing Scars: ')
	submit = SubmitField(u'Generate Scars')
