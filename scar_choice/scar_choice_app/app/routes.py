from flask import Flask, render_template, request
from forms import ScarChoice
import scar_choice as sc

app = Flask(__name__)

app.secret_key = 'development key'
 
@app.route('/')
def home():
	form = ScarChoice()
	if request.method == 'POST':
		return render_template('home.html', success=True)

	elif request.method == 'GET':
		return render_template('home.html', form=form)

if __name__ == '__main__':
	app.run(debug=True)
