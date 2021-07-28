from flask import Flask, render_template, request, redirect
import sqlite3
from flask import g
from sqlite3 import Error
from data_reader import extract_data, add_info, test
import os

app = Flask(__name__)

DATABASE = './database.db'

def get_db():
    db = getattr(g, '_database', None)
    if db is None:
        db = g._database = sqlite3.connect(DATABASE)
    return db

def init_db():
    with app.app_context():
        db = get_db()
        with app.open_resource('schema.sql', mode='r') as f:
            db.cursor().executescript(f.read())
        db.commit()

@app.teardown_appcontext
def close_connection(exception):
    db = getattr(g, '_database', None)
    if db is not None:
        db.close()

@app.route("/")
def home(name=None):
	conn = get_db()
	name = test('/home/alyne/Documents/Thesis/examples_for_alyne/MBXM037256', conn)
	return render_template("home.html", name=name)

@app.route("/add-data", methods=["GET","POST"])
def add_data(error=None,genes=None):
	if request.method == "POST":
		target_dir = request.form.get('target_dir')
		if not os.path.isdir(target_dir):
			return render_template("add_sample_data.html", error="Path does not exist")
		conn = get_db()
		sample, genes = extract_data(target_dir,conn)
		return render_template("add_sample_data.html", sample=sample, genes=genes)
	return render_template("add_sample_data.html")

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=105)
