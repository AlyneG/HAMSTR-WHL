from flask import Flask, render_template, request, redirect, g, url_for
import sqlite3
from sqlite3 import Error
from data_reader import *
import os
import json

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
def home(database=None):
	conn = get_db()
	database = DATABASE
	return render_template("home.html", database=database)

@app.route("/add-data", methods=["GET","POST"])
def add_data(error=None,genes=None):
	if request.method == "POST":
		target_dir = request.form.get('target_dir')
		if not os.path.isdir(target_dir):
			return render_template("add_sample_data.html", error="Path does not exist")
		conn = get_db()
		sample, genes = extract_data(target_dir,conn,['unphased','longshot','sniffles'])
		return render_template("add_sample_data.html", sample=sample, genes=genes)
	return render_template("add_sample_data.html")

@app.route("/search", methods=["POST","GET"])
def search():
	if request.method == 'POST':
		sample = request.form['sample']
		if sample:
			# find how many matches there are
			conn = get_db()
			matches = find_samples(sample, conn)
			num_matches = len(matches)
			# if none, display error
			if(num_matches == 0):
				error = "No results found matching the term "+sample
				return render_template("search_page.html", error = error)
			# if only one, return that sample page
			if(num_matches == 1):
				return redirect(url_for('get_sample', sample = matches[0]))
			# else display options with links
			if(num_matches > 1):
				return render_template("search_page.html", sample = sample, num_matches = num_matches, matches = matches)
		else:
			error = "Please enter a search term"
			return render_template("search_page.html", error = error)
	else:
		return render_template("search_page.html")

@app.route("/sample/<sample>", methods=["POST","GET"])
def get_sample(sample):
	error = None
	conn = get_db()
	results = get_sample_results(sample, conn)
	if not results:
		error = "No results found for "+sample
	return render_template("sample_page.html", sample=sample, results=results, error=error)

@app.route("/<sample>/<gene>", methods=["POST","GET"])
def get_gene(sample,gene):
	error = None
	conn = get_db()
	results = get_sample_gene_results(sample, gene, conn)
	if not results:
		error = "No results found for "+sample+" "+gene
	return render_template("gene_page.html", sample=sample, gene=gene, results=results, error=error)

@app.route("/<result>", methods=["POST","GET"])
def get_result(result):
	error = None
	motif_lo = 0
	motif_hi = 0
	count = 0
	lower_bound = 0
	upper_bound = 0
	conn = get_db()
	result_info = get_sample_result(result, conn)
	if not result_info:
		error = "No result entry found for "+result
	else:
		motif_lo = result_info[7]
		motif_hi = result_info[8]
		count = result_info[5]
		if motif_lo is not None and motif_hi is not None:
			lower_bound = motif_lo if motif_lo < count else count
			upper_bound = motif_hi if motif_hi > count else count
	return render_template("result_page.html", result=result, result_info=result_info, error=error, motif_lo=json.dumps(motif_lo), motif_hi=json.dumps(motif_hi), count=json.dumps(count), lower_bound=json.dumps(lower_bound), upper_bound=json.dumps(upper_bound))

if __name__ == '__main__':
    app.run(debug=True)
