import os
import re
import csv
import sqlite3
from sqlite3 import Error


def create_connection(db_file):
	conn = None
	try:
		conn = sqlite3.connect(db_file)
	except Error as e:
		print(e)
	return conn


def extract_data(target_dir, conn, phasers):
	if not os.path.isdir(target_dir):
		conn.close()
		exit()
	else:
		genes = []
		for dirpath in os.listdir(target_dir):
			#find each gene directory
			check = re.search("((\w+)_chr\d+_\d+_\d+_(\w+))\.assembly_dir$",dirpath)
			if(check):
				#obtain sample number and gene
				info = {
					"sample" : check.group(2),
					"gene" : check.group(3)
				}
				found = False
				for phaser in phasers:
					info["phaser"] = phaser
					if phaser == 'sniffles' or phaser == 'longshot':
						alleles = [1,2]
						for allele in alleles:
							file = phaser+"_hap"+str(allele)
							info["allele"] = allele
							found = add_info(target_dir,conn,file,allele,check,info)
					else:
						allele = 0
						info["allele"] = allele
						found = add_info(target_dir,conn,phaser,allele,check,info)
				if found:
					genes.append(info["gene"])
		return info["sample"], genes

def add_motif_match(cur,info,motif):
	count = info["count"]
	info["motif_id"] = motif[0]
	pathMin = int(motif[2])
	pathMax = int(motif[3])
	if pathMin is not None and pathMax is not None:
		if count >= pathMin and count <= pathMax: #need to handle when there is no max
			info["cat"] = "Likely Pathogenic"
		elif count < pathMin or count > pathMax:
			info["cat"] = "Likely Non-Pathogenic"
		else: #change criteria for an abnormal classification
			info["cat"] = "Abnormal"
	elif pathMin is not None:
		if count >= pathMin:
			info["cat"] = "Likely Pathogenic"
		else:
			info["cat"] = "Likely Non-Pathogenic"
	elif pathMax is not None:
		if count <= pathMax:
			info["cat"] = "Likely Pathogenic"
		else:
			info["cat"] = "Likely Non-Pathogenic"
	else:
		info["cat"] = "Unknown"
	cur.execute("INSERT INTO Result (sample,gene,allele,phaser,pattern,motif,numRepeats,category) VALUES (:sample,:gene_num,:allele,:phaser,:found_pattern,:motif_id,:count,:cat)", info)

def add_info(target_dir,conn,file,allele,check,info):
	gene_found = False
	cur = conn.cursor()
	#find tsv data file for unphased consensus
	tsv_filepath = target_dir+'/'+check.group(0)+'/'+check.group(1)+'_'+file+'.fa.annotated_str.tsv'
	if os.path.isfile(tsv_filepath):
		gene_found = True #need to consider the case the gene is found but no data in the file
		#find the gene in the database
		cur.execute("SELECT id, name FROM Gene WHERE name LIKE :gene;", info)
		gene_data = cur.fetchone()
		info["gene_num"] = gene_data[0]
		#extract data from tsv file
		tsv_file = open(tsv_filepath)
		read_tsv = csv.reader(tsv_file, delimiter='\t')
		for row in read_tsv:
			if(row[0] == 'motif'):
				continue
			#get found pattern, reverse inverse and count
			info["pattern"] = row[0].split('::')[0]
			info["rev_inv_pattern"] = row[0].split('::')[1]
			info["count"] = int(row[2])
			#get all motifs associated with the gene
			cur.execute("SELECT id, pattern, pathMin, pathMax, gene from Motif WHERE gene =:gene_num", info)
			motifs = cur.fetchall()
			found = 0
			for motif in motifs:
				#generate all possibile cyclical variations of the motif
				variations = generate_variations(motif[1])
				#find a match
				if info["pattern"] in variations:
					found = 1
					info["found_pattern"] = info["pattern"]
					add_motif_match(cur,info,motif)
				if info["rev_inv_pattern"] in variations:
					found = 1
					info["found_pattern"] = info["rev_inv_pattern"]
					add_motif_match(cur,info,motif)
			#if no match found, insert separately
			if found == 0:
				info["cat"] = "Unknown"
				cur.execute("INSERT INTO Result (sample,gene,allele,phaser,pattern,numRepeats,category) VALUES (:sample,:gene_num,:allele,:phaser,:pattern,:count,:cat)", info)
		tsv_file.close()
		conn.commit()
	return gene_found

def get_sample_results(sample, conn):
	cur = conn.cursor()
	cur.execute("SELECT g.name, r.pattern, m.pattern, m.pathMin, m.pathMax, r.numRepeats, r.allele, r.phaser, r.category FROM Result r LEFT JOIN Gene g ON r.gene = g.id LEFT JOIN Motif m ON r.motif = m.id WHERE r.sample LIKE ?", (sample,))
	results = cur.fetchall()
	return results

def generate_variations(pattern):
	patterns = []
	length = len(pattern)
	bases = list(pattern)
	i = 1
	while i <= length:
		patterns.append("".join(bases[i:]+bases[:i]))
		i = i + 1
	return patterns

def find_samples(sample, conn):
	cur = conn.cursor()
	cur.execute("SELECT DISTINCT r.sample FROM Result r WHERE r.sample LIKE ?", ('{}%'.format(sample),))
	results = cur.fetchall()
	return results

def main():
	conn = create_connection('./database.db')
	target_dir = '/home/alyne/Documents/Thesis/examples_for_alyne/MBXM037256'
	extract_data(target_dir, conn)
	#get_sample_results('MBXM037256', conn)
	conn.close()

if __name__ == '__main__':
	main()
