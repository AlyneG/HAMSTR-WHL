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


def extract_data(target_dir, conn):
	if not os.path.isdir(target_dir):
		conn.close()
		exit()
	else:
		genes = []
		cur = conn.cursor()
		for dirpath in os.listdir(target_dir):
			#find each gene directory
			check = re.search("((\w+)_chr\d+_\d+_\d+_(\w+))\.assembly_dir$",dirpath)
			if(check):
				#obtain sample number and gene
				info = {
					"sample" : check.group(2),
					"gene" : check.group(3)
				}
				#find tsv data file for unphased consensus
				unphased_tsv_filepath = target_dir+'/'+check.group(0)+'/'+check.group(1)+'_unphased.fa.annotated_str.tsv'
				if os.path.isfile(unphased_tsv_filepath):
					#find the gene in the database
					cur.execute("SELECT id, name FROM Gene WHERE name LIKE :gene;", info)
					gene_data = cur.fetchone()
					info["gene_num"] = gene_data[0]
					#extract data from tsv file
					tsv_file = open(unphased_tsv_filepath)
					read_tsv = csv.reader(tsv_file, delimiter='\t')
					for row in read_tsv:
						if(row[0] == 'motif'):
							continue
						#get found motif and count
						info["pattern"] = row[0].split('::')[0] #do we only want the forward? or should we get both?
						info["reverse_pattern"] = row[0].split('::')[1]
						info["count"] = int(row[2])
						#check if motif is in database
						cur.execute("SELECT id, pathMin, pathMax from Motif WHERE gene =:gene_num AND pattern LIKE :pattern", info)
						match_motifs = cur.fetchall()
						found = 0
						if match_motifs:
							found = 1
							add_info(cur,info,match_motifs)
							print("Matching")
							print(match_motifs)
						#check if reverse motif is in database
						cur.execute("SELECT id, pathMin, pathMax from Motif WHERE gene =:gene_num AND pattern LIKE :reverse_pattern", info)
						match_motifs = cur.fetchall()
						if match_motifs:
							found = 1
							info["pattern"] = info["reverse_pattern"]
							add_info(cur,info,match_motifs)
							print("Matching reverse")
							print(match_motifs)
						if found == 0:
							info["cat"] = "Unknown"
							cur.execute("INSERT INTO Result (sample,gene,pattern,numRepeats,category) VALUES (:sample,:gene_num,:pattern,:count,:cat)", info)
					tsv_file.close()
					conn.commit()
					genes.append(info["gene"])
		return info["sample"], genes

def add_info(cur,info, match_motifs):
	for motif in match_motifs:
		count = info["count"]
		info["motif_id"] = motif[0]
		pathMin = int(motif[1])
		pathMax = int(motif[2])
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
		cur.execute("INSERT INTO Result (sample,gene,pattern,motif,numRepeats,category) VALUES (:sample,:gene_num,:pattern,:motif_id,:count,:cat)", info)

def get_sample_results(sample, conn):
	cur = conn.cursor()
	cur.execute("SELECT g.name, r.pattern, m.pathMin, m.pathMax, r.numRepeats, r.category FROM Result r LEFT JOIN Gene g ON r.gene = g.id LEFT JOIN Motif m ON r.motif = m.id WHERE r.sample LIKE ?", (sample,))
	results = cur.fetchall()
	return results


def main():
	conn = create_connection('./database.db')
	#target_dir = '/home/alyne/Documents/Thesis/examples_for_alyne/MBXM037256'
	#extract_data(target_dir, conn)
	get_sample_results('MBXM037256', conn)
	conn.close()

if __name__ == '__main__':
	main()
