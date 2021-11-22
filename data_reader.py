import os
import re
import csv
import pandas as pd
import sqlite3
from sqlite3 import Error
import json
from collections import OrderedDict


# create a connection to the database
def create_connection(db_file):
    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except Error as e:
        print(e)
    return conn


# extract the result data from the given directory
def extract_data(target_dir, conn):
    if not os.path.isdir(target_dir):
        conn.close()
        exit()
    else:
        with open('config.json', 'r') as f:
            cfg = json.load(f, object_pairs_hook=OrderedDict)
        genes = []
        for dirpath in os.listdir(target_dir):
            # find each gene directory
            check = re.search(r"((\w+)_chr\d+_\d+_\d+_(\w+))\.assembly_dir$",
                              dirpath)
            if check:
                # obtain sample number and gene
                info = {
                    "gene": check.group(3)
                }
                sample = check.group(2)
                # check if sample already exists in the database
                sample_id = get_sample_id(conn, sample)
                if sample_id is None:
                    cur = conn.cursor()
                    cur.execute("INSERT INTO Sample (sample) VALUES (?);",
                                (sample,))
                    info["sample"] = get_sample_id(conn, sample)
                else:
                    info["sample"] = sample_id
                found = False
                for phaser in cfg["phasers"]:
                    info["phaser"] = phaser
                    if phaser in cfg["phased"]:
                        alleles = [1, 2]
                        for allele in alleles:
                            file = phaser + "_hap" + str(allele)
                            info["allele"] = allele
                            found = add_info(target_dir, conn, file, allele,
                                             check, info, cfg)
                    else:
                        allele = 0
                        info["allele"] = allele
                        found = add_info(target_dir, conn, phaser, allele,
                                         check, info, cfg)
                if found:
                    genes.append(info["gene"])
        return sample, genes


# assign the category based on the matching motif's pathogenicity range
def add_motif_match(cur, info, motif):
    count = info["count"]
    info["motif_id"] = motif[0]
    pathMin = int(motif[2])
    pathMax = int(motif[3])
    if pathMin is not None and pathMax is not None:
        if count >= pathMin and count <= pathMax:
            info["cat"] = "Likely Pathogenic"
        elif count < pathMin or count > pathMax:
            info["cat"] = "Likely Non-Pathogenic"
        else:
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
    cur.execute(
        "INSERT INTO Result (sample,gene,allele,phaser,pattern,motif," +
        "numRepeats,category,match,score,sequence) VALUES (:sample," +
        ":gene_num,:allele,:phaser,:found_pattern,:motif_id,:count,:cat," +
        ":match,:score,:found_sequence)", info)


# add the information extracted from the directory
def add_info(target_dir, conn, file, allele, check, info, cfg):
    gene_found = False
    cur = conn.cursor()
    # find tsv data file for unphased consensus
    tsv_filepath = target_dir + '/' + check.group(0) + '/' + check.group(1)
    tsv_filepath = tsv_filepath + '_' + file + '.fa.annotated_str.tsv'

    min_len = cfg["min_pattern_length"]
    max_len = cfg["max_pattern_length"]

    if os.path.isfile(tsv_filepath):
        gene_found = True
        # find the gene in the database
        cur.execute("SELECT id, name FROM Gene WHERE name LIKE :gene;", info)
        gene_data = cur.fetchone()
        info["gene_num"] = gene_data[0]
        # extract data from tsv file
        tsv_file = open(tsv_filepath)
        read_tsv = csv.reader(tsv_file, delimiter='\t')
        empty_file = True
        for row in read_tsv:
            if row[0] == 'motif':
                continue
            # get found pattern, reverse inverse and count
            info["pattern"] = row[0].split('::')[0]
            if len(
              info["pattern"]) < min_len or len(info["pattern"]) > max_len:
                continue
            empty_file = False
            info["rev_inv_pattern"] = row[0].split('::')[1]
            info["count"] = float(row[3])
            info["match"] = int(row[4])
            info["score"] = int(row[6])
            info["sequence"] = row[11].split('::')[0]
            info["rev_inv_sequence"] = row[11].split('::')[1]
            # get all motifs associated with the gene
            cur.execute("SELECT id, pattern, pathMin, pathMax, gene FROM " +
                        "Motif WHERE gene =:gene_num", info)
            motifs = cur.fetchall()
            found = 0
            for motif in motifs:
                # generate all possibile cyclical variations of the motif
                variations = generate_variations(motif[1])
                # find a match
                if info["pattern"] in variations:
                    found = 1
                    info["found_pattern"] = info["pattern"]
                    info["found_sequence"] = info["sequence"]
                    add_motif_match(cur, info, motif)
                if info["rev_inv_pattern"] in variations:
                    found = 1
                    info["found_pattern"] = info["rev_inv_pattern"]
                    info["found_sequence"] = info["rev_inv_sequence"]
                    add_motif_match(cur, info, motif)
            # if no match found, insert separately
            if found == 0:
                info["cat"] = "Unknown"
                cur.execute(
                    "INSERT INTO Result (sample, gene, allele, phaser, " +
                    "pattern, numRepeats, category, match, score, sequence) "
                    "VALUES (:sample,:gene_num,:allele,:phaser,:pattern," +
                    ":count,:cat,:match,:score,:sequence)", info)
        if empty_file is True:
            info["cat"] = "Unknown"
            cur.execute("INSERT INTO Result (sample, gene, allele, phaser, " +
                        "category) VALUES (:sample,:gene_num,:allele," +
                        ":phaser,:cat)", info)

        tsv_file.close()
        conn.commit()
    return gene_found


# checks if the sample should be added
def check_add_sample(conn, target_dir):
    if not os.path.isdir(target_dir):
        conn.close()
        exit()
    else:
        for dirpath in os.listdir(target_dir):
            check = re.search(r"((\w+)_chr\d+_\d+_\d+_(\w+))\.assembly_dir$",
                              dirpath)
            if check:
                sample = check.group(2)
                if get_sample_id(conn, sample) is not None:
                    error = "is already in the database and cannot be readded."
                    return "Sample " + sample + " " + error
                return None
    return "The folder does not appear to be in the right format."


# generates the various configurations of a motif/pattern
# e.g. CAG -> CAG, AGC, GCA
def generate_variations(pattern):
    patterns = []
    length = len(pattern)
    bases = list(pattern)
    i = 0
    while i < length:
        patterns.append("".join(bases[i:]+bases[:i]))
        i = i + 1
    return patterns


# returns the number of samples in the database
def get_sample_count(conn):
    cur = conn.cursor()
    cur.execute("SELECT COUNT(*) FROM Sample;")
    count = cur.fetchone()
    return count


# returns the sample ID in the database from the name of a sample
def get_sample_id(conn, sample):
    cur = conn.cursor()
    cur.execute("SELECT id FROM Sample WHERE sample LIKE ?;", (sample,))
    id = cur.fetchone()
    if id is not None:
        return id[0]
    return id


# returns the gene ID in the database from the name of a gene
def get_gene_id(conn, gene):
    cur = conn.cursor()
    cur.execute("SELECT id FROM Gene WHERE name LIKE ?;", (gene,))
    id = cur.fetchone()
    if id is not None:
        return id[0]
    return id


# returns the results associated with a specific sample
def get_sample_results(sample, conn):
    sample = get_sample_id(conn, sample)
    if sample is None:
        return None
    cur = conn.cursor()
    cur.execute("SELECT r.id, g.name, r.pattern, m.pattern, m.pathMin, " +
                "m.pathMax, r.numRepeats, r.allele, r.phaser, r.category, " +
                "s.diagnosis, s.status FROM Result r " +
                "LEFT JOIN Gene g ON r.gene = g.id " +
                "LEFT JOIN Motif m ON r.motif = m.id " +
                "LEFT JOIN Status s ON s.gene = g.id AND s.sample = r.sample" +
                " WHERE r.sample LIKE ? " +
                "ORDER BY g.name ASC", (sample,))
    results = cur.fetchall()
    return results


# returns the data for a single result
def get_sample_result(result, conn):
    cur = conn.cursor()
    cur.execute("SELECT s.sample, g.name, r.pattern, r.phaser, r.allele, " +
                "r.numRepeats, m.pattern, m.pathMin, m.pathMax, r.category, " +
                "r.match, r.score, r.sequence FROM Result r " +
                "LEFT JOIN Gene g ON r.gene = g.id " +
                "LEFT JOIN Motif m ON r.motif = m.id LEFT JOIN Sample s " +
                "ON r.sample = s.id WHERE r.id LIKE ?", (result,))
    result = cur.fetchone()
    return result


# returns all results associated with the given gene and sample
def get_sample_gene_results(sample, gene, conn):
    sample = get_sample_id(conn, sample)
    if sample is None:
        return None
    cur = conn.cursor()
    cur.execute("SELECT r.id, r.pattern, m.pattern, r.numRepeats, m.pathMin," +
                " m.pathMax, r.match, r.allele, r.phaser, r.category " +
                "FROM Result r LEFT JOIN Gene g " +
                "ON r.gene = g.id LEFT JOIN Motif m ON " +
                "r.motif = m.id WHERE r.sample LIKE ? AND g.name LIKE ?",
                (sample, gene,))
    results = cur.fetchall()
    return results


# deletes the sample from the database
def delete_sample(conn, sample):
    sample = get_sample_id(conn, sample)
    if sample is None:
        return False
    cur = conn.cursor()
    cur.execute("DELETE FROM Result WHERE sample LIKE ?", (sample,))
    cur.execute("DELETE FROM Sample WHERE id LIKE ?", (sample,))
    conn.commit()
    return True


# returns the data for a specific gene
def get_gene_info(gene, conn):
    cur = conn.cursor()
    cur.execute("SELECT g.name, g.chromosome, g.phenotype," +
                "g.inheritanceMode, g.startCoordHG38, g.endCoordHG38 " +
                "FROM Gene g WHERE g.name LIKE ?", (gene,))
    gene_info = cur.fetchone()
    return gene_info


# returns the sample names that start with the given query
def find_samples(sample, conn):
    cur = conn.cursor()
    cur.execute("SELECT sample FROM Sample WHERE sample LIKE ?",
                ('{}%'.format(sample),))
    results = cur.fetchall()
    return results


# returns the gene names that start with the given query
def find_gene(gene, conn):
    cur = conn.cursor()
    cur.execute("SELECT name FROM Gene WHERE name LIKE ?",
                ('{}%'.format(gene),))
    results = cur.fetchall()
    return results


# returns the notes, status and diagnosis for a sample's gene
def get_info(conn, sample, gene):
    sample_id = get_sample_id(conn, sample)
    if sample_id is None:
        return None
    gene_id = get_gene_id(conn, gene)
    if gene_id is None:
        return None
    cur = conn.cursor()
    cur.execute("SELECT id, note FROM Note WHERE gene LIKE ? AND " +
                "sample LIKE ?", (gene_id, sample_id,))
    notes = cur.fetchall()
    cur.execute("SELECT status, diagnosis FROM Status WHERE gene LIKE ? AND " +
                "sample LIKE ?", (gene_id, sample_id,))
    info = cur.fetchone()
    status = None
    diagnosis = None
    if info is not None:
        status = info[0]
        diagnosis = info[1]
    return status, diagnosis, notes


# updates the status of a sample's gene
def update_status(conn, sample, gene, status, diagnosis):
    sample_id = get_sample_id(conn, sample)
    if sample_id is None:
        return None
    gene_id = get_gene_id(conn, gene)
    if gene_id is None:
        return None
    cur = conn.cursor()
    cur.execute("DELETE FROM Status WHERE gene LIKE ? AND sample LIKE ?",
                (gene_id, sample_id,))
    # make sure that succeeded
    cur.execute("INSERT INTO Status (sample, gene, status, diagnosis) VALUES" +
                " (?, ?, ?, ?)", (sample_id, gene_id, status, diagnosis))
    conn.commit()


# adds a note to the sample's gene
def add_note(conn, sample, gene, note):
    sample_id = get_sample_id(conn, sample)
    if sample_id is None:
        return None
    gene_id = get_gene_id(conn, gene)
    if gene_id is None:
        return None
    cur = conn.cursor()
    cur.execute("INSERT INTO Note (sample, gene, note) VALUES (?, ?, ?)",
                (sample_id, gene_id, note))
    conn.commit()


# deletes the specified note from the sample's gene
def delete_note(conn, note):
    cur = conn.cursor()
    cur.execute("DELETE FROM Note WHERE id LIKE ?", (note,))
    conn.commit()


# saves all notes in the database to a tsv file
def notes_to_tsv(conn):
    db_df = pd.read_sql_query(
            "SELECT s.sample, g.name, n.note FROM Note n " +
            "LEFT JOIN Sample s on s.id = n.sample " +
            "LEFT JOIN Gene g on g.id = n.gene", conn)
    db_df.to_csv('notes.tsv', sep="\t", index=False)


# saves all assigned statuses and diagnoses to a tsv file
def statuses_to_tsv(conn):
    db_df = pd.read_sql_query(
            "SELECT sp.sample, g.name, s.status, s.diagnosis FROM Status s " +
            "LEFT JOIN Sample sp on sp.id = s.sample " +
            "LEFT JOIN Gene g on g.id = s.gene", conn)
    db_df.to_csv('statuses.tsv', sep="\t", index=False)


# def main():
#     conn = create_connection('./database.db')
#     conn.close()
#
#
# if __name__ == '__main__':
#     main()
