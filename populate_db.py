import sqlite3
from sqlite3 import Error
import csv
from data_reader import get_gene_id


def create_connection(db_file):
    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except Error as e:
        print(e)
    return conn


def main():
    conn = create_connection('./database.db')
    cur = conn.cursor()
    cur.execute("DELETE FROM Gene;")
    cur.execute("DELETE FROM Motif;")

    genes = open("data/genes.tsv")
    read_genes = csv.reader(genes, delimiter="\t")
    for row in read_genes:
        cur.execute("INSERT INTO Gene (name, chromosome, phenotype, inheritanceMode, startCoordHG38, endCoordHG38) VALUES (?, ?, ?, ?, ?, ?)", (row[0], row[1], row[2], row[3], row[4], row[5]))
    genes.close()

    motifs = open("data/motifs.tsv")
    read_motifs = csv.reader(motifs, delimiter="\t")
    for row in read_motifs:
        gene_id = get_gene_id(conn, row[0])
        if gene_id is not None:
            cur.execute("INSERT INTO Motif (gene, pattern, size, pathMin, pathMax) VALUES (?, ?, ?, ?, ?)", (gene_id, row[1], row[2], row[3], row[4]))
    motifs.close()
    # put in a try catch?
    conn.commit()
    conn.close()


if __name__ == '__main__':
    main()
