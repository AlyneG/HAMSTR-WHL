DROP TABLE IF EXISTS Gene;
DROP TABLE IF EXISTS Motif;
DROP TABLE IF EXISTS Phaser;
DROP TABLE IF EXISTS Sample;
DROP TABLE IF EXISTS Result;
DROP TABLE IF EXISTS Note;
DROP TABLE IF EXISTS Status;


CREATE TABLE Gene (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  name VARCHAR(10) UNIQUE NOT NULL,
  chromosome VARCHAR(10) NOT NULL,
  phenotype TEXT,
  inheritanceMode VARCHAR(2) CHECK(inheritanceMode IN ('AD', 'AR', 'XL', 'CMPLX')) NULL DEFAULT NULL,
  startCoordHG38 VARCHAR(12),
  endCoordHG38 VARCHAR(12)
);

CREATE TABLE Motif (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  gene INTEGER NOT NULL,
  pattern VARCHAR(10) NOT NULL,
  size INTEGER,
  pathMin INTEGER,
  pathMax INTEGER,
  FOREIGN KEY (gene) REFERENCES Gene(id)
);

CREATE TABLE Phaser (
	id INTEGER PRIMARY KEY,
	name TEXT
);

CREATE TABLE Sample (
	id INTEGER PRIMARY KEY,
	sample VARCHAR(20)
);

CREATE TABLE Result (
	id INTEGER PRIMARY KEY AUTOINCREMENT,
	sample INTEGER NOT NULL,
	gene INTEGER NOT NULL,
	allele INTEGER,
	phaser INTEGER NOT NULL,
	pattern TEXT,
	motif INTEGER,
	numRepeats REAL,
	match INTEGER,
	score INTEGER,
	category TEXT CHECK(category IN ('Likely Pathogenic', 'Likely Non-Pathogenic', 'Unknown')) NULL DEFAULT NULL,
	sequence TEXT,
	FOREIGN KEY (sample) REFERENCES Sample(id),
	FOREIGN KEY (gene) REFERENCES Gene(id),
	FOREIGN KEY (motif) REFERENCES Motif(id),
	FOREIGN KEY (phaser) REFERENCES Phaser(id)
);

CREATE TABLE Note (
	id INTEGER PRIMARY KEY AUTOINCREMENT,
	sample INTEGER,
	gene INTEGER,
	note TEXT,
	FOREIGN KEY (sample) REFERENCES Sample(id),
	FOREIGN KEY (gene) REFERENCES Gene(id)
);

CREATE TABLE Status (
	id INTEGER PRIMARY KEY AUTOINCREMENT,
	sample INTEGER,
	gene INTEGER,
	status TEXT CHECK(status IN ('In Progress', 'Complete', 'Error')),
	diagnosis TEXT CHECK(diagnosis IN ('Pathogenic', 'Non-Pathogenic', 'Unknown')),
	FOREIGN KEY (sample) REFERENCES Sample(id),
	FOREIGN KEY (gene) REFERENCES Gene(id),
	CONSTRAINT UC_Status UNIQUE (sample,gene)
);
