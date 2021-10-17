DROP TABLE IF EXISTS Gene;
DROP TABLE IF EXISTS Motif;
DROP TABLE IF EXISTS Result;
DROP TABLE IF EXISTS Note;


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

CREATE TABLE Result (
	id INTEGER PRIMARY KEY AUTOINCREMENT,
	sample VARCHAR(20) NOT NULL,
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
	FOREIGN KEY (gene) REFERENCES Gene(id),
	FOREIGN KEY (motif) REFERENCES Motif(id),
	FOREIGN KEY (phaser) REFERENCES Phaser(id)
);

CREATE TABLE Note (
	id INTEGER PRIMARY KEY AUTOINCREMENT,
	result INTEGER,
	status TEXT CHECK(status IN ('In Progress', 'Complete', 'Error')) NOT NULL DEFAULT 'In Progress',
	diagnosis TEXT CHECK(diagnosis IN ('Pathogenic', 'Non-Pathogenic', 'Unknown')) NOT NULL DEFAULT 'Unknown',
	notes TEXT,
	FOREIGN KEY (result) REFERENCES Result(id)
);
