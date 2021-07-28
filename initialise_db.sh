cat ./schema.sql | sqlite3 database.db;
python3 populate_db.py;
