### About

This tool is used to search a proteome data file for upstream regulator matches and output a CSV containing the matches.

### Setup

1) Create your python virtual environment
  ```
  python3.13 -m venv yourvenvname
  ```
2) Activate the venv
  ```
  source ./yourvenvname/bin/activate
  ```
3) Install requirements
  ```
  pip install -r requirements.txt
  ```

### How-To
```
usage: URAProteomeMatcher.py [-h] -pd [PROTEOME_DATA_FILE] -ur [UPSTREAM_REGULATOR_FILE] [-ug [UPSTREAM_REGULATOR_GROUP_FILE]] [-ipn] [-ll [LOG_LEVEL]] [-lf] [-of [OUT_FILE]]

arguments

options:
  -h, --help            show this help message and exit
  -pd, --proteome-data-file [PROTEOME_DATA_FILE]
                        REQUIRED: Local file path to proteome data
  -ur, --upstream-regulator-file [UPSTREAM_REGULATOR_FILE]
                        REQUIRED: Local file path to upstream regulator data
  -ug, --upstream-regulator-group-file [UPSTREAM_REGULATOR_GROUP_FILE]
                        OPTIONAL: Local file path to json gene group file
  -ipn, --ignore-protein-name-matches
                        OPTIONAL: Only match Genes, ignore Protein_Names matches
  -ll, --log-level [LOG_LEVEL]
                        OPTIONAL: Log level, defaults to INFO
  -lf, --log-to-file    OPTIONAL: Log to file, defaults to False
  -of, --out-file [OUT_FILE]
                        OPTIONAL: Local file to output results to
```

Run with proteome_data file & upstream regulator list file:
```
python ./URAProteomeMatcher.py -pd ./test_files/test_proteome.csv -ur ./test_files/test_ur_list.txt
```

Run with a group data file as well:
```
python ./URAProteomeMatcher.py -pd ./test_files/test_proteome.csv -ur ./test_files/test_ur_list.txt -ug ./test_files/test_groups.json
```

Set the log level to debug (shows table data in script run):
```
python ./URAProteomeMatcher.py -pd ./test_files/test_proteome.csv -ur ./test_files/test_ur_list.txt -ll debug
```
