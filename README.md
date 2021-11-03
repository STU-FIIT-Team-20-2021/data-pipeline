# Usage:

merger.py automatically attempts to read **'new_pubchem.csv'** and **'new_swiss.csv'** from **./data** and then outputs **./merger/merge.csv**

Altnernatively, use **merger.py -h** to show arguments that allow passing of specific files on input.

---
populate_chem.py attempts to read data from **./merger/merge.csv** then output it to **./chem/new_attrib.csv**