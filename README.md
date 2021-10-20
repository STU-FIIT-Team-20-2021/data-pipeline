# Usage:

**Data folder must contain:**
- merged.csv (up to date data)
- new_pubchem.csv (new pubchem scrape)
- new_swiss.csv (new swissadme scrape)

The new_pubchem and new_swiss csv files must have the same ordering of new compounds.
Once files are in the right place, just run **merger.py** to update merged.csv

### populate_chem.py

Populates the new_merged.csv with chemical information using rdkit (cycles, charges, etc.)