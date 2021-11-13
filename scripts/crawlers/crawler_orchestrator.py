import os
import shutil
from scripts.crawlers.dependencies import pubchem
from scripts.crawlers.dependencies import swissadme


def main(file):
    prefs = {'download.default_directory': os.path.dirname(os.getcwd()) + f'\\cache', 'directory_upgrade': True}

    shutil.copy(f'{file}', 'crawlers/cache/chem.csv')

    for crawl in [pubchem.crawl, swissadme.crawl]:
        crawl(prefs)


if __name__ == '__main__':
    main("crawlers/cache/chem.csv")
