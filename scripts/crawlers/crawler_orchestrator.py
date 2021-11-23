"""Manages the execution of web crawling."""
import os
import shutil

from scripts.crawlers.dependencies import pubchem
from scripts.crawlers.dependencies import swissadme


def main(file):
    """
    Accesses SwissADME and PubChem to retrieve smile codes and descriptors.
    Uses two separate crawlers for each executed in squence, each producing a separate output.

    :param file: The input file for the crawler, by default data/input_data/chem.csv
    :return None
    """
    prefs = {"download.default_directory": os.path.dirname(os.getcwd()) + "\\cache",
             "directory_upgrade": True}

    shutil.copy(f"{file}", "crawlers/cache/chem.csv")

    for crawler in [pubchem.PubchemCrawler, swissadme.SwissadmeCrawler]:
        crawler.crawl(prefs)


if __name__ == "__main__":
    main("crawlers/cache/chem.csv")
