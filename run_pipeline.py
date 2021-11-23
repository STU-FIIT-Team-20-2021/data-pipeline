from scripts import populate_chem, cleaner, merger
import scripts.crawlers.crawler_orchestrator as crawler
import argparse
import shutil


def run_pipe():
    """
    Runs the pipeline in sequence of steps: Crawler, merger, chem populator, cleaner.

    The chem populator also merges existing data with new data before it's sent to the cleaner.

    Multiple arguments are supported, as described in README.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-inputfile', type=str, help="Input file location.")
    parser.add_argument('-outputfile', type=str, help="Output file destination.")
    parser.add_argument('-noscrape', action='store_false', help="Skip scraping/crawler step.")
    args = parser.parse_args()

    if not args.noscrape:
        crawler.main(args.inputfile if args.inputfile else "./data/input_data/input.csv")
    merger.merge()
    populate_chem.populate()
    cleaner.main()

    if args.outputfile:
        shutil.copy("./data/final/phototox.csv", args.outputfile)


if __name__ == "__main__":
    run_pipe()
