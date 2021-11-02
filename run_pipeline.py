import cleaner
import merger
import populate_chem
import crawlers.crawler_orchestrator as crawler
import argparse
import shutil

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-inputfile', type=str, help="Input file location.")
    parser.add_argument('-outputfile', type=str, help="Output file destination.")
    args = parser.parse_args()

    crawler.main(args.inputfile if args.inputfile else "./data/input_data/chem.csv")

    if args.outputfile:
        shutil.copy("./data/final/chem.csv", args.outputfile)
