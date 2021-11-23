"""Manages the execution of the pipeline."""
import argparse
import shutil
import sys

from scripts import populate_chem, cleaner, merger
import scripts.crawlers.crawler_orchestrator as crawler


def run_pipe():
    """
    Runs the pipeline in sequence of steps: Crawler, merger, chem populator, cleaner.

    The chem populator also merges existing data with new data before it's sent to the cleaner.

    Multiple arguments are supported, as described in README.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-inputfile", type=str, help="Input file location.")
    parser.add_argument("-outputfile", type=str, help="Output file destination.")
    parser.add_argument("-noscrape", action="store_false", help="Skip scraping/crawler step.")
    parser.add_argument("-executestep", type=int,
                        help="Execute only one step of the pipeline by index (from 0): "
                             "[crawler, merger, chem_populator, cleaner]")
    parser.add_argument("-args", type=str,
                        help="Pass custom positional arguments to function. Used with -executestep")
    parser.add_argument("-kwargs", type=str,
                        help="Pass custom keyword arguments to function. Used with -executestep")

    args = parser.parse_args()
    steps = [crawler.main, merger.merge, populate_chem.populate, cleaner.main]

    posargs = ()
    kwargs = {}

    if args.args:
        posargs = args.args.split(',')

    if args.kwargs:
        keys = args.kwargs.split(',')
        keys = map(lambda x: x.split('='), keys)
        for keypair in keys:
            kwargs[keypair[0]] = keypair[1]

    if args.executestep:
        if args.args:
            steps[args.executestep](*posargs)
            sys.exit(0)
        elif args.kwargs:
            steps[args.executestep](**kwargs)
            sys.exit(0)
        else:
            steps[args.executestep]()
            sys.exit(0)

    if not args.noscrape:
        crawler.main(args.inputfile if args.inputfile else "./data/input_data/input.csv")
    merger.merge()
    populate_chem.populate()
    cleaner.main()

    if args.outputfile:
        shutil.copy("./data/final/phototox.csv", args.outputfile)


if __name__ == "__main__":
    run_pipe()
