"""CLI interface for monviso_reloaded project.

Be creative! do whatever you want!

- Install click or typer and create a CLI app
- Use builtin argparse
- Start a web application
- Import things from your .base module
"""


# def main():  # pragma: no cover
#     """
#     The main function executes on commands:
#     `python -m monviso_reloaded` and `$ monviso_reloaded `.

#     This is your program's entry point.

#     You can change this function to do whatever you want.
#     Examples:
#         * Run a test suite
#         * Run a server
#         * Do some other stuff
#         * Run a command line application (Click, Typer, ArgParse)
#         * List all available tasks
#         * Run an application (Flask, FastAPI, Django, etc.)
#     """
#     print("This will do something")

from monviso_reloaded.utils.utils import parse_input
from monviso_reloaded.utils.utils import get_parameters
from monviso_reloaded.utils.utils import make_gene_directories
from monviso_reloaded.utils.utils import add_arguments
from monviso_reloaded.utils.utils import check_arguments
from monviso_reloaded.utils.utils import print_parameters

import argparse
import sys



# mut_file = "/mnt/d/work/monviso/mutations.txt"
# param_file = "/mnt/d/work/monviso/parameters.dat"
#output_dir = "/mnt/d/work/monviso/parameters.dat"

def main(argv=None) -> None:
    """
    Main function

    :param argv: argv

    :return: None
    """
    parser = argparse.ArgumentParser()
    add_arguments(parser)
    args, unparsed = parser.parse_known_args(argv)
    check_arguments(args)
    parameters = get_parameters(args.par_file) if args.par_file else None
    print_parameters(args, parameters)
    # blocks, protein_list = parse_input(mut_file)

    # print(parameters)

    # make_gene_directories(blocks, parameters["OUT_DIR"])

def init() -> None:
    if __name__ == "__main__":
        sys.exit(main())


init()
