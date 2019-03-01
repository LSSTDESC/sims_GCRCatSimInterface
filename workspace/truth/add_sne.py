from ExtractSNeParams import add_unsprinkled_sne_params
import os
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--param_db', type=str, default=None,
                        help='Path to the sqlite database to which you '
                        'wish to add unsprinkled supernova parameters')

    args = parser.parse_args()

    if args.param_db is None:
        raise RuntimeError("Must specify param_db")

    add_unsprinkled_sne_params(args.param_db)
