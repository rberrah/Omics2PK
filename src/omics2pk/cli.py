"""Offline CLI; never uploads study data."""
import argparse
import json
from .exchange import ExchangeError, inspect_exchange, prepare_sweet


def main(argv=None):
    parser = argparse.ArgumentParser(prog="omics2pk")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("validate-exchange").add_argument("exchange_dir")
    sweet = sub.add_parser("prepare-sweet")
    sweet.add_argument("exchange_dir")
    sweet.add_argument("out_dir")
    args = parser.parse_args(argv)
    try:
        out = (inspect_exchange(args.exchange_dir) if args.command == "validate-exchange"
               else prepare_sweet(args.exchange_dir, args.out_dir))
    except (ExchangeError, KeyError, TypeError, ValueError) as exc:
        parser.exit(2, f"INVALID: {exc}\n")
    print(json.dumps(out, indent=2, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
