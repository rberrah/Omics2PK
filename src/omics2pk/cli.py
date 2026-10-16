"""Offline CLI; never uploads study data."""
import argparse
import json
from .exchange import ExchangeError, inspect_exchange, prepare_sweet
from .sweet_features import SWEETOutputError, summarize_networks


def main(argv=None):
    parser = argparse.ArgumentParser(prog="omics2pk")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("validate-exchange").add_argument("exchange_dir")
    sweet = sub.add_parser("prepare-sweet")
    sweet.add_argument("exchange_dir")
    sweet.add_argument("out_dir")
    summary = sub.add_parser("summarize-sweet")
    summary.add_argument("edges_dir")
    summary.add_argument("patient_list")
    summary.add_argument("out_dir")
    args = parser.parse_args(argv)
    try:
        if args.command == "validate-exchange":
            out = inspect_exchange(args.exchange_dir)
        elif args.command == "prepare-sweet":
            out = prepare_sweet(args.exchange_dir, args.out_dir)
        else:
            out = summarize_networks(args.edges_dir, args.patient_list, args.out_dir)
    except (ExchangeError, SWEETOutputError, KeyError, TypeError, ValueError) as exc:
        parser.exit(2, f"INVALID: {exc}\n")
    print(json.dumps(out, indent=2, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
