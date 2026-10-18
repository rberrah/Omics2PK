"""Run Omics2PK on authorized ORIGINAL HLC/Synapse ZIPs, locally and privately.

NO original human records are uploaded, committed, or output. For research
use only under original participant privacy and Synapse sharing terms.

Example from a local/research-workspace original dataset directory:
  python -m scripts.benchmark_authorized_hlc_original \
     --expression /private/curatedExpressionLiver.zip \
     --phenotype /private/curatedPhenotypesynapse.zip \
     --genotype /private/curatedGenotype.zip --fullarray

No files are saved by default; stdout is privacy-safe aggregate scores and
original study provenance only. User can use this in a hosted research
workspace too; never supply original source or tokens to public CI.
"""
from __future__ import annotations

import argparse
import json

from omics2pk.hlc_real_cyp_validation import (
    OriginalHLCQualificationError,benchmark_original_hlc_expression
)


def main(argv=None):
    parser=argparse.ArgumentParser(
        description="Private original HLC 10 human microsomal CYP activities: target RNA and multi-gene controls")
    parser.add_argument("--expression",required=True,
                        help="AUTHORIZED local original Synapse expression ZIP, never GitHub")
    parser.add_argument("--phenotype",required=True,
                        help="AUTHORIZED local original Synapse phenotype ZIP, never GitHub")
    parser.add_argument("--genotype",default=None,
                        help="Optional authorized Synapse original genotype ZIP")
    parser.add_argument("--fullarray",action="store_true",
                        help="Expensive optional 40k-probe PCA8 train-only control, 5 CYP endpoints")
    parser.add_argument("--permutations",type=int,default=199,
                        help="Within age-sex stratum RNA negative-control permutations (default 199)")
    args=parser.parse_args(argv)
    if args.permutations<1 or args.permutations>999:
        raise SystemExit("Permutations must be between 1 and 999")
    try:
        results=benchmark_original_hlc_expression(
            args.expression,args.phenotype,permutations=args.permutations)
        if args.genotype is not None:
            from omics2pk.hlc_real_genomics_validation import (
                benchmark_original_hlc_rna_plus_pgx)
            results["germline_genotype_subset"]=benchmark_original_hlc_rna_plus_pgx(
                args.expression,args.phenotype,args.genotype)
        if args.fullarray:
            from omics2pk.hlc_real_fullarray_validation import (
                benchmark_original_hlc_fullarray)
            results["full_transcriptome_PCA8"]=benchmark_original_hlc_fullarray(
                args.expression,args.phenotype)
    except (OriginalHLCQualificationError,OSError,ValueError) as err:
        raise SystemExit("HLC original source scientific qualification denied: "
                         +str(err)) from None
    # Never print file paths, participant IDs, genotype calls or individual
    # original CYP assays, only aggregate predictions/errors.
    print(json.dumps(results,indent=2,sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
