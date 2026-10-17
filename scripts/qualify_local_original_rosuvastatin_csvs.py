"""Run private local original rosuvastatin source qualification with no data uploads.

No original CSV/person rows printed, logged, committed or read from GitHub.
Study Silveira et al 2021, DOI 10.3389/fphar.2021.752960.
"""
from __future__ import annotations
import argparse
import json
from omics2pk.rosuvastatin_original_local import (
    SourceGateError, summarize_local_sources, qualify_exact_subject_link,
)

def main():
    p=argparse.ArgumentParser(description="Local original rosuvastatin MTBLS3128 clinical omics-PK source gate")
    p.add_argument("--source-file",action="append",default=[],help="Original original-author processed CSV (repeat for all five)")
    p.add_argument("--omics",help="Original ONE-ROW-PER-HUMAN predose feature CSV for one LC-MS modality")
    p.add_argument("--pk",help="Original author quantitative 40-human PK outcome CSV")
    p.add_argument("--omics-id-column",help="EXACT original omics person/sample ID column")
    p.add_argument("--pk-id-column",help="EXACT original PK biological-person ID column")
    p.add_argument("--auc-column",help="EXACT original AUC column, quantified units must be checked manually")
    p.add_argument("--cmax-column",help="EXACT original Cmax column, quantified units must be checked manually")
    p.add_argument("--original-crosswalk",help="OPTIONAL author-authenticated one-to-one sample↔person mapping CSV")
    p.add_argument("--crosswalk-sample-column",default="omics_sample_id")
    p.add_argument("--crosswalk-person-column",default="pk_subject_id")
    p.add_argument("--expected-original-subjects",type=int,default=40)
    args=p.parse_args()
    try:
        result={"local_upload_inventory": summarize_local_sources(args.source_file)} if args.source_file else {}
        if args.omics or args.pk:
            if not all((args.omics,args.pk,args.omics_id_column,args.pk_id_column,
                        args.auc_column,args.cmax_column)):
                p.error("Strict structural join requires explicit --omics --pk and all four exact original column names")
            result["pair_qualification"]=qualify_exact_subject_link(
                args.omics,args.pk,omics_id_col=args.omics_id_column,
                pk_id_col=args.pk_id_column,auc_col=args.auc_column,cmax_col=args.cmax_column,
                expected_subjects=args.expected_original_subjects,
                author_crosswalk_path=args.original_crosswalk,
                crosswalk_sample_col=args.crosswalk_sample_column,
                crosswalk_person_col=args.crosswalk_person_column,
            )
        if not result:
            p.error("Supply at least one --source-file or an explicit --omics and --pk pair")
    except (SourceGateError,OSError,UnicodeError,ValueError) as error:
        p.exit(2,"ABSTAIN: "+str(error)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=="__main__":
    main()
