"""Diagnose matrix semantics without recording expression values per donor/cell."""
import csv
import gzip
import hashlib
import io
import json
import math
import re
import tarfile

from scripts.qualify_peng_phh_source import GEO_TAR, load
EXPECTED_GEO_TAR_SHA = "7f114365ba5e977e685f9cee108def64deb0a31c96ea13d95e19c0b42c2b01e8"

def audit():
    binary=load(GEO_TAR,40_000_000)
    if hashlib.sha256(binary).hexdigest()!=EXPECTED_GEO_TAR_SHA:
        raise ValueError("GEO SHA mismatch")
    out=[]
    with tarfile.open(fileobj=io.BytesIO(binary)) as tar:
        for member in tar.getmembers():
            if not member.isfile(): continue
            fp=tar.extractfile(member)
            csv_text=io.TextIOWrapper(gzip.GzipFile(fileobj=fp),encoding="utf-8-sig",newline="")
            reader=csv.reader(csv_text)
            header=next(reader)
            n_features=len(header)-1
            candidate_cyps=("CYP1A2","CYP2B6","CYP2C9","CYP2C19","CYP2D6","CYP3A4")
            listed_cyps=[g for g in candidate_cyps if g in [t.upper() for t in header[1:]]]
            n_rows=0
            counts={"zeros":0,"strict_negative":0,"positive":0,"nonfinite":0,"nonnumeric":0,"noninteger_positive":0}
            gene_header=header[0][:80]
            gene_examples=[]
            offending_classes=[]
            for row in reader:
                n_rows+=1
                if n_rows<=3: gene_examples.append(row[0][:60])
                for v in row[1:]:
                    try: f=float(v)
                    except ValueError:
                        counts["nonnumeric"]+=1
                        continue
                    if not math.isfinite(f):
                        counts["nonfinite"]+=1
                    elif f<0:
                        counts["strict_negative"]+=1
                    elif f==0:
                        counts["zeros"]+=1
                    else:
                        counts["positive"]+=1
                        if abs(f-round(f))>1e-8:
                            counts["noninteger_positive"]+=1
            out.append({
                "filename":member.name.rsplit("/",1)[-1],
                "n_barcodes":n_rows,"n_markers":n_features,
                "source_gene_columns":header[1:],
                "cyp_targets_present":listed_cyps,
                "all_six_cyp_targets_present":len(listed_cyps)==6,
                "first_header":gene_header,"first_three_genes":gene_examples,
                "n_value_categories":counts,
                "count_matrix_all_nonnegative_integer":counts["strict_negative"]==0 and counts["nonfinite"]==0 and counts["nonnumeric"]==0 and counts["noninteger_positive"]==0,
            })
    return {"source": "GSE289636 exact frozen archive", "n_samples":len(out),"sample_structural_semantics":out,
            "note":"No expression values, cell barcodes or donor-level feature vectors released"}
if __name__=="__main__":
    print(json.dumps(audit(),indent=2,sort_keys=True))
