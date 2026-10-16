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
from scripts.validate_peng_rna_activity import EXPECTED_GEO_TAR_SHA

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
            n_cells=len(header)-1
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
                "n_cells":n_cells,"n_rows":n_rows,
                "first_header":gene_header,"first_three_genes":gene_examples,
                "n_value_categories":counts,
                "count_matrix_all_nonnegative_integer":counts["strict_negative"]==0 and counts["nonfinite"]==0 and counts["nonnumeric"]==0 and counts["noninteger_positive"]==0,
            })
    return {"source": "GSE289636 exact frozen archive", "n_samples":len(out),"sample_structural_semantics":out,
            "note":"No expression values, cell barcodes or donor-level feature vectors released"}
if __name__=="__main__":
    print(json.dumps(audit(),indent=2,sort_keys=True))
