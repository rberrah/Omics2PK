#!/usr/bin/env python3
"""Read only sheet structure of original NCI-60 SWATH supplemental protein workbook.
Original primary source: PMID 31711193? Cite DOI via article instead of assuming PMID.
No source protein intensities, cell-line response rows, or phenotype values printed.
"""
from __future__ import annotations
import io, os, urllib.request, hashlib, zipfile
from openpyxl import load_workbook
URLS=(
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC6889472/bin/mmc2.xlsx",
    "https://pmc.ncbi.nlm.nih.gov/articles/instance/6889472/bin/mmc2.xlsx",
    "https://europepmc.org/articles/PMC6889472/bin/mmc2.xlsx",
    "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6889472/supplementaryFiles",
)
def get_source():
    for url in URLS:
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"Omics2PK research source audit (citation: PMC6889472; contact via GitHub)"})
            with urllib.request.urlopen(req,timeout=55) as h:
                b=h.read(50_000_001)
                final=h.geturl()
            print("FETCH",url,"resolved",final,"bytes",len(b),"sha256",hashlib.sha256(b).hexdigest())
            if len(b)>4_000_000 and zipfile.is_zipfile(io.BytesIO(b)):
                return b,final
        except Exception as exc:
            print("ABSTAIN_CANDIDATE_UNAVAILABLE",url,type(exc).__name__,str(exc)[:120])
    raise SystemExit("ABSTAIN_NO_OFFICIAL_PUBLISHED_PROTEOME_XLSX_AVAILABLE")
def main():
    blob,url=get_source()
    wb=load_workbook(io.BytesIO(blob),read_only=True,data_only=True)
    for ws in wb.worksheets:
        print("SHEET",repr(ws.title),"rows",ws.max_row,"cols",ws.max_column)
        for idx,row in enumerate(ws.iter_rows(min_row=1,max_row=3,values_only=True),start=1):
            print("HEADER_ROW",ws.title,idx,repr([str(x)[:54] if x is not None else "" for x in row[:18]]))
    wb.close()
    print("SOURCE_STRUCTURAL_AUDIT_ONLY",url)
if __name__=="__main__":main()
