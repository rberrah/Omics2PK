#!/usr/bin/env python3
"""Read only sheet structure of original NCI-60 SWATH supplemental protein workbook.
Original primary source: PMID 31711193? Cite DOI via article instead of assuming PMID.
No source protein intensities, cell-line response rows, or phenotype values printed.
"""
from __future__ import annotations
import io, os, urllib.request, hashlib, zipfile
from openpyxl import load_workbook
URLS=(
    "https://ars.els-cdn.com/content/image/1-s2.0-S2589004219304407-mmc2.xlsx",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC6889472/bin/mmc2.xlsx",
    "https://pmc.ncbi.nlm.nih.gov/articles/instance/6889472/bin/mmc2.xlsx",
    "https://europepmc.org/articles/PMC6889472/bin/mmc2.xlsx",
    "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6889472/supplementaryFiles",
)
def get_source():
    for url in URLS:
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"Omics2PK research source audit (citation: PMC6889472; contact via GitHub)"})
            with urllib.request.urlopen(req,timeout=150) as h:
                b=h.read(180_000_001)
                final=h.geturl()
            print("FETCH",url,"resolved",final,"bytes",len(b),"sha256",hashlib.sha256(b).hexdigest())
            if len(b)>4_000_000 and zipfile.is_zipfile(io.BytesIO(b)):
                if b[:2]==b"PK" and final.lower().endswith(".xlsx"):
                    return b,final
                with zipfile.ZipFile(io.BytesIO(b)) as archive:
                    names=archive.namelist()
                    print("ARCHIVE_ENTRIES", [(name,archive.getinfo(name).file_size) for name in names[:22]])
                    if len(names)==1 and names[0].lower().endswith(".xlsx"):
                        return archive.read(names[0]),final
                    target=[name for name in names if name.lower().endswith(("mmc2.xlsx","mmc2.xls"))]
                    if len(target)==1:
                        return archive.read(target[0]),final+"#"+target[0]
                    if len(b)<40_000_000 and b[:2]==b"PK":
                        return b,final
                print("NO_SINGLE_OFFICIAL_MMC2_XLSX_IN_ARCHIVE")
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
