"""Print compressed, nonidentified derived PUBLIC ARIANA 47x assay table.

This helper emits no original IDs and no drug-derived outcome leakage columns.
Original source independently published by the ARIANA authors on GitHub.
Only for reproducible private-workspace packaging, not patient advice.
"""
import base64,csv,gzip,io,json
from scripts.audit_ariana_original_source import (_download,_read_shared,_get_workbook_sheets,_read_sheet)
from hashlib import sha256
from zipfile import ZipFile
ORIGINAL_SHA="7397b13e43cd546b93f83a20b2dc84b6d49e97633173d56711c39973beedfe37"
HEADERS=["allele_1","allele_2","cyp2d6_as","cyp2d6_cnv","metoprolol_ng_ml","daily_dose_mg","ln_oh_solanidine_ratio","ln_ssda_solanidine_ratio","age_years","male_sex","bmi","weight_kg","cyp2d6_inhibitor"]
ORIGINAL_COLS=[1,2,3,4,5,6,9,10,11,12,13,14,15]
def export():
    b=_download()
    if sha256(b).hexdigest()!=ORIGINAL_SHA:raise ValueError("Original author xlsx changed")
    with ZipFile(io.BytesIO(b)) as z:
        sheets=_get_workbook_sheets(z)
        if len(sheets)!=1:raise ValueError("Original source layout changed")
        rows=_read_sheet(z,sheets[0][1],_read_shared(z))
    if len(rows)!=48:raise ValueError("Original 47 original patient rows changed")
    out=io.StringIO()
    writer=csv.writer(out,lineterminator="\n")
    writer.writerow(HEADERS)
    for row in rows[1:]:
        writer.writerow([row.get(i,"") for i in ORIGINAL_COLS])
    contents=out.getvalue().encode("utf-8")
    return {"status":"PUBLIC_SOURCE_ORIGINAL_DEIDENTIFIED_ARIANA47_SELECTED_FIELDS",
      "n_independent_humans":47,"n_columns":len(HEADERS),
      "original_author_xlsx_sha256":ORIGINAL_SHA,
      "derived_csv_sha256":sha256(contents).hexdigest(),
      "excludes_identifiers":True,
      "excludes_drug_derived_metoprolol_ratio":True,
      "gzip_base64_csv":base64.b64encode(gzip.compress(contents,mtime=0)).decode()}
if __name__=="__main__":print(json.dumps(export(),sort_keys=True))
