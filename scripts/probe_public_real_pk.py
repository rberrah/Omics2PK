#!/usr/bin/env python3
"""One-time bounded probes of original public article XML access (no sample values)."""
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError
import json
import xml.etree.ElementTree as ET

URLS = {
  "tucatinib_epmc": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13095627/fullTextXML",
  "imatinib_epmc": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12990092/fullTextXML",
  "peng_epmc": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12532551/fullTextXML",
  "tucatinib_pmc_aws": "https://pmc-oa-opendata.s3.amazonaws.com/metadata/PMC13095627.1.json",
  "peng_pmc_aws": "https://pmc-oa-opendata.s3.amazonaws.com/metadata/PMC12532551.1.json",
}
def main():
    result={}
    for name, url in URLS.items():
        try:
            with urlopen(Request(url,headers={"User-Agent":"Omics2PK-public-research/0.1 (no patient data)","Accept":"application/xml, application/json;q=0.9, text/xml;q=0.9"}),timeout=25) as rsp:
                raw = rsp.read(14_000_000)
                info = {"status":rsp.status,"bytes":len(raw),"prefix":raw[:35].decode("utf-8","replace")}
                if name.endswith("epmc") and raw.lstrip().startswith(b"<?xml"):
                    root=ET.fromstring(raw)
                    tables=root.findall(".//table-wrap")
                    info["table_count"]=len(tables)
                    info["tables_with_clearance"]=sum("clearance" in " ".join(t.itertext()).lower() for t in tables)
                    info["tables_with_tucatinib"]=sum("tucatinib" in " ".join(t.itertext()).lower() for t in tables)
                    info["n_table_row_elements"]=sum(len(t.findall(".//tr")) for t in tables)
                if name.endswith("aws") and raw.lstrip().startswith(b"{"):
                    d=json.loads(raw);info["media_count"]=len(d.get("media_urls",[]))
                    info["supplement_docx_count"]=sum(u.endswith(".docx") for u in d.get("media_urls",[]))
                    info["license"]=str(d.get("license",""))[:120]
                result[name]=info
        except Exception as exc:
            result[name]={"error_type":type(exc).__name__,"error_msg":str(exc)[:180]}
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=="__main__":
    main()
