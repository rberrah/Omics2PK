"""Diagnostic safe metadata-only probes for original public GNPS job endpoints."""
import json
from urllib.request import Request, urlopen
from urllib.error import HTTPError,URLError

T="b3966445abad4658a3cdc63c8198a666"
tests=(
 ("gnps_original_download_GET",f"https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task={T}&view=download_cytoscape_data",None),
 ("gnps_original_download_POST",f"https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task={T}&view=download_cytoscape_data",b""),
 ("gnps_original_download_POST_form","https://gnps.ucsd.edu/ProteoSAFe/DownloadResult",f"task={T}&view=download_cytoscape_data".encode()),
 ("gnps_status",f"https://gnps.ucsd.edu/ProteoSAFe/status_json.jsp?task={T}",None),
 ("gnps_download_one_original_file",f"https://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task={T}&block=main&file=qiime2_output/qiime2_table.qza",None),
 ("massive_original_download",f"https://massive.ucsd.edu/ProteoSAFe/DownloadResult?task={T}&view=download_cytoscape_data",None),
)
for name,url,body in tests:
    req=Request(url,data=body,headers={
        "User-Agent":"Omics2PK-research-source-audit/1.0",
        "Range":"bytes=0-31",
        "Accept":"*/*"})
    try:
        with urlopen(req,timeout=18) as r:
            status=r.status
            mime=r.headers.get("Content-Type","")
            allow=r.headers.get("Allow")
            # Bound read to 64 bytes, no contents ever logged
            head=r.read(64)
            result={"test":name,"code":status,"mime":mime,
                    "allow":allow,"first64_read_bytes":len(head)}
    except HTTPError as e:
        result={"test":name,"code":e.code,"allow":e.headers.get("Allow"),
                "mime":e.headers.get("Content-Type")}
    except (URLError,OSError,ValueError,TimeoutError) as e:
        result={"test":name,"error_class":type(e).__name__}
    print(json.dumps(result,sort_keys=True))
