#!/usr/bin/env python3
from __future__ import annotations
import argparse, warnings
from pathlib import Path
from Bio import SeqIO, BiopythonParserWarning

def parse_summary(path):
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        records=list(SeqIO.parse(str(path),"genbank"))
    return {
        "records":len(records),
        "bp":sum(len(r.seq) for r in records),
        "warnings":[str(w.message) for w in ws if issubclass(w.category,BiopythonParserWarning)]
    }

def normalize_locus_line(line):
    raw=line.rstrip("\r\n")
    if not raw.startswith("LOCUS"):
        return line,None

    parts=raw.split()
    # Already normal enough: LOCUS <name> <length> bp ...
    if len(parts)>=4 and parts[2].isdigit() and parts[3].lower()=="bp":
        return line,None

    # Known Prokka/SPAdes glued layout:
    # LOCUS NODE_1_length_467231_cov_77.474813467231 bp DNA linear
    if len(parts)>=3 and parts[0]=="LOCUS" and parts[2].lower()=="bp":
        glued=parts[1]
        marker="_length_"
        covmarker="_cov_"
        if marker in glued and covmarker in glued:
            left,after_len=glued.split(marker,1)
            emb_len,after_cov=after_len.split(covmarker,1)
            if emb_len.isdigit() and after_cov.endswith(emb_len):
                cov=after_cov[:-len(emb_len)]
                if cov and all(ch.isdigit() or ch=="." for ch in cov):
                    name=f"{left}{marker}{emb_len}{covmarker}{cov}"
                    rest=" ".join(parts[2:])
                    return f"LOCUS       {name} {emb_len} {rest}\n","SPADES_GLUE"

    return line,None

def validate_or_normalize(src,out,audit):
    src=Path(src); out=Path(out); audit=Path(audit)
    if not src.is_file():
        raise SystemExit(f"Input GenBank not found: {src}")

    original_error=""
    try:
        before=parse_summary(src)
        audit.write_text(
            "GENBANK INPUT AUDIT v0.8.2\n"
            f"input={src}\nstatus=VALID_AS_IS\n"
            f"records={before['records']}\ntotal_bp={before['bp']}\n"
            f"parser_warnings={len(before['warnings'])}\noriginal_modified=NO\n",
            encoding="utf-8"
        )
        return src,before,0
    except Exception as e:
        original_error=repr(e)

    fixed=0
    with open(src,encoding="utf-8",errors="replace") as fi, \
         open(out,"w",encoding="utf-8",newline="\n") as fo:
        for line in fi:
            nl,kind=normalize_locus_line(line)
            if kind: fixed+=1
            fo.write(nl)

    if fixed==0:
        raise SystemExit(
            "GenBank parsing failed and no recognized safe LOCUS repair was found.\n"
            f"Original parse error: {original_error}"
        )

    after=parse_summary(out)
    raw=src.read_text(encoding="utf-8",errors="replace")
    terminators=sum(1 for x in raw.splitlines() if x.strip()=="//")
    if terminators and after["records"]!=terminators:
        raise SystemExit(
            f"Normalization validation failed: parsed records={after['records']} "
            f"but flatfile records={terminators}"
        )

    audit.write_text(
        "GENBANK INPUT AUDIT v0.8.2\n"
        f"input={src}\nstatus=NORMALIZED_INTERNAL_COPY\n"
        f"normalized_copy={out}\nrecords={after['records']}\n"
        f"total_bp={after['bp']}\nlocus_lines_fixed={fixed}\n"
        f"parser_warnings_after={len(after['warnings'])}\n"
        f"original_parse_error={original_error}\noriginal_modified=NO\n",
        encoding="utf-8"
    )
    return out,after,fixed

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("genbank")
    ap.add_argument("--out",required=True)
    ap.add_argument("--audit",required=True)
    a=ap.parse_args()
    used,summary,n=validate_or_normalize(a.genbank,a.out,a.audit)
    print("GenBank input validation: PASS")
    print("Input used:",used)
    print("Records:",summary["records"])
    print("Total bp:",summary["bp"])
    print("LOCUS lines normalized:",n)
    print("Original input modified: NO")

if __name__=="__main__":
    main()
