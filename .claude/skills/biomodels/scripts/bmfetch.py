#!/usr/bin/env python3
"""bmfetch - thin CLI over the BioModels REST API (https://www.ebi.ac.uk/biomodels).

Every endpoint used here was verified live (see ../SKILL.md for the probe log).
Stdlib only, so it runs under the repo venv or bare python3.

Subcommands
-----------
  search   <query> [-n N] [--facets]     keyword / field search, prints id | name
  info     <ID>                          metadata (name, publication, files, curation)
  files    <ID>                          file manifest with sizes and sha256
  get      <ID> [-o DIR] [--name FILE]   download the main SBML (or a named file)
  probe    <UNIPROT[,UNIPROT...]>        which models annotate these proteins
"""
import argparse, hashlib, json, os, sys, time, urllib.error, urllib.parse, urllib.request

BASE = "https://www.ebi.ac.uk/biomodels"
UA = {"User-Agent": "bmfetch/1.0 (+Ezh2_CyclinD1_Sphase)", "Accept": "application/json"}


def _get(path, params=None, raw=False, retries=3, timeout=90):
    """GET with redirect following, JSON decode, and retry on transient failure.

    NOTE: the API 301-redirects http->https and MODEL* -> BIOMD* accessions, so
    redirects must be followed (urllib does this by default; curl needs -L).
    """
    url = f"{BASE}{path}"
    if params:
        url += "?" + urllib.parse.urlencode(params)
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
                body = r.read()
            return body if raw else json.loads(body)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            last = e
            # search is slow and intermittently times out; back off and retry
            time.sleep(2 * (attempt + 1))
    raise SystemExit(f"bmfetch: {url} failed after {retries} tries: {last}")


def cmd_search(a):
    d = _get("/search", {"query": a.query, "numResults": a.n, "format": "json"})
    print(f"# matches: {d.get('matches')}  (showing {len(d.get('models', []))})")
    for m in d.get("models", []):
        print(f"{m['id']}\t{m.get('format','?'):5s}\t{m.get('name','')}")
    if a.facets:
        for f in d.get("facets", []):
            vals = [f"{v.get('value')}({v.get('count')})" for v in f.get("facetValues", [])][:10]
            print(f"# facet {f.get('id')}: {', '.join(vals)}")


def cmd_info(a):
    d = _get(f"/{a.id}", {"format": "json"})
    pub = d.get("publication", {}) or {}
    print(f"id         : {d.get('submissionId')} / {d.get('publicationId')}")
    print(f"name       : {d.get('name')}")
    print(f"format     : {(d.get('format') or {}).get('name')}")
    print(f"curation   : {d.get('curationStatus')}")
    appr = d.get("modellingApproach")
    print(f"approach   : {appr.get('name') if isinstance(appr, dict) else appr}")
    print(f"publication: {pub.get('title')}")
    print(f"             {pub.get('journal')} {pub.get('year')}  doi:{pub.get('link','')}")
    for f in (d.get("files", {}) or {}).get("main", []):
        print(f"main file  : {f.get('name')} ({f.get('fileSize')} bytes)")


def cmd_files(a):
    d = _get(f"/model/files/{a.id}", {"format": "json"})
    for group in ("main", "additional"):
        for f in d.get(group, []):
            print(f"{group:10s} {f.get('name'):55s} {f.get('fileSize'):>9s}  {f.get('sha256sum','')}")


def _main_filename(mid):
    d = _get(f"/model/files/{mid}", {"format": "json"})
    main = d.get("main", [])
    if not main:
        raise SystemExit(f"bmfetch: {mid} has no main file")
    return main[0]["name"], main[0].get("sha256sum", "")


def cmd_get(a):
    name, want_sha = (a.name, "") if a.name else _main_filename(a.id)
    blob = _get(f"/model/download/{a.id}", {"filename": name}, raw=True)
    os.makedirs(a.o, exist_ok=True)
    dest = os.path.join(a.o, name)
    with open(dest, "wb") as fh:
        fh.write(blob)
    got = hashlib.sha256(blob).hexdigest()
    ok = "" if not want_sha else ("  sha256 OK" if got == want_sha else "  !! SHA256 MISMATCH")
    print(f"{dest}  {len(blob)} bytes  sha256={got}{ok}")
    return dest


def cmd_probe(a):
    for acc in a.accessions.split(","):
        acc = acc.strip()
        if not acc:
            continue
        d = _get("/search", {"query": f"UNIPROT:{acc}", "numResults": a.n, "format": "json"})
        hits = d.get("models", [])
        print(f"\n## {acc}  -> {d.get('matches')} model(s)")
        for m in hits:
            print(f"   {m['id']}\t{m.get('name','')[:70]}")


def main():
    p = argparse.ArgumentParser(prog="bmfetch", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("search"); s.add_argument("query"); s.add_argument("-n", type=int, default=20)
    s.add_argument("--facets", action="store_true"); s.set_defaults(fn=cmd_search)

    s = sub.add_parser("info"); s.add_argument("id"); s.set_defaults(fn=cmd_info)
    s = sub.add_parser("files"); s.add_argument("id"); s.set_defaults(fn=cmd_files)

    s = sub.add_parser("get"); s.add_argument("id"); s.add_argument("-o", default=".")
    s.add_argument("--name", default=None); s.set_defaults(fn=cmd_get)

    s = sub.add_parser("probe"); s.add_argument("accessions"); s.add_argument("-n", type=int, default=10)
    s.set_defaults(fn=cmd_probe)

    a = p.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
