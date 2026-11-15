"""Download images/assets referenced by saved html pages."""
import os, re, time, urllib.parse, pathlib
import requests
from bs4 import BeautifulSoup

BASE = "https://docs.coreform.com/cifa/latest/"
PREFIX = "/cifa/latest/"
OUT_ROOT = pathlib.Path(r"D:\PythonCodeStudy\CAD-CAE-COPILOT\research-coreform-IGA\user-guide")

session = requests.Session()
session.headers.update({"User-Agent": "Mozilla/5.0 research-crawler-assets"})

html_files = sorted(OUT_ROOT.rglob("*.html"))
print(f"html files: {len(html_files)}")

def page_url_for(html_path):
    rel = html_path.relative_to(OUT_ROOT).as_posix()  # e.g. example-problems/foo/foo.html
    if rel.endswith(".html"):
        rel = rel[:-5]
    if rel.endswith("/index") or rel == "index":
        rel = rel[:-5] if rel != "index" else ""
    # url path portion after latest/
    return urllib.parse.urljoin(BASE, rel)

ASSET_EXT = (".png", ".jpg", ".jpeg", ".svg", ".gif", ".webp", ".mp4", ".webm",
             ".pdf", ".stp", ".step", ".zip", ".csv", ".tsv", ".jnl", ".inp", ".cae", ".ico", ".css", ".js")

targets = {}  # asset_url -> local_path
for hp in html_files:
    page_url = page_url_for(hp)
    try:
        soup = BeautifulSoup(hp.read_text(encoding="utf-8"), "html.parser")
    except Exception as e:
        print(f"parse fail {hp}: {e}")
        continue
    refs = []
    for t in soup.find_all(["img", "source", "video", "script"]):
        for attr in ("src", "srcset", "poster"):
            v = t.get(attr)
            if v:
                # srcset may have multiple
                if attr == "srcset":
                    for part in v.split(","):
                        refs.append(part.strip().split(" ")[0])
                else:
                    refs.append(v)
    for t in soup.find_all("link"):
        h = t.get("href")
        if h:
            refs.append(h)
    for a in soup.find_all("a", href=True):
        h = a["href"]
        low = urllib.parse.urldefrag(h)[0].lower()
        if low.endswith(ASSET_EXT):
            refs.append(h)
    for ref in refs:
        if ref.startswith(("data:", "mailto:", "tel:", "javascript:", "#")):
            continue
        abs_url = urllib.parse.urljoin(page_url, ref)
        nofrag = urllib.parse.urldefrag(abs_url)[0]
        if not nofrag.startswith(BASE):
            continue
        # images + docs only: skip jnl/zip/site code, skip generated csv (already have md)
        low = urllib.parse.urlparse(nofrag).path.lower()
        if "/downloads/" in low or "/generated-downloads/" in low:
            continue
        if low.endswith((".jnl", ".zip", ".js", ".css")):
            continue
        if not low.endswith((".png", ".jpg", ".jpeg", ".svg", ".gif", ".webp", ".mp4", ".webm", ".pdf", ".stp", ".step", ".ico")):
            continue
        rel = urllib.parse.urlparse(nofrag).path.split(PREFIX, 1)[1].strip("/")
        local = OUT_ROOT / urllib.parse.unquote(rel)
        targets[nofrag] = local

print(f"unique assets: {len(targets)}")

ok = fail = skip = 0
for i, (url, local) in enumerate(sorted(targets.items()), 1):
    if local.exists() and local.stat().st_size > 0:
        skip += 1
        continue
    local.parent.mkdir(parents=True, exist_ok=True)
    try:
        r = session.get(url, timeout=30)
        if r.status_code == 404:
            print(f"[{i}] 404 {url}")
            fail += 1
            continue
        r.raise_for_status()
        local.write_bytes(r.content)
        ok += 1
        if i % 20 == 0:
            print(f"[{i}/{len(targets)}] ok={ok} fail={fail} skip={skip} last={local.name}")
        time.sleep(0.15)
    except Exception as e:
        print(f"[{i}] FAIL {url}: {e}")
        fail += 1

print(f"DONE ok={ok} fail={fail} skip={skip} total={len(targets)}")
