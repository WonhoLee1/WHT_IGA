"""Coreform IGA docs crawler - saves under research-coreform-IGA/user-guide"""
import os, re, time, urllib.parse
import requests
from bs4 import BeautifulSoup

BASE = "https://docs.coreform.com/cifa/latest/"
PREFIX = "/cifa/latest/"
OUT_ROOT = r"D:\PythonCodeStudy\CAD-CAE-COPILOT\research-coreform-IGA\user-guide"

session = requests.Session()
session.headers.update({"User-Agent": "Mozilla/5.0 research-crawler"})

visited = set()
to_visit = [BASE]
saved = []

def url_to_path(url):
    parsed = urllib.parse.urlparse(url)
    rel = parsed.path
    # strip prefix
    if PREFIX in rel:
        rel = rel.split(PREFIX, 1)[1]
    rel = rel.strip("/")
    if not rel or rel == "":
        rel = "index"
    if rel.endswith(".html"):
        rel = rel[:-5]
    # safe path
    parts = [re.sub(r'[<>:"|?*]', '_', p) for p in rel.split("/") if p not in ("", ".", "..")]
    if not parts:
        parts = ["index"]
    return os.path.join(OUT_ROOT, *parts)

def fetch(url):
    r = session.get(url, timeout=30)
    r.raise_for_status()
    return r.text

os.makedirs(OUT_ROOT, exist_ok=True)
max_pages = 200

while to_visit and len(visited) < max_pages:
    url = to_visit.pop(0)
    # normalize: drop fragment, ensure trailing handling
    url_nofrag = urllib.parse.urldefrag(url)[0]
    if url_nofrag in visited:
        continue
    if not url_nofrag.startswith(BASE):
        continue
    print(f"[{len(visited)+1}] GET {url_nofrag}")
    try:
        html = fetch(url_nofrag)
    except Exception as e:
        print(f"  FAIL {e}")
        visited.add(url_nofrag)
        continue
    visited.add(url_nofrag)

    soup = BeautifulSoup(html, "html.parser")

    # save files
    base_path = url_to_path(url_nofrag)
    os.makedirs(os.path.dirname(base_path) if "." in os.path.basename(base_path) or "/" in base_path else base_path if False else os.path.dirname(base_path) or OUT_ROOT, exist_ok=True)
    # handle directory case
    if os.path.basename(base_path) == "":
        base_path = os.path.join(base_path, "index")
    html_path = base_path + ".html"
    md_path = base_path + ".md"
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html)
    # markdown-ish: title + text
    title = soup.title.get_text(strip=True) if soup.title else url_nofrag
    # main content
    main = soup.find("main") or soup.find("article") or soup.body
    text = main.get_text("\n", strip=True) if main else soup.get_text("\n", strip=True)
    # links list for provenance
    with open(md_path, "w", encoding="utf-8") as f:
        f.write(f"# {title}\n\nSource: {url_nofrag}\n\n---\n\n{text}\n")
    saved.append((url_nofrag, md_path))

    # discover links
    for a in soup.find_all("a", href=True):
        href = a["href"]
        if href.startswith(("mailto:", "tel:", "javascript:")):
            continue
        abs_url = urllib.parse.urljoin(url_nofrag, href)
        abs_nofrag = urllib.parse.urldefrag(abs_url)[0]
        # only same docs tree, html pages
        if not abs_nofrag.startswith(BASE):
            continue
        # skip assets (pdf, images, css, js)
        if re.search(r'\.(pdf|png|jpg|jpeg|svg|css|js|zip|inp|cae|mp4)(\?.*)?$', abs_nofrag, re.I):
            continue
        if abs_nofrag not in visited and abs_nofrag not in to_visit:
            to_visit.append(abs_nofrag)
    time.sleep(0.3)

print(f"\nDONE visited={len(visited)} saved={len(saved)}")
for u, p in saved:
    print(f"  {u} -> {p}")
