"""Build combined Coreform IGA user-guide: one .md + one self-contained .html (stdlib only)."""
import re, html, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
UG = ROOT / "user-guide"
EXCLUDE = ("generated-downloads", "site_libs", ".omc", "__pycache__")

GROUPS = [
    "index.md",
    "abaqus_keyword_support_guide.md",
    "install-config/",
    "theory/",
    "user-guide/",
    "cifa-capabilities/",
    "example-problems/",
    "verification-manual/",
]

def collect():
    files = [p for p in UG.rglob("*.md")
             if not any(e in p.parts for e in EXCLUDE)]
    def key(p):
        rel = p.relative_to(UG).as_posix()
        for i, g in enumerate(GROUPS):
            if rel == g or rel.startswith(g):
                return (i, rel)
        return (len(GROUPS), rel)
    return sorted(files, key=key)

MD_OUT = ROOT / "user-guide_combined.md"
HTML_OUT = ROOT / "user-guide_combined.html"

def build_md(files):
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    out = [f"# Coreform IGA User Guide — Combined (md)\n",
           f"- Built: {now} | Sources: {len(files)} md files under `user-guide/` "
           "(excludes `generated-downloads/`, `site_libs/`, `.omc/`)\n",
           "## Contents\n"]
    for i, p in enumerate(files, 1):
        rel = p.relative_to(UG).as_posix()
        anchor = re.sub(r"[^a-z0-9]+", "-", rel.lower()).strip("-")
        out.append(f"{i}. [{rel}](#src-{anchor})")
    out.append("\n---\n")
    for p in files:
        rel = p.relative_to(UG).as_posix()
        anchor = re.sub(r"[^a-z0-9]+", "-", rel.lower()).strip("-")
        txt = p.read_text(encoding="utf-8", errors="replace")
        out.append(f'\n<a id="src-{anchor}"></a>\n## SOURCE: `{rel}`\n\n' + txt.strip() + "\n")
    MD_OUT.write_text("\n".join(out), encoding="utf-8")
    return MD_OUT

# --- minimal md -> html (headings, tables, lists, code fences, links, images, inline) ---
def inline(s):
    s = html.escape(s)
    s = re.sub(r"!\[([^\]]*)\]\(([^)]+)\)", r'<img alt="\1" src="\2" loading="lazy">', s)
    s = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2">\1</a>', s)
    s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
    s = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", s)
    return s

def md_to_html(md):
    lines = md.split("\n")
    h, i, n = [], 0, len(lines)
    in_code, code = False, []
    while i < n:
        ln = lines[i]
        if ln.strip().startswith("```"):
            if in_code:
                h.append("<pre><code>" + html.escape("\n".join(code)) + "</code></pre>")
                code, in_code = [], False
            else:
                in_code = True
            i += 1; continue
        if in_code:
            code.append(ln); i += 1; continue
        # table
        if "|" in ln and i + 1 < n and re.match(r"^\s*\|?[\s:\-|]+\|?\s*$", lines[i + 1]):
            hdr = [c.strip() for c in ln.strip().strip("|").split("|")]
            i += 2; rows = []
            while i < n and "|" in lines[i] and lines[i].strip():
                rows.append([c.strip() for c in lines[i].strip().strip("|").split("|")])
                i += 1
            t = ["<table><thead><tr>" + "".join(f"<th>{inline(c)}</th>" for c in hdr) + "</tr></thead><tbody>"]
            t += ["<tr>" + "".join(f"<td>{inline(c)}</td>" for c in r) + "</tr>" for r in rows]
            h.append("\n".join(t) + "</tbody></table>"); continue
        m = re.match(r"^(#{1,6})\s+(.*)", ln)
        if m:
            lv = len(m.group(1)); t = m.group(2).strip()
            slug = re.sub(r"[^a-z0-9]+", "-", re.sub(r"<[^>]+>", "", inline(t)).lower()).strip("-")
            h.append(f'<h{lv} id="{slug}">{inline(t)}</h{lv}>'); i += 1; continue
        if re.match(r"^\s*([-*+]|\d+\.)\s+", ln):
            items = []
            while i < n and re.match(r"^\s*([-*+]|\d+\.)\s+", lines[i]):
                items.append("<li>" + inline(re.sub(r"^\s*([-*+]|\d+\.)\s+", "", lines[i])) + "</li>")
                i += 1
            h.append("<ul>" + "".join(items) + "</ul>"); continue
        if not ln.strip():
            i += 1; continue
        if ln.strip().startswith("<a id=") or ln.strip().startswith("---"):
            h.append(ln); i += 1; continue
        h.append(f"<p>{inline(ln)}</p>"); i += 1
    if code:
        h.append("<pre><code>" + html.escape("\n".join(code)) + "</code></pre>")
    return "\n".join(h)

def build_html(files):
    toc = []
    sections = []
    for p in files:
        rel = p.relative_to(UG).as_posix()
        anchor = "sec-" + re.sub(r"[^a-z0-9]+", "-", rel.lower()).strip("-")
        toc.append(f'<li><a href="#{anchor}">{html.escape(rel)}</a></li>')
        txt = p.read_text(encoding="utf-8", errors="replace")
        sections.append(f'<section id="{anchor}"><h2 class="src">SOURCE: <code>{html.escape(rel)}</code></h2>\n'
                        + md_to_html(txt) + "\n<p><a href=\"#top\">↑ top</a></p></section>")
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    css = ("body{font-family:'Segoe UI',Arial,sans-serif;max-width:1100px;margin:0 auto;padding:24px;line-height:1.6}"
           "nav.toc{background:#f4f6f8;padding:16px;border-radius:8px}table{border-collapse:collapse;margin:12px 0}"
           "th,td{border:1px solid #bbb;padding:4px 8px}pre{background:#111;color:#eee;padding:12px;overflow:auto;border-radius:6px}"
           "code{background:#eef;padding:0 4px;border-radius:4px}pre code{background:none;padding:0}"
           "section{border-top:2px solid #888;margin-top:32px;padding-top:8px}img{max-width:100%}")
    doc = f"""<!DOCTYPE html><html lang="ko"><head><meta charset="utf-8">
<title>Coreform IGA User Guide — Combined</title><style>{css}</style></head>
<body id="top"><h1>Coreform IGA User Guide — Combined (html)</h1>
<p>Built: {now} | Sources: {len(files)} md files under <code>user-guide/</code>
(excludes <code>generated-downloads/</code>, <code>site_libs/</code>, <code>.omc/</code>).
Images/scripts referenced by original relative paths.</p>
<nav class="toc"><h2>Contents</h2><ol>{''.join(toc)}</ol></nav>
{''.join(sections)}
</body></html>"""
    HTML_OUT.write_text(doc, encoding="utf-8")
    return HTML_OUT

if __name__ == "__main__":
    files = collect()
    print("sources:", len(files))
    md = build_md(files); print("wrote", md, md.stat().st_size)
    ht = build_html(files); print("wrote", ht, ht.stat().st_size)
