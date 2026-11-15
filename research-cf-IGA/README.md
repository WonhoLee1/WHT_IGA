# research-coreform-IGA

Source: https://docs.coreform.com/cifa/latest/
Collected: 2026-11-14
Counts: 111 pages (.html 원본 + .md 텍스트), 그림/첨부 857개 (png 483, svg 304, pdf 70), css 10 + js 19, 수치데이터 csv/tsv 76 + stp 1, user-guide/ 하위 저장

## 구조 (user-guide/ 기준 상대경로)
- `index.md` — Welcome (전체 목차)
- `install-config/` — Installation & Configuration 6편
- `user-guide/` — User Guide 본체 (overview, workflow, capabilities, using-cifa-cae, response-spectrum 3편, roadmap, release-notes 5편, scripting-guide 4편, python-api 9편)
- `theory/` — Theory + consistent-unit-systems
- `example-problems/` — Example 5편 (thin-plate, adaptivity, glass-to-metal-seal, cook-membrane, door-lock)
- `verification-manual/` — Verification 12편 (linear-elastic 7, contact 1, thermal-strain 1, shock-vibration 1, implicit 1, heat-transfer 3)
- `abaqus_keyword_support_guide.md`, `cifa-capabilities/index.md`
- `generated-downloads/`, `example-problems/.../data/` — CSV/TSV 수치 데이터 (검증용)
- `crawl.py` — 수집 스크립트 (상위 폴더)

## 메모
- `.jnl` 다운로드는 404로 실패 (사이트 결함) — md/html 본문은 정상.
- `?card=` 쿼리 페이지는 동일 index로 수렴 — 정상.
- 원문 그대로 보존, 인코딩 UTF-8.
