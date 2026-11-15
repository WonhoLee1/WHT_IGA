# 06 — 데이터 버스 mesh_data.sqlite + bezier_grid.mdpa (LIT_NOTE)

## 문헌/계획 상관
- master-plan §2.1: elements(상태/적분점수/C̃^e BLOB),
  quadrature_points(좌표/가중치/법선), model_metadata 스키마.
- bezier_grid.mdpa: Kratos IgaApplication 표준 포맷, C̃^e 희소 주입.
- master-plan §5.4: 텍스트 변환 없는 Rust 고속 직렬화가 본개발 목표 —
  본 PoC는 파이썬 sqlite3로 스키마·라운드트립 계약을 먼저 고정.

## 오픈소스 상관
- **SQLite** (Public Domain, https://sqlite.org) + **rusqlite**(본개발):
  가변 C̃^e·적분점 BLOB 버스.
- **Kratos** ModelPart/mdpa 파서: 노드·요소·NodalData 블록 규격 준용.

## 핵심 알고리즘 (본 PoC)
1. 스키마 생성 → 요소별 C̃^e float64 BLOB 패킹/언패킹.
2. 절단셀 적분점 테이블 기록 → 요소별 조회.
3. 최소 mdpa 작성 → 파서로 노드/요소 수 역확인.

## 본 PoC 범위/한계
- 합성 3요소. 대용량 스트리밍·MPI 파서는 본개발.
