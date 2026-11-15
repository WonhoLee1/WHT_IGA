# -*- coding: utf-8 -*-
"""Generate papers/ MD files from 문헌종합 catalog + verified arXiv metadata."""
import pathlib
OUT = pathlib.Path("research-coreform-IGA/papers")
OUT.mkdir(parents=True, exist_ok=True)

FULL = {  # arXiv id -> (pdf file, pages, size)
 "2602.18130": ("pdf/2602_18130.pdf", 52, 9292030),
 "2509.09236": ("pdf/2509_09236.pdf", 19, 1451921),
 "2106.13728": ("pdf/2106_13728.pdf", 23, 9471140),
 "2107.09024": ("pdf/2107_09024.pdf", 24, 12237904),
 "2105.01340_fix": None,
 "2205.01340": ("pdf/2205_01340.pdf", 31, 4560903),
 "2206.03626": ("pdf/2206_03626.pdf", 23, 1202245),
 "2401.12649": ("pdf/2401_12649.pdf", 19, 7125578),
 "2506.18811": ("pdf/2506_18811.pdf", 33, 13095682),
 "2603.25279": ("pdf/2603_25279.pdf", 31, 867084),
}

P = [
dict(n=1, slug="immersed-iga-deep-report", title="Immersed IGA 심층 연구 보고서 (Bézier Extraction, FCM/CutFEM/IBM 비교, Algoim/무적분/Ghost/AgFEM, 오픈소스 생태계)",
 auth="Gemini Notebook Research Synthesis", year="2026", venue="Technical Report (프로젝트 지식베이스)",
 src="프로젝트 지식베이스", status="메타데이터만 — 원문은 프로젝트 내부 합성 보고서, 별도 PDF 없음",
 core="기하-해석 분리, Bézier Extraction으로 상용 솔버 재활용, FCM/CutFEM/IBM 비교, Algoim·무적분·Ghost·AgFEM·오픈소스(ngsxfem, GridapEmbedded, deal.II) 집대성",
 impl=["Bézier Extraction: NURBS 패치 = Bernstein 기저 × 추출행렬 C, 요소 강성은 표준 Bernstein 루틴으로 조립 후 C로 변환 (솔버 수정 불필요).","워크플로우: B-Rep/복셀 입력 → 배경 직교격자 침강 → 절단셀 분류 → 적분(Algoim/옥트리/무적분) → 안정화(Ghost/AgFEM) → 표준 조립/해석.","적응 세분화는 경계 교차셀 중심 트리 구조로 국소 해상도 확보."],
 link="Coreform IGA 전체 아키텍처의 개념 지도 — 아래 P02~P30 각 모듈의 상위 설계도.",
 repro=["추출행렬 C가 단위분할을 유지하는지 단위 테스트","동일 모델을 경계적합 격자 결과와 변형에너지 비교"]),
dict(n=2, slug="schillinger2012-design-through-analysis", title="An Isogeometric design-through-analysis methodology based on adaptive hierarchical refinement of NURBS, immersed boundary methods, and T-spline CAD surfaces",
 auth="Schillinger, Dedè, Scott, Evans, Borden, Rank, Hughes", year="2012", venue="CMAME 249–252, pp.116–150",
 src="ResearchGate Pub.256918134", doi="10.1016/j.cma.2012.03.017", status="메타데이터만 — 유료 저널, 기관 접근 필요",
 core="계층적 B-스플라인 + FCM 결합, T-스플라인 산업 CAD(프로펠러·휠) 직접 해석, Design-through-Analysis 최초 정립",
 impl=["계층적 세분화: 트리 구조로 경계 피처 셀만 국소 세분, 선형독립성 유지.","FCM 지시함수: 물리영역 내 α=1, 허구영역 α=ε(≈1e-8~1e-10)로 약재료 처리 후 표준 적분.","T-스플라인 표면을 배경셀과 교차 판정하여 내/외/절단 분류."],
 link="Coreform 계열(Scott→Coreform)의 원류 — 트리밍된 CAD 직접 해석 철학의 출발점.",
 repro=["허구영역 α 감도 스터디","세분화 레벨별 변형에너지 수렴 곡선"]),
dict(n=3, slug="messmer2022-trimmed-solids", title="Efficient CAD-integrated isogeometric analysis of trimmed solids",
 auth="Meßmer, Teschemacher, Leidinger, Wüchner, Bletzinger", year="2022", venue="CMAME 400, 115584",
 src="TUM mediaTUM 1767547", doi="10.1016/j.cma.2022.115584", status="메타데이터만 — 유료 저널",
 core="B-Rep 솔리드 → 트리밍 3변량 B-스플라인 배경격자 매핑 파이프라인, QuESo 전처리, 선형 FEM 대비 12배 적은 자유도로 동등 정확도",
 impl=["QuESo(Kratos) 전처리: B-Rep 면과 배경셀 교차 → 절단셀 적분점/가중치 생성.","트리밍된 체적은 경계 B-Rep 서피스 적분으로 체적항을 변환하는 발산정리 경로와 병행 가능.","자유도당 효율 지표(변형에너지 오차/자유도)로 검증."],
 link="B-Rep 직접 입력 파이프라인의 실용 구현 레퍼런스.",
 repro=["동일 B-Rep에서 자유도-오차 곡선 재현","QuESo 버전·파라미터 기록"]),
dict(n=4, slug="divi2021-scan-based", title="Scan-based immersed isogeometric analysis (PhD thesis)",
 auth="Sai Chandana Divi (Verhoosel/Auricchio/Reali/van Brummelen)", year="2021/2022", venue="TU/e & Univ. Pavia PhD thesis, ISBN 978-90-386-5469-0",
 src="Univ. Pavia IRIS", status="메타데이터만 — 학위논문 공개본 기관 저장소 확인 필요",
 core="복셀 스캔→절단셀 Immersed IGA, Strang 제1보조정리 기반 오차추정 적응 옥트리, 오일러지표 위상오류 방지 THB 세분화",
 impl=["스캔 복셀 레벨셋을 배경 THB-스플라인 격자에 침강.","셀별 오차지표 기반 옥트리 적분 깊이 적응 조절.","필터링으로 생기는 고립 위상 성분은 오일러지표 검사로 제거 후 해석."],
 link="이미지 기반 해석(의료·적층제조 CT)의 이론 근거.",
 repro=["복셀 해상도별 수렴 스터디","위상검사 on/off 비교"]),
dict(n=5, slug="loibl2026-quadrature-comparison", title="Comparative study of different quadrature methods for cut elements",
 auth="Loibl, Teixeira, Toprak, Shishkina, Miao, Kiendl, Marussig, Kummer", year="2026", venue="Arch. Comput. Methods Eng. / arXiv:2602.18130",
 src="arXiv math.NA", arxiv="2602.18130", status="FULL-TEXT 확보 — pdf/2602_18130.pdf (52pp, ~9.3MB, arXiv 검증: 2026-02-20)",
 core="Algoim·BoSSS·FCMLab·QUGaR·Gridap/STLCutters·mlhp·ngsxfem·Nutils 망라 2D/3D 벤치마크: 정확도·적분점수·효율 비교, 모멘트피팅/고차매핑/높이함수 축소 정량 비교",
 impl=["벤치마크 축: 기하 난이도(곡률·박판·특이점) × 요구 정확도 × 적분점 수.","모멘트피팅: 단항식 모멘트 일치 가중치 해법, 점 수는 적지만 조건수 관리 필요.","높이함수 축소(Saye): 재귀 차원축소 + 1D Gauss-Legendre/Tanh-Sinh 텐서곱, 기계정밀도 경로.","실무 선택표: 곡면 B-Rep이면 축소법, 다면체 근사면이면 옥트리/피팅."],
 link="적분 백엔드 선택의 근거표 — Coreform 메시 생성 옵션 해석에 직접 활용.",
 repro=["2D 원/3D 구 벤치마크 재현","적분점수-오차 파레토 곡선"]),
dict(n=6, slug="antolin2022-quadrature-free", title="Quadrature-free immersed isogeometric analysis",
 auth="Antolin, Hirschler", year="2022", venue="Eng. with Computers 38(5):4475–4499 / arXiv:2107.09024",
 src="arXiv + ResearchGate 360184654 + EPFL Infoscience", arxiv="2107.09024", doi="10.1007/s00366-022-01644-3",
 status="FULL-TEXT 확보 — pdf/2107_09024.pdf (24pp, ~12.2MB, arXiv 검증: 2021-07-19)",
 core="절단셀 적분점 제로: Galerkin 텐서를 일관 다항식 공간으로 변환 → 발산정리 2회 적용 → 3D 체적=2D B-Rep 면적분=1D 경계 선적분, 닫힌해 기계정밀도",
 impl=["피적분 다항식 차수 p에서 체적항을 경계 선적분으로 축소하는 변환 테이블 사전 구축.","B-Rep 트림 곡면의 파라미터 경계 루프를 따라 1D Gauss 적분만 수행.","구현 조건: 커널이 다항식이어야 하며, 비다항식 재료·하중은 투영 후 적용."],
 link="Coreform 무적분 경로의 이론 원형 — 적분점 폭발이 문제될 때 대안.",
 repro=["구/원기둥 체적·질량행렬 기계정밀도 테스트","다항식 차수별 변환 정확도"]),
dict(n=7, slug="saye-algoim", title="Algoim: Algorithms for implicitly defined geometry, level set methods, and Voronoi implicit interface methods",
 auth="Robert I. Saye (LBNL/UC Berkeley)", year="2015/2022", venue="SISC 2015; JCP 2022 + C++ 오픈소스 (algoim.github.io)",
 src="Algoim 공식 포털", doi="10.1137/140966290; 10.1016/j.jcp.2021.110720", status="메타데이터만 + 코드 공개 — 깃허브 클론으로 확보 가능",
 core="초직사각형 내 레벨셋 절단영역을 다가지 높이함수 그래프로 재구성 → 재귀 차원축소 → 1D 적분 텐서곱, 고차 가중치 기계정밀도; deal.II·PETSc 백엔드 채택",
 impl=["셀 내 레벨셋 부호로 위상 분류 → 높이함수 방향 자동 선택.","특이 방향은 재귀 세분 후 축소 반복.","출력: 점 좌표+가중치 테이블 → 기존 FEM 조립 루프에 그대로 주입."],
 link="절단셀 적분의 사실상 표준 라이브러리 — Coreform 메시 적분과 상호 검증용.",
 repro=["레벨셋 구/토러스 적분 오차","deal.II 튜토리얼 스텝 재현"]),
dict(n=8, slug="quadrature-on-implicit-regions-julia", title="QuadratureOnImplicitRegions (Julia package)",
 auth="Julia 전산역학 커뮤니티", year="2023", venue="Julia Packages Registry",
 src="Julia Packages", status="메타데이터만 + 코드 공개 — Pkg.add로 확보 가능",
 core="Julia JIT 환경에서 레벨셋/암시적 영역 고차 적분룰 자동 생성, 비적합 해석 라이브러리 연동 인터페이스",
 impl=["Saye류 축소법을 Julia로 래핑, 다중 디스패치로 기하 유형별 분기.","GridapEmbedded 등과의 파이프라인 예제 포함."],
 link="빠른 프로토타이핑용 — 아이디어 검증 후 C++/Coreform 경로로 이식.",
 repro=["Julia 재현 스크립트 버전 고정","적분룰 차수 감도"]),
dict(n=9, slug="thesis-curved-geometry", title="High-Order Spatial Discretization and Numerical Integration Schemes for Curved Geometries (PhD)",
 auth="ProQuest 학위논문", year="미상", venue="ProQuest Dissertations",
 src="ProQuest", status="미확보 — 서지정보 불완전(저자·연도 미상), ProQuest 접근 필요",
 core="곡면 경계 비적합 FEM의 기하 포화오차 차단용 고차 이산화·서브셀 적분",
 impl=["곡면 근사 차수 < 해 근사 차수이면 수렴 정체 → 경계셀 기하 차수 상향이 처방.","서브셀 분할 깊이와 기하 차수를 함께 보고."],
 link="곡면 정밀도 문제 발생 시 점검표.",
 repro=["기하 차수 고정-해 차수 상향 실험"]),

dict(n=10, slug="boundary-conformal-finite-elasticity", title="Immersed isogeometric analysis with boundary conformal quadrature for finite deformation elasticity",
 auth="PMC/Comput. Mech. 연구진", year="2023", venue="Comput. Mech. / PMC12399451",
 src="PMC OA", status="메타데이터만 — PMC ID로 전문 확보 가능 (이번 배치 미수집)",
 core="유한변형 초탄성에서 경계정합 적분으로 적분오차·체적잠김 방지",
 impl=["절단셀 기하를 경계 추종 서브분할로 정합 후 적분.","비압축성 근방은 혼합정식화와 병행 검토."],
 link="대변형·고무류 해석 시 적분 처방.",
 repro=["Neo-Hookean 패치테스트","압축 잠김 여부 확인"]),
dict(n=11, slug="burman2010-ghost-penalty", title="Ghost penalty",
 auth="Erik Burman", year="2010", venue="C. R. Math. 348(21–22):1217–1220",
 src="Numdam OA", doi="10.1016/j.crma.2010.10.006", status="메타데이터만 — Numdam OA로 확보 가능 (이번 배치 미수집)",
 core="Ghost Penalty 원조: 절단셀-내부셀 공유면의 법선 고차도함수 도약 페널티 → 조건수 상한 κ≤Ch^{-2p}, 강압성 증명",
 impl=["페널티항: Σ_F γ h^{2k-1} ‖⟦∂_n^k u⟧‖²_F, k=1..p.","γ는 0.01~0.1 범주에서 시작, 조건수-정확도 트레이드.","구현 위치: 면 루프에서 도약행렬만 추가, 나머지는 표준 조립."],
 link="소형 절단셀 안정화의 출발점 — 모든 후속 기법의 베이스라인.",
 repro=["극소 체적분율(1e-6) 조건수 테스트","γ 스윕"]),
dict(n=12, slug="badia2022-linking-ghost-agfem", title="Linking ghost penalty and aggregated unfitted methods",
 auth="Badia, Neiva, Verdugo", year="2022", venue="CMAME 388, 114232 / arXiv:2106.13728",
 src="arXiv math.NA", arxiv="2106.13728", doi="10.1016/j.cma.2021.114232",
 status="FULL-TEXT 확보 — pdf/2106_13728.pdf (23pp, ~9.5MB, arXiv 검증: 2021-06-19)",
 core="Ghost↔AgFEM 수학적 연결 최초 규명: 페이스기반 페널티의 무한대 극한 잠김 규명 → Locking-Free Ghost Penalty 정립",
 impl=["집합(aggregate) 구속: 불량 절단 DOF를 내부셀 DOF의 외삽으로 제거(강한 안정화).","새 페널티는 해-집합확장 거리 제어 → γ→∞에서도 잠김 없음.","구현: 셀 집합 그래프 구축 → 구속행렬 적용 → 축소 시스템 해석."],
 link="γ 튜닝 부담을 없애는 실무 해법.",
 repro=["γ→∞ 극한 수렴 테스트","집합 패턴 시각화"]),
dict(n=13, slug="burman2025-locking-free-design", title="On the design of locking free ghost penalty stabilization and the relation to CutFEM with discrete extension",
 auth="Burman, Hansbo, Larson", year="2025/2026", venue="Numer. Math. 158:249–275 / arXiv:2205.01340",
 src="Springer OA + arXiv", arxiv="2205.01340", doi="10.1007/s00211-025-01502-6",
 status="FULL-TEXT 확보 — pdf/2205_01340.pdf (31pp, ~4.6MB, arXiv 검증: 2022-05-03)",
 core="이산확장 CutFEM 공간–페널티 커널 관계 일반화, 임의 고차에서도 잠김 없는 함수형(DOF기반) Ghost Penalty 설계 완성",
 impl=["셀 내부 DOF 기반 벌점(함수형)으로 설계 자유도 확보.","차수 p별 스케일링 h^{2k-1} 유지, 상수 γ의 p 의존성 표기.","이산확장 연산자와의 동치성으로 안정성 증명 재사용."],
 link="고차(p≥3) 안정화 설계 지침.",
 repro=["p=1..4 조건수-오차 매트릭스"]),
dict(n=14, slug="gross2020-cutfem-preconditioner", title="Analysis of optimal preconditioners for CutFEM",
 auth="Groß, Reusken", year="2020", venue="RWTH IGPM Report 519",
 src="RWTH 저장소", status="메타데이터만 — 기관 리포트 공개본 확인 필요",
 core="Ghost 적용 CutFEM의 최적 전처리기 스펙트럼 수렴 증명: 배경기저+절단인접 국소공간 부분공간분해 → 격자·절단위치 독립 다중격자 최적수렴",
 impl=["두 부분공간 가산 슈바르츠/다중격자 스무더 구성.","국소공간은 절단면 인접 셀 집합으로 한정 → 비용 절감."],
 link="대규모 문제 반복해법 수렴 확보용.",
 repro=["절단위치 스위프 반복횟수","격자 독립성 확인"]),
dict(n=15, slug="ghost-penalty-survey", title="Ghost Penalty Stabilization in Unfitted FEM (survey)",
 auth="Emergent Mind synthesis", year="2026", venue="Topic synthesis",
 src="Emergent Mind KB", status="2차 종합자료 — 원문 대체 불가, 입문용으로만 사용",
 core="페이스도약·구배투영·질량스케일링·AgFEM 결합 등 이론·FSI 응용 종합",
 impl=["명시동역학 질량스케일링: 소형셀 질량행렬 하한 보정.","애플리케이션별 페널티 선택 가이드."],
 link="입문·교육용.",
 repro=["서베이 인용 원문 대조 필수"]),
dict(n=16, slug="badia2022-space-time-agfem", title="Space-time AgFEM",
 auth="Badia, Neiva, Verdugo", year="2022", venue="CMAME / arXiv:2206.03626",
 src="arXiv math.NA", arxiv="2206.03626",
 status="FULL-TEXT 확보 — pdf/2206_03626.pdf (23pp, ~1.2MB, arXiv 검증: 2022-06-08, 제목 주의: 실제 Space-time unfitted FEM)",
 core="시공간 도메인 전체에 AgFEM 확장: 시간슬랩별 공간 이산확장 → 이동경계도 무조건 안정·최적수렴, 페널티 튜닝 불필요",
 impl=["시간슬랩마다 그 순간의 절단 위상으로 집합 재구성.","시공간 업윈딩과 결합 시 안정성 유지.","이동경계 문제는 슬랩 경계에서 확장 연산자 전달."],
 link="이동경계·ALE 대안으로서 시공간 비적합법.",
 repro=["이동 원 벤치마크","슬랩 크기 감도"]),
dict(n=17, slug="badia2024-moving-explicit-geometry", title="Space-time unfitted finite elements on moving explicit geometry representations",
 auth="Badia, Martorell, Verdugo", year="2024", venue="JCP 511, 113127 / arXiv:2401.12649",
 src="arXiv comp-phys", arxiv="2401.12649", doi="10.1016/j.jcp.2024.113127",
 status="FULL-TEXT 확보 — pdf/2401_12649.pdf (19pp, ~7.1MB, arXiv 검증: 2024-01-23)",
 core="명시적 다면체 B-Rep 이동계면에 시공간 비적합 FEM + 고차 안정화 파이프라인",
 impl=["B-Rep↔배경격자 교차의 폴리토프 분할 → 2단계 적분메시.","이동 시 위상변화(합병/분리) 이벤트의 슬랩 경계 처리."],
 link="STL/B-Rep 직접 입력 동적 문제의 기하 커널.",
 repro=["이동계면 질량보존 오차","위상변화 케이스"]),
dict(n=18, slug="thomas2022-u-splines", title="U-splines: splines over unstructured meshes",
 auth="Thomas, Engvall, Schmidt, Tew, Scott", year="2022", venue="CMAME 403, 115747 / Coreform·BYU",
 src="Coreform 기술자료", doi="10.1016/j.cma.2022.115747", status="메타데이터만 — 유료 저널, Coreform 라이브러리 요약본 확인 권장",
 core="비구조격자 U-스플라인: h·p·연속성 국소 가변, T-접합·혼합셀 지원, 국소 영공간 해법으로 분할단위·양성·국소선형독립 증명",
 impl=["셀별 요구 연속성(r,θ) → 베지어 계수 구속식 → 국소 영공간 기저 추출.","추출행렬은 희소 저장, 요소 루프는 Bernstein 루틴 재활용(P01 연계).","산업 B-Rep → 해석적합 U-스플라인 변환이 Coreform 전처리 핵심."],
 link="Coreform IGA 기저 기술 그 자체 — 최우선 정독 대상.",
 repro=["단위분할·양성 수치 검증","특이점 패치테스트"]),
dict(n=19, slug="thomas2024-extraordinary-spectrum", title="Spectrum analysis of C0, C1, and G1 constructions for extraordinary points",
 auth="Thomas, Scott et al.", year="2024", venue="CMAME",
 src="ResearchGate 382998545", status="메타데이터만 — 유료 저널",
 core="특이점 근방 연속성 구조의 스펙트럼 비교: Kirchhoff-Love 셸 등 고차 PDE의 진동·오염 방지 기저 선택 원리",
 impl=["고유값 스펙트럼의 이상치(outlier) 모드 탐지 → 연속성 구조별 비교.","G1이 기하 엄밀성, C1이 스펙트럼 청정성에서 강점인 트레이드."],
 link="셸·고차방정식 해석 시 기저 선택 근거.",
 repro=["자유진동 스펙트럼 재현"]),
dict(n=20, slug="martorell2023-intersecting-meshes", title="Intersecting surface meshes and Cartesian background meshes for embedded simulations",
 auth="Martorell, Badia, Verdugo", year="2022/2023", venue="UPC UPCommons Report",
 src="UPCommons", status="메타데이터만 — 기관 리포트 공개본 확인 필요",
 core="STL 표면메시×데카르트 배경격자 교차의 폴리토프 셀 분할·2단계 적분메시 체계화",
 impl=["삼각형-셀 교차 → 폴리토프 분해 → 테트라 분할 후 표준 적분.","2단계: 1차 위상분류 메시 + 2차 고정밀 적분메시."],
 link="STL 직접 입력 경로의 기하 커널.",
 repro=["복잡 STL 교차 견고성 테스트"]),
dict(n=21, slug="hoang2018-immersogeometric-flow", title="Isogeometric and immersogeometric analysis of incompressible flow problems (PhD)",
 auth="Tuong Hoang (Verhoosel/van Brummelen/Reali)", year="2018", venue="TU/e & IUSS PhD, ISBN 978-90-386-4449-3",
 src="TU/e Pure", status="메타데이터만 — 학위논문 공개본 확인 필요",
 core="비압축성 NS/Stokes Immersogeometric: 속도-압력 동일 고차 B-스플라인 + 스켈레톤 안정화 + 속도 Ghost Penalty, 다공성 유동 실증",
 impl=["스켈레톤 안정화: 요소 골격 위 고차도함수 도약 penalize → inf-sup 우회.","등차수(Equal-order) 속도-압력 공간 사용 가능.","다공성: 고체 골격은 FCM 지시함수로 처리."],
 link="유동·다공성 확장 시 안정화 처방.",
 repro=["원주 유동 벤치마크","압력 진동 여부"]),
dict(n=22, slug="rim2025-enriched-ibm-cht", title="Enriched immersed boundary method for convective conjugate heat transfer with level set-based surface reconstruction",
 auth="Rim, Zhao, Yan", year="2025", venue="AIMS ACSE 3:1–23",
 src="AIMS Press", doi="10.3934/acse.2025002", status="메타데이터만 — OA 저널이나 이번 배치 PDF 확보 실패(500), 재시도 필요",
 core="열유속 양립성 문제 해결: 절단셀 온도 자유도를 유체/고체에 복제 부여(EIBM) + 레벨셋 표면재구성",
 impl=["절단셀 Heaviside 풍부화: T = H·T_f + (1-H)·T_s 형태 자유도 복제.","계면에 Nitsche/페널티로 온도·유속 연속성 약부과.","레벨셋 재구성 표면에서 열전달계수 평가."],
 link="공액열전달(CHT) 문제의 절단셀 처방.",
 repro=["평판 대류 벤치마크 Nusselt 비교"]),
dict(n=23, slug="frei2026-eulerian-fsi", title="Numerical Analysis of a Cut Finite Element Approach for Fully Eulerian Fluid-Structure Interaction with Fixed Interface",
 auth="Frei, Knoke, Steinbach, Wenske, Wick", year="2026", venue="arXiv:2603.25279 [math.NA]",
 src="arXiv", arxiv="2603.25279",
 status="FULL-TEXT 확보 — pdf/2603_25279.pdf (31pp, ~0.9MB, arXiv 검증: 2026-03-26)",
 core="유체+구조 단일 오일러 배경격자 CutFEM 완전오일러 FSI: 체적절단비 적응 가중 Ghost Penalty로 계면접촉·연성 강건화",
 impl=["가중 페널티: 절단비 η→0 셀로 갈수록 가중 상향.","고정계면 가정 하 오차수렴 분석 → 이동계면 확장 전제.","유체·고체 서로 다른 구성식을 동일 격자 위 Heaviside 분기로 처리."],
 link="FSI 확장 시 오일러 기술 경로.",
 repro=["고정계면 수렴율 재현"]),
dict(n=24, slug="cahaly2025-pcic-vof", title="PCIC: Cylindrical Volume Moment Calculation and Interface Reconstruction for Sub-Grid Modeling in VOF Methods",
 auth="Cahaly, Wasquel, Zou, Desjardins, Evrard", year="2026", venue="JCP / arXiv:2506.18811",
 src="arXiv comp-phys", arxiv="2506.18811",
 status="FULL-TEXT 확보 — pdf/2506_18811.pdf (33pp, ~13.1MB, arXiv 검증: 2025-06-23)",
 core="VOF 아격자 계면재구성용 원통형 체적모멘트 고정밀 산출 PCIC",
 impl=["원통 좌표 체적모멘트의 해석식+수치적분 하이브리드.","아격자 계면 다항식 재구성 → 곡률·법선 정밀도 향상."],
 link="다상유동 계면모멘트 계산의 기하 커널(연계 참고).",
 repro=["구/원통 모멘트 해석해 대조"]),
dict(n=25, slug="ruizbaier-darcy-divfree", title="Divergence-free unfitted finite element discretisations for the Darcy problem",
 auth="Ruiz-Baier, de Prenter, Larson et al.", year="2024", venue="Monash/Authors Monograph",
 src="저자 저장소", status="메타데이터만 — 모노그래프 공개본 확인 필요",
 core="Darcy 문제 발산보존 비적합 이산화: H(div) 적합공간+AgFEM 확장으로 무진동 체적보존",
 impl=["Raviart-Thomas류 H(div) 공간을 절단셀에 확장.","압력-유속 inf-sup은 AgFEM 구속과 병행."],
 link="다공성·침투류 체적보존 처방.",
 repro=["발산오차 기계정밀도 테스트"]),
dict(n=26, slug="teixeira2025-iga-topopt", title="Isogeometric Topology Optimization Based on Topological Derivatives",
 auth="Teixeira, Krenn, Gangl, Marussig", year="2025", venue="Struct. Multidisc. Optim. / arXiv:2509.09236",
 src="arXiv", arxiv="2509.09236",
 status="FULL-TEXT 확보 — pdf/2509_09236.pdf (19pp, ~1.5MB, arXiv 검증: 2025-09-11)",
 core="위상도함수+IGA 배경격자 재격자 없는 위상최적화: Algoim 체적비 산출 + 절단강성 완충전이",
 impl=["위상도함수 필드로 구멍 생성 위치 결정 → 레벨셋 업데이트.","절단셀 강성은 체적비 η의 완충함수 E(η)로 전이(0/1 불연속 완화).","Algoim으로 η 고정밀 산출이 수렴 핵심."],
 link="Coreform 결과→경량설계 연계 시 방법론.",
 repro=["캔틸레버 컴플라이언스 벤치마크"]),
dict(n=27, slug="bath2025-qmc-helmholtz", title="Quasi-Monte Carlo methods for UQ of Helmholtz wave propagation and scattering",
 auth="Bath PureHost 연구진", year="2025", venue="PureHost Monograph",
 src="Bath PureHost", status="메타데이터만 — 기관 모노그래프 확인 필요",
 core="Helmholtz 고주파 파동·산란의 매질/경계 불확실성을 QMC로 효율 정량화하는 고차 FEM",
 impl=["QMC 격자점(랭크-1 격자)으로 몬테카를로 N^{-1/2}보다 빠른 수렴.","고주파 오염(pollution) 대비 고차 공간이산화 병행."],
 link="진동·음향 UQ 확장 참고.",
 repro=["산란 단면적 분산 수렴"]),
dict(n=28, slug="zonta2022-landau-petsc", title="The Landau Collision Integral in the Particle Basis in the PETSc Library",
 auth="Zonta, Pusztay, Hirvijoki, Adams, Knepley", year="2022/2023", venue="Phys. Plasmas 29(12):123906",
 src="ResearchGate 371786455", doi="10.1063/5.0105182", status="메타데이터만 — 유료 저널",
 core="Landau 충돌적분자의 보존형(FE/입자) 이산화 + Algoim 다변수 다항식 암묵영역 적분을 PETSc 백엔드 통합",
 impl=["에너지·운동량·양성 보존 이산화.","Algoim 적분 커널을 PETSc 분산 조립에 연결."],
 link="Algoim이 플라즈마까지 쓰인다는 성숙도 증거.",
 repro=["보존량 장시간 추적"]),
dict(n=29, slug="antolin2021-qfree-preprint", title="Quadrature-free Immersed Isogeometric Analysis (arXiv preprint)",
 auth="Antolin, Hirschler", year="2021", venue="arXiv:2107.09024",
 src="arXiv", arxiv="2107.09024", status="P06과 동일 원문 — pdf/2107_09024.pdf로 확보済, 별도 파일 미생성",
 core="P06의 원본 프리프린트 (내용 동일)",
 impl=["P06 참조."],
 link="P06 참조.",
 repro=["P06 참조"]),
dict(n=30, slug="antolin2022-qfree-infoscience", title="Quadrature-free immersed isogeometric analysis (EPFL Infoscience)",
 auth="Antolin, Hirschler", year="2022", venue="EPFL Infoscience OA",
 src="EPFL 기관저장소", arxiv="2107.09024", status="P06과 동일 원문 — pdf/2107_09024.pdf로 확보済, 별도 파일 미생성",
 core="P06의 EPFL 정식 등재 OA 버전 (내용 동일, 구현 지침 포함)",
 impl=["P06 참조 + B-Rep 선적분 루틴 구현 지침."],
 link="P06 참조.",
 repro=["P06 참조"]),
]

for p in P:
    n = p["n"]
    lines = []
    A = lines.append
    A(f"# P{n:02d} — {p['title']}")
    A("")
    A(f"- 저자: {p['auth']}")
    A(f"- 연도/게재: {p['year']} / {p['venue']}")
    if "doi" in p: A(f"- DOI: {p['doi']}")
    if "arxiv" in p: A(f"- arXiv: {p['arxiv']}")
    A(f"- 출처: {p['src']}")
    A(f"- 확보상태: {p['status']}")
    A("")
    A("## 핵심 수치 기여")
    A("")
    A(p["core"])
    A("")
    A("## 핵심 구현 내용 추출")
    A("")
    for b in p["impl"]:
        A(f"- {b}")
    A("")
    A("## Coreform IGA 연관성")
    A("")
    A(p["link"])
    A("")
    A("## 재현 체크리스트")
    A("")
    for r in p["repro"]:
        A(f"- [ ] {r}")
    (OUT / f"P{n:02d}_{p['slug']}.md").write_text("\n".join(lines), encoding="utf-8")

# index
idx = ["# Papers — 문헌 확보 현황 및 핵심 구현 인덱스", "",
 "출처: ../문헌종합.md (번호 P01–P30은 해당 문서의 1–30번과 일치)", "",
 "## 확보 요약", "",
 "- FULL-TEXT 9편: P05, P06(=P29, P30), P12, P13, P16, P17, P23, P24, P26 → `pdf/` (arXiv API 검증済)",
 "- 메타데이터+공개코드: P07(Algoim), P08(Julia)",
 "- 메타데이터만(유료·기관저장소): P02, P03, P04, P09, P10, P11, P14, P15(2차), P18, P19, P20, P21, P22(OA 재시도), P25, P27, P28",
 "- P29/P30은 P06과 동일 원문이므로 별도 파일 없이 P06 PDF로 커버", "",
 "## 파일 목록", ""]
for p in P:
    idx.append(f"- P{p['n']:02d} `{p['slug']}` — {p['year']} — {p['status'][:60]}")
idx += ["", "## 구현 스택별 읽기 순서 (상세는 implementation_synthesis.md)",
 "- 기하입력: P02 → P03 → P20 → P04",
 "- 적분: P07 → P05 → P06",
 "- 안정화: P11 → P12 → P13 → P14(+P16 시공간)",
 "- 기저: P18 → P19",
 "- 다중물리: P21 → P22 → P23 → P25 → P24",
 "- 최적화·UQ: P26 → P27 → P28",
 "- 시공간 이동계면: P16 → P17"]
(OUT / "00_index.md").write_text("\n".join(idx), encoding="utf-8")

syn = """# 핵심 구현 내용 종합 (Implementation Synthesis)

## 전체 스택 (기하 → 적분 → 안정화 → 기저 → 해법)
1. **기하 입력**: B-Rep/STL/복셀을 배경 직교격자에 침강 (P02, P03, P20, P04).
   분류(내/외/절단) → 절단셀만 고정밀 처리.
2. **절단셀 적분** (택1/혼합): Algoim 높이함수 축소 (P07, 비교표 P05),
   무적분 발산정리 2회 축소 (P06), 옥트리/모멘트피팅 (P05, P09).
3. **소형셀 안정화**: Ghost 페널티 면도약항 (P11),
   Locking-Free 확장 (P12, P13), AgFEM 집합구속, 전처리기 부분공간분해 (P14),
   시공간 슬랩 확장 (P16, P17).
4. **스플라인 기저**: U-스플라인 국소 영공간 추출행렬 (P18),
   특이점 연속성 스펙트럼 선택 (P19), Bézier Extraction으로 상용솔버 재활용 (P01).
5. **다중물리**: 스켈레톤 안정화 유동 (P21), EIBM 열유속 복제자유도 (P22),
   완전오일러 FSI 가중페널티 (P23), Darcy 발산보존 (P25), VOF PCIC 모멘트 (P24).
6. **설계·불확실성**: 위상도함수+Algoim 체적비 (P26), QMC Helmholtz UQ (P27),
   Landau-PETSc 보존형 적분 (P28).

## 최소 재현 파라미터셋
- FCM 허구강성 α=1e-8~1e-10, Ghost γ=0.01~0.1 시작 후 조건수 스윕 (P02, P11).
- 적분: 곡면 B-Rep=축소법, 다면체=옥트리/피팅 (P05).
- 고차(p≥3)는 함수형 Ghost (P13), 이동계면은 슬랩별 집합 재구성 (P16).
- 체적비 η 전이는 완충함수 E(η)로 불연속 완화 (P26).

## Coreform IGA 매핑
- 전처리(B-Rep→배경격자→절단분류)는 P03+P20 파이프라인.
- 적분 백엔드는 P07 표준 + P06 무적분 대안.
- 안정화는 P12/P13 Locking-Free 계열.
- 기저는 P18 U-스플라인 + P01 Bézier Extraction.
"""
(OUT / "implementation_synthesis.md").write_text(syn, encoding="utf-8")
print(f"wrote {len(P)} papers + index + synthesis")
