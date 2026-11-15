# 07 — 프로브 역매핑 + 듀얼 후처리 스케치 (LIT_NOTE)

## 문헌/계획 상관
- master-plan Track A/B/C: VTKHDF(ParaView) / .fem+.hvascii→hvtrans→.h3d /
  probe_results.json.
- Coreform API 패러다임 §6.1: coreform.probe (점/선/면 프로브, XY 플로팅).

## 오픈소스 상관
- **splinepy** 역매핑(Point Inversion) → 본 PoC는 텐서곱 패치 Newton 역매핑 직접 구현.
- **h5py** (BSD): control_points.hdf5 변위벡터 저장 → ParaView VTKHDF의 전단계.
- **ryu** (본개발 Rust): 고속 부동소수점 포맷 .hvascii 작성 → 본 PoC는 텍스트 규격만 고정.

## 핵심 알고리즘 (본 PoC)
1. Newton 역매핑: 야코비 2x2 해석식으로 (u,v) 복원, 기계정밀도 확인.
2. hdf5 변위 라운드트립.
3. 최소 hvascii 형식 블록 작성·파싱.

## 본 PoC 범위/한계
- 단일 패치, 합성 변위장. 서브메시 투영·hvtrans 연동은 본개발.
