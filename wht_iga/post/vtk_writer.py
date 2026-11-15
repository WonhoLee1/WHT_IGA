"""ParaView-compatible VTKHDF writer and reader for WHT_IGA.

Writes ``<job>.vtkhdf`` following the VTK/ParaView VTKHDF 2.1
``UnstructuredGrid`` specification (including ``Steps`` time-series metadata
for single- or multi-frame output) alongside the custom
``/WHT_IGA/control_points/<frame>/{U, NT11, state}`` group for exact
control-point recovery.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import h5py
import numpy as np

from wht_iga.post.proxy_mesh import ProxyMesh


def _as_frame_list(
    frames_or_mesh: ProxyMesh | Sequence[ProxyMesh],
) -> list[ProxyMesh]:
    if isinstance(frames_or_mesh, ProxyMesh):
        return [frames_or_mesh]
    frames = list(frames_or_mesh)
    if not frames:
        raise ValueError("At least one ProxyMesh frame is required")
    return frames


def write_vtkhdf(
    path: str | Path,
    frames_or_mesh: ProxyMesh | Sequence[ProxyMesh],
    control_points: dict[int, dict[str, np.ndarray]] | None = None,
    times: Sequence[float] | None = None,
    block_names: Sequence[str] | None = None,
) -> str:
    """Write a ParaView-compatible VTKHDF file with ``/WHT_IGA/control_points``.

    Parameters
    ----------
    path:
        Output ``.vtkhdf`` file path.
    frames_or_mesh:
        Single :class:`ProxyMesh` or sequence of per-frame :class:`ProxyMesh`
        instances.
    control_points:
        Optional mapping ``{frame_idx: {var_name: ndarray}}`` written under
        ``/WHT_IGA/control_points/<frame_idx>/<var_name>``.
    times:
        Optional physical time values per frame (defaults to ``0.0, 1.0, ...``).
    block_names:
        Optional list of block names stored under ``/WHT_IGA/block_names``.

    Returns
    -------
    str
        Resolved string path of the written ``.vtkhdf`` file.
    """
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    frames = _as_frame_list(frames_or_mesh)
    n_frames = len(frames)

    if times is None:
        time_vals = np.arange(n_frames, dtype=np.float64)
    else:
        time_vals = np.asarray(list(times), dtype=np.float64)
        if len(time_vals) != n_frames:
            raise ValueError(f"times length {len(time_vals)} != n_frames {n_frames}")

    n_pts_arr = np.array([f.n_points for f in frames], dtype=np.int64)
    n_cells_arr = np.array([f.n_cells for f in frames], dtype=np.int64)
    n_conn_arr = np.array([len(f.connectivity) for f in frames], dtype=np.int64)

    points = np.vstack([f.points for f in frames]) if sum(n_pts_arr) > 0 else np.zeros((0, 3), dtype=np.float64)
    connectivity = (
        np.concatenate([f.connectivity for f in frames])
        if sum(n_conn_arr) > 0
        else np.zeros(0, dtype=np.int64)
    )
    offsets = (
        np.concatenate([f.offsets for f in frames if f.offsets is not None])
        if n_frames > 0
        else np.zeros(1, dtype=np.int64)
    )
    types = (
        np.concatenate([f.cell_types for f in frames])
        if sum(n_cells_arr) > 0
        else np.zeros(0, dtype=np.uint8)
    )

    pt_keys = sorted(set.intersection(*(set(f.point_data.keys()) for f in frames))) if frames else []
    cell_keys = sorted(set.intersection(*(set(f.cell_data.keys()) for f in frames))) if frames else []

    pt_offsets = np.zeros(n_frames, dtype=np.int64)
    cell_offsets = np.zeros(n_frames, dtype=np.int64)
    conn_offsets = np.zeros(n_frames, dtype=np.int64)
    if n_frames > 1:
        np.cumsum(n_pts_arr[:-1], out=pt_offsets[1:])
        np.cumsum(n_cells_arr[:-1], out=cell_offsets[1:])
        np.cumsum(n_conn_arr[:-1], out=conn_offsets[1:])

    with h5py.File(out_path, "w") as h5:
        root = h5.create_group("VTKHDF")
        root.attrs.create("Version", np.array([2, 1], dtype=np.int64))
        root.attrs.create("Type", np.bytes_("UnstructuredGrid"))

        root.create_dataset("NumberOfPoints", data=n_pts_arr, dtype=np.int64)
        root.create_dataset("NumberOfCells", data=n_cells_arr, dtype=np.int64)
        root.create_dataset("NumberOfConnectivityIds", data=n_conn_arr, dtype=np.int64)
        root.create_dataset("Points", data=points, dtype=np.float64)
        root.create_dataset("Connectivity", data=connectivity, dtype=np.int64)
        root.create_dataset("Offsets", data=offsets, dtype=np.int64)
        root.create_dataset("Types", data=types, dtype=np.uint8)

        pd_grp = root.create_group("PointData")
        for key in pt_keys:
            arr = np.concatenate([f.point_data[key] for f in frames], axis=0)
            pd_grp.create_dataset(key, data=arr, dtype=np.float64)

        cd_grp = root.create_group("CellData")
        for key in cell_keys:
            arr = np.concatenate([f.cell_data[key] for f in frames], axis=0)
            cd_grp.create_dataset(key, data=arr, dtype=np.float64)

        steps = root.create_group("Steps")
        steps.attrs.create("NSteps", np.int64(n_frames))
        steps.create_dataset("Values", data=time_vals, dtype=np.float64)
        steps.create_dataset("PartOffset", data=np.arange(n_frames, dtype=np.int64), dtype=np.int64)
        steps.create_dataset("NumberOfParts", data=np.ones(n_frames, dtype=np.int64), dtype=np.int64)
        steps.create_dataset("PointOffsets", data=pt_offsets, dtype=np.int64)
        steps.create_dataset("CellOffsets", data=cell_offsets, dtype=np.int64)
        steps.create_dataset("ConnectivityIdOffsets", data=conn_offsets, dtype=np.int64)

        pd_off_grp = steps.create_group("PointDataOffsets")
        for key in pt_keys:
            pd_off_grp.create_dataset(key, data=pt_offsets, dtype=np.int64)

        cd_off_grp = steps.create_group("CellDataOffsets")
        for key in cell_keys:
            cd_off_grp.create_dataset(key, data=cell_offsets, dtype=np.int64)

        wht_grp = h5.create_group("WHT_IGA")
        if block_names is not None:
            dt = h5py.string_dtype(encoding="utf-8")
            wht_grp.create_dataset("block_names", data=np.asarray(list(block_names), dtype=object), dtype=dt)

        cp_root = wht_grp.create_group("control_points")
        for frame_idx, vars_dict in (control_points or {}).items():
            f_grp = cp_root.create_group(str(int(frame_idx)))
            for var_name, arr in vars_dict.items():
                arr_np = np.asarray(arr)
                f_grp.create_dataset(var_name, data=arr_np)

    return str(out_path)


def read_vtkhdf(path: str | Path) -> dict:
    """Read a VTKHDF file and ``/WHT_IGA`` metadata/control points for inspection.

    Parameters
    ----------
    path:
        Path to ``.vtkhdf`` file.

    Returns
    -------
    dict
        Dictionary containing ``version``, ``type``, ``times``, ``frames``
        (list of :class:`ProxyMesh`), ``block_names``, and ``control_points``.
    """
    in_path = Path(path)
    with h5py.File(in_path, "r") as h5:
        root = h5["VTKHDF"]
        version = np.asarray(root.attrs["Version"], dtype=np.int64)
        raw_type = root.attrs["Type"]
        vtk_type = raw_type.decode("ascii") if isinstance(raw_type, (bytes, np.bytes_)) else str(raw_type)

        n_pts = np.asarray(root["NumberOfPoints"], dtype=np.int64)
        n_cells = np.asarray(root["NumberOfCells"], dtype=np.int64)
        n_conn = np.asarray(root["NumberOfConnectivityIds"], dtype=np.int64)
        all_pts = np.asarray(root["Points"], dtype=np.float64)
        all_conn = np.asarray(root["Connectivity"], dtype=np.int64)
        all_offsets = np.asarray(root["Offsets"], dtype=np.int64)
        all_types = np.asarray(root["Types"], dtype=np.uint8)

        pd_all = {k: np.asarray(v, dtype=np.float64) for k, v in root["PointData"].items()}
        cd_all = {k: np.asarray(v, dtype=np.float64) for k, v in root["CellData"].items()}

        if "Steps" in root and "Values" in root["Steps"]:
            times = np.asarray(root["Steps"]["Values"], dtype=np.float64)
        else:
            times = np.arange(len(n_pts), dtype=np.float64)

        frames: list[ProxyMesh] = []
        p_off = 0
        c_off = 0
        conn_off = 0
        off_idx = 0
        for f_idx in range(len(n_pts)):
            np_f = int(n_pts[f_idx])
            nc_f = int(n_cells[f_idx])
            nconn_f = int(n_conn[f_idx])

            pts_f = all_pts[p_off : p_off + np_f]
            conn_f = all_conn[conn_off : conn_off + nconn_f]
            offsets_f = all_offsets[off_idx : off_idx + nc_f + 1]
            types_f = all_types[c_off : c_off + nc_f]
            pd_f = {k: v[p_off : p_off + np_f] for k, v in pd_all.items()}
            cd_f = {k: v[c_off : c_off + nc_f] for k, v in cd_all.items()}

            frames.append(
                ProxyMesh(
                    points=pts_f,
                    connectivity=conn_f,
                    cell_types=types_f,
                    point_data=pd_f,
                    cell_data=cd_f,
                    offsets=offsets_f,
                )
            )
            p_off += np_f
            c_off += nc_f
            conn_off += nconn_f
            off_idx += nc_f + 1

        control_points: dict[int, dict[str, np.ndarray]] = {}
        block_names: list[str] = []
        if "WHT_IGA" in h5:
            wht = h5["WHT_IGA"]
            if "block_names" in wht:
                block_names = [
                    x.decode("utf-8") if isinstance(x, (bytes, np.bytes_)) else str(x)
                    for x in wht["block_names"][:]
                ]
            if "control_points" in wht:
                for f_key, f_grp in wht["control_points"].items():
                    control_points[int(f_key)] = {k: np.asarray(v) for k, v in f_grp.items()}

    return {
        "version": version,
        "type": vtk_type,
        "times": times,
        "frames": frames,
        "block_names": block_names,
        "control_points": control_points,
    }
