#!/usr/bin/env python3
"""
Load SpatialData and export:
1) cell polygons CSV
2) transcript points CSV (as WKT POINT in the `polygon` column)
"""

import argparse
import os
import sys
from pathlib import Path


def _resolve_reader():
    errors = []

    try:
        os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")
        Path(os.environ["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)

        import spatialdata  # type: ignore

        for attr in ("read_zarr", "read"):
            reader = getattr(spatialdata, attr, None)
            if callable(reader):
                return reader

        spatial_data_cls = getattr(spatialdata, "SpatialData", None)
        reader = getattr(spatial_data_cls, "read", None)
        if callable(reader):
            return reader

        errors.append(
            "spatialdata imported but has no callable read_zarr, read, or SpatialData.read"
        )
    except Exception as exc:
        errors.append(f"import spatialdata failed: {type(exc).__name__}: {exc}")

    raise ImportError(
        "Could not find a SpatialData reader (read_zarr / read / SpatialData.read). "
        + "; ".join(errors)
    )


def _pick_table(sdata, table_name: str | None):
    if not hasattr(sdata, "tables") or sdata.tables is None or len(sdata.tables) == 0:
        raise ValueError("No tables found in SpatialData object.")
    if table_name:
        if table_name not in sdata.tables:
            raise ValueError(
                f"Table '{table_name}' not found. Available: {list(sdata.tables.keys())}"
            )
        return sdata.tables[table_name]
    return next(iter(sdata.tables.values()))


def _pick_shapes(sdata, layer_name: str | None):
    if not hasattr(sdata, "shapes") or sdata.shapes is None or len(sdata.shapes) == 0:
        raise ValueError("No shapes found in SpatialData object.")
    if layer_name:
        if layer_name not in sdata.shapes:
            raise ValueError(
                f"Shapes layer '{layer_name}' not found. Available: {list(sdata.shapes.keys())}"
            )
        return sdata.shapes[layer_name]
    return next(iter(sdata.shapes.values()))


def _pick_points(sdata, layer_name: str | None):
    if not hasattr(sdata, "points") or sdata.points is None or len(sdata.points) == 0:
        raise ValueError("No points found in SpatialData object.")
    if layer_name:
        if layer_name not in sdata.points:
            raise ValueError(
                f"Points layer '{layer_name}' not found. Available: {list(sdata.points.keys())}"
            )
        return sdata.points[layer_name]
    return next(iter(sdata.points.values()))


def _coerce_int(series):
    try:
        return series.astype(int)
    except Exception:
        return series.astype(str).str.extract(r"(\d+)", expand=False).astype(int)


def _to_pandas_frame(obj):
    if hasattr(obj, "compute"):
        return obj.compute()
    return obj


def _string_join_key(series):
    return series.astype(str)


def _shape_ids(gdf, id_column: str | None):
    import pandas as pd

    if id_column:
        if id_column not in gdf.columns:
            raise ValueError(
                f"id column '{id_column}' not found in shapes. Available: {list(gdf.columns)}"
            )
        return pd.Series(gdf[id_column])
    if "object" in gdf.columns:
        return pd.Series(gdf["object"])
    return pd.Series(gdf.index, index=gdf.index)


def _transcripts_to_wkt(points_obj, buffer_radius: float = 0.0):
    import pandas as pd

    df = _to_pandas_frame(points_obj)
    if hasattr(df, "to_pandas"):
        df = df.to_pandas()

    features_df = df.copy().reset_index(drop=False)
    features_df = features_df.rename(columns={"index": "points_index"})

    if "geometry" in df.columns:
        geometry = df["geometry"]
        if buffer_radius > 0:
            geometry = geometry.buffer(buffer_radius)
        try:
            wkt = geometry.to_wkt()
        except Exception:
            wkt = geometry.apply(
                lambda geom: geom.wkt if geom is not None else None
            )
        if "geometry" in features_df.columns:
            features_df = features_df.drop(columns=["geometry"])
        out = pd.concat(
            [
                pd.Series(range(1, len(wkt) + 1), name="object"),
                pd.Series(wkt.values, name="polygon"),
                features_df.reset_index(drop=True),
            ],
            axis=1,
        )
        return out

    x_col = None
    y_col = None
    candidates = [("x", "y"), ("global_x", "global_y"), ("center_x", "center_y")]
    for cx, cy in candidates:
        if cx in df.columns and cy in df.columns:
            x_col, y_col = cx, cy
            break
    if x_col is None or y_col is None:
        raise ValueError(
            f"Could not find transcript coordinate columns. Available columns: {list(df.columns)}"
        )

    if buffer_radius > 0:
        from shapely.geometry import Point

        wkt = df.apply(
            lambda row: Point(row[x_col], row[y_col]).buffer(buffer_radius).wkt,
            axis=1,
        )
    else:
        wkt = "POINT (" + df[x_col].astype(str) + " " + df[y_col].astype(str) + ")"
    out = pd.concat(
        [
            pd.Series(range(1, len(df) + 1), name="object"),
            pd.Series(wkt.values, name="polygon"),
            features_df.reset_index(drop=True),
        ],
        axis=1,
    )
    return out


def main():
    parser = argparse.ArgumentParser(
        description="Export SpatialData cells CSV and transcripts CSV."
    )
    parser.add_argument("spatialdata_path", help="Path to SpatialData zarr store")
    parser.add_argument("--shapes-layer", dest="shapes_layer", default=None)
    parser.add_argument("--points-layer", dest="points_layer", default=None)
    parser.add_argument("--table", dest="table", default=None)
    parser.add_argument(
        "--id-column",
        dest="id_column",
        default=None,
        help="Column in shapes to use as cell id (otherwise uses shapes object column or index)",
    )
    parser.add_argument("--out", default=None, help="Cells CSV output path")
    parser.add_argument(
        "--out-transcripts",
        dest="out_transcripts",
        default=None,
        help="Transcripts CSV output path",
    )
    parser.add_argument(
        "--transcript-buffer",
        dest="transcript_buffer",
        type=float,
        default=0.0,
        help="Buffer radius for transcript points. Default 0 keeps POINT WKT; positive values export circular POLYGON WKT.",
    )
    args = parser.parse_args()

    if args.transcript_buffer < 0:
        raise ValueError("--transcript-buffer must be greater than or equal to 0.")

    stem = Path(args.spatialdata_path).resolve().name
    if stem.endswith(".zarr"):
        stem = stem[: -len(".zarr")]
    if stem == "":
        stem = "spatialdata"

    if args.out is None:
        args.out = f"{stem}.csv"
    if args.out_transcripts is None:
        args.out_transcripts = f"{stem}_transcripts.csv"

    read_sdata = _resolve_reader()
    sdata = read_sdata(args.spatialdata_path)

    gdf = _pick_shapes(sdata, args.shapes_layer)
    points = _pick_points(sdata, args.points_layer)
    table = _pick_table(sdata, args.table)

    if not hasattr(gdf, "geometry"):
        raise TypeError("Shapes layer is not a GeoDataFrame with geometry column.")

    try:
        wkt = gdf.geometry.to_wkt()
    except Exception:
        wkt = gdf.geometry.apply(lambda geom: geom.wkt if geom is not None else None)

    import pandas as pd

    ids = _shape_ids(gdf, args.id_column)
    n = len(wkt)

    obs_df = table.obs.copy()
    obs_df["obs_index"] = obs_df.index
    obs_df["_shape_join_key"] = _string_join_key(
        pd.Series(obs_df.index, index=obs_df.index)
    )
    shape_index = pd.DataFrame(
        {
            "_shape_join_key": _string_join_key(ids).to_numpy(),
            "_shape_order": range(n),
        }
    )
    obs_df = shape_index.merge(obs_df, on="_shape_join_key", how="left").drop(
        columns=["_shape_join_key"]
    ).sort_values("_shape_order")
    obs_df = obs_df.drop(columns=["_shape_order"]).set_index(ids.index)

    unmatched_obs_rows = obs_df["obs_index"].isna().sum()
    if len(table.obs) != n or unmatched_obs_rows:
        print(
            f"WARNING: length mismatch shapes={len(wkt)} obs={len(table.obs)} ids={len(ids)}; "
            f"exporting all shape rows and leaving table columns empty for {unmatched_obs_rows} unmatched shapes.",
            file=sys.stderr,
        )

    base_ids = _coerce_int(ids)
    if base_ids.duplicated().any() or base_ids.isna().any():
        object_ids = pd.Series(range(1, n + 1), index=ids.index, name="object")
    else:
        object_ids = base_ids.rename("object")

    out_df = pd.concat(
        [object_ids, wkt.rename("polygon"), obs_df], axis=1
    )
    out_df.to_csv(args.out, index=False)

    transcripts_df = _transcripts_to_wkt(points, args.transcript_buffer)
    transcripts_df.to_csv(args.out_transcripts, index=False)

    print(f"Wrote {len(out_df)} rows to {args.out}")
    print(f"Wrote {len(transcripts_df)} rows to {args.out_transcripts}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
