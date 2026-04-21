#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2024 Wellcome Sanger Institute

import fire
from dask_image.imread import imread
from spatialdata.models import Image2DModel
from spatialdata import SpatialData, rasterize
from spatialdata.models import ShapesModel, TableModel, PointsModel, Labels2DModel
from geopandas import GeoDataFrame
from xarray import DataArray
from spatialdata.transformations.transformations import Identity
import pandas as pd
import anndata
from shapely import from_wkt, from_geojson, MultiPoint, MultiPolygon
from shapely.geometry import Polygon
from shapely.geometry.collection import GeometryCollection
import numpy as np
from skimage.segmentation import expand_labels
from skimage.measure import regionprops_table, regionprops, find_contours

from collections.abc import Mapping
from types import MappingProxyType
from typing import Any
import logging

logging.basicConfig(level=logging.INFO)

logger = logging.getLogger(__name__)

VERSION = "0.1.0"


def load_shapemodel(cells: str):
    with open(cells, "r") as f:
        if cells.endswith(".wkt"):
            multipoly = from_wkt(f.read())
        elif cells.endswith(".geojson"):
            multipoly = from_geojson(f.read())
    if isinstance(multipoly, MultiPolygon) or isinstance(multipoly, GeometryCollection):
        ids, polys = zip(
            *{i + 1: poly for i, poly in enumerate(list(multipoly.geoms))}.items()
        )
        df = pd.DataFrame({"instance_id": ids, "geometry": polys})
        df["instance_id"] = df["instance_id"].astype(int)
        df.index.name = None
        return ShapesModel.parse(GeoDataFrame(df))
    else:
        raise ValueError("Please provide a wkt file with multipolygon geometry")


def main(
    transcripts: str,
    cells: str,
    registered_image: str,
    out_name: str,
    pixelsize: int = 1,
    y_col: str = "y_int",
    x_col: str = "x_int",
    feature_col: str = "feature_name",
    multiscale_image: bool = True,
    raw_image_channels_to_save: list = [0],
    expansion_in_pixels: int = -1,
    image_models_kwargs: Mapping[str, Any] = MappingProxyType({}),
    imread_kwargs: Mapping[str, Any] = MappingProxyType({}),
    save_label_img: bool = False,
    cell_props: list = ["label", "bbox"],
    # cell_props:list = ['label', 'area', 'intensity_mean', "centroid", "axis_major_length", "axis_minor_length"]
):

    raw_image = imread(registered_image, **imread_kwargs)[raw_image_channels_to_save]
    raw_image = DataArray(raw_image, dims=("c", "y", "x"))

    raw_image_parsed = Image2DModel.parse(
        raw_image,
        scale_factors=[2, 2, 2, 2] if multiscale_image else None,
        transformations={"global": Identity()},
        **image_models_kwargs,
    )

    if transcripts.endswith(".csv"):
        spots = pd.read_csv(transcripts, header=0, sep=",")
    elif transcripts.endswith(".tsv"):
        spots = pd.read_csv(transcripts, header=0, sep="\t")
    elif transcripts.endswith(".wkt"):
        # Assuming that the wkt file contains a multipoint geometry
        with open(transcripts, "r") as f:
            multispots = from_wkt(f.read())
        if not isinstance(multispots, MultiPoint):
            raise ValueError("Please provide a wkt file with multipoint geometry")
        spots = pd.DataFrame(
            [(geom.y, geom.x) for geom in multispots.geoms], columns=[y_col, x_col]
        )
        spots[feature_col] = "spot"
    else:
        raise ValueError("Format not recognized. Please provide a csv, tsv or wkt file")

    required_columns = {y_col, x_col, feature_col}
    missing_columns = required_columns.difference(spots.columns)
    if missing_columns:
        raise ValueError(
            f"Missing required columns in transcripts file: {sorted(missing_columns)}"
        )

    logger.info("Load cell polygons from file")
    cell_shape = load_shapemodel(cells)

    sdata = SpatialData(
        shapes={"cell_shapes": cell_shape},
    )
    logger.info("Rasterizing cell shapes")
    cell_labels = rasterize(
        sdata["cell_shapes"],
        ["y", "x"],
        min_coordinate=[0, 0],
        max_coordinate=[raw_image.shape[-2], raw_image.shape[-1]],
        target_coordinate_system="global",
        target_unit_to_pixels=1.0,
    )
    sdata["raw_image"] = raw_image_parsed

    logger.info("Assigning spots to cells")
    if expansion_in_pixels > 0:
        lab_img = expand_labels(np.array(cell_labels.data), expansion_in_pixels)
        mask = np.squeeze(lab_img).astype(np.int32)
        expanded_polys = {}
        for prop in regionprops(mask):
            r0, c0, r1, c1 = prop.bbox
            r0p = max(0, r0 - 1); c0p = max(0, c0 - 1)
            r1p = min(mask.shape[0], r1 + 1); c1p = min(mask.shape[1], c1 + 1)
            contours = find_contours((mask[r0p:r1p, c0p:c1p] == prop.label).astype(np.uint8), level=0.5)
            if not contours:
                continue
            contour = max(contours, key=len)
            expanded_polys[prop.label] = Polygon([(c[1] + c0p, c[0] + r0p) for c in contour])
        ids = sorted(expanded_polys.keys())
        expanded_df = GeoDataFrame({"instance_id": ids, "geometry": [expanded_polys[i] for i in ids]})
        sdata["cell_shapes_expanded"] = ShapesModel.parse(expanded_df)
    else:
        lab_img = np.array(cell_labels.data)
    props_dict = regionprops_table(
        np.squeeze(lab_img).astype(np.int32),
        intensity_image=np.array(raw_image).transpose(1, 2, 0),
        properties=cell_props,
    )
    props_df = pd.DataFrame(props_dict)
    props_df = props_df[props_df["label"] != 0]
    props_df.to_csv(out_name.replace(".sdata", "_cell_props.csv"))

    if save_label_img:
        sdata["cell_labels"] = Labels2DModel.parse(
            np.squeeze(lab_img).astype(np.int32),
            scale_factors=[2, 2, 2, 2],
            transformations={"global": Identity()},
        )

    cell_ids = lab_img[
        0,
        (spots[y_col] / pixelsize).astype(int),
        (spots[x_col] / pixelsize).astype(int),
    ]
    spots["cell_id"] = cell_ids.astype(int)

    logger.info("Create count matrix")
    count_matrix = spots.pivot_table(
        index="cell_id", columns=feature_col, aggfunc="size", fill_value=0
    )
    count_matrix = count_matrix.drop(count_matrix[count_matrix.index == 0].index)
    count_matrix["num_cells"] = np.max(lab_img)
    count_matrix["num_spots"] = spots.shape[0]
    count_matrix.to_csv(out_name.replace(".sdata", "_count_matrix.csv"))

    props_df_intersect = props_df[props_df["label"].isin(count_matrix.index)]

    logger.info("Construct anndata object")
    adata = anndata.AnnData(X=count_matrix.values, obs=props_df_intersect)

    REGION = "cell_labels"
    REGION_KEY = "region"
    instance_key = "instance_id"
    adata.obs[instance_key] = adata.obs.index.astype(int)
    adata.var_names = count_matrix.columns
    adata.var_names_make_unique()
    adata.obs[REGION_KEY] = REGION

    table = TableModel.parse(
        adata, region=REGION, region_key=REGION_KEY, instance_key=instance_key
    )

    points = PointsModel.parse(
        spots,
        coordinates={"x": x_col, "y": y_col},
        feature_key=feature_col,
        # instance_key=instance_key,
        transformations={"global": Identity()},
    )
    sdata["transcripts"] = points
    sdata["table"] = table
    sdata.write(out_name)


if __name__ == "__main__":
    options = {"run": main, "version": VERSION}
    fire.Fire(options)
