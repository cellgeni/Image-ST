#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
from typing import Optional, Sequence

import click
import zarr

from vitessce import (
    CoordinationLevel as CL,
    SpatialDataWrapper,
    VitessceConfig,
    ViewType as vt,
    get_initial_coordination_scope_prefix,
)


def build_config(wrapper: SpatialDataWrapper, description: Optional[str] = None):
    """Construct a VitessceConfig for the given SpatialDataWrapper."""
    vc = VitessceConfig(schema_version="1.0.16", name=description)

    dataset = vc.add_dataset(name="ISS decoding").add_object(wrapper)

    spatial = vc.add_view("spatialBeta", dataset=dataset)
    feature_list = vc.add_view(vt.FEATURE_LIST, dataset=dataset)
    layer_controller = vc.add_view("layerControllerBeta", dataset=dataset)
    [obs_color_encoding_scope] = vc.add_coordination("obsColorEncoding")
    obs_color_encoding_scope.set_value("geneSelection")
    [feature_selection_scope] = vc.add_coordination("featureSelection")
    feature_selection_scope.set_value(None)

    vc.link_views_by_dict(
        [spatial, layer_controller],
        {
            "imageLayer": CL(
                [
                    {
                        "photometricInterpretation": "BlackIsZero",
                        "imageChannel": CL(
                            [
                                {
                                    "spatialTargetC": 0,
                                    "spatialChannelColor": [255, 255, 255],
                                    "spatialChannelWindow": [0, 4000],
                                }
                            ]
                        ),
                    }
                ]
            ),
        },
        scope_prefix=get_initial_coordination_scope_prefix("A", "image"),
    )

    vc.link_views_by_dict(
        [spatial, layer_controller],
        {
            "segmentationLayer": CL(
                [
                    {
                        "segmentationChannel": CL(
                            [
                                {
                                    "obsColorEncoding": "geneSelection",
                                    "spatialChannelOpacity": 0.5,
                                }
                            ]
                        ),
                    }
                ]
            ),
        },
        scope_prefix=get_initial_coordination_scope_prefix("A", "obsSegmentations"),
    )
    obs_sets = vc.add_view(vt.OBS_SETS, dataset=dataset)
    vc.link_views_by_dict(
        [spatial, layer_controller, feature_list, obs_sets],
        {
            "obsType": "cell",
            "featureSelection": feature_selection_scope,
            "obsColorEncoding": obs_color_encoding_scope,
        },
    )

    vc.link_views_by_dict(
        [spatial, layer_controller],
        {
            "segmentationLayer": CL(
                [
                    {
                        "segmentationChannel": CL(
                            [
                                {
                                    "obsColorEncoding": obs_color_encoding_scope,
                                    "featureSelection": feature_selection_scope,
                                }
                            ]
                        ),
                    }
                ]
            ),
        },
        scope_prefix=get_initial_coordination_scope_prefix("A", "obsSegmentations"),
    )
    vc.layout((layer_controller / obs_sets / feature_list) | spatial)
    return vc


@click.command()
@click.option(
    "--sdata-store", required=True, help="Path to the SpatialData zarr store."
)
@click.option(
    "--image-path", required=True, help="Path to the image data within the zarr store."
)
@click.option(
    "--table-path", required=True, help="Path to the table data within the zarr store."
)
@click.option(
    "--obs-feature-matrix-path",
    required=True,
    help="Path to the observation feature matrix within the zarr store.",
)
@click.option(
    "--obs-segmentations-path",
    required=True,
    help="Path to the observation segmentations within the zarr store.",
)
@click.option("--coordinate-system", default="global", help="Coordinate system to use.")
@click.option(
    "--description", default=None, help="Description for the Vitessce configuration."
)
@click.option("--output-file", default=None, help="Output JSON file name.")
@click.option(
    "--data_http_url",
    default="https://webatlas.sanger.ac.uk/s3/",
    help="https base URL for data file.",
)
@click.option(
    "--obs_embedding_paths",
    multiple=True,
    help="Paths to observation embeddings within the zarr store.",
    default=["obsm/X_umap", "obsm/X_pca"],
)
@click.version_option(version='0.0.1', prog_name="spatialdata_vitessce_config")
def cli(
    sdata_store: str,
    image_path: str,
    table_path: str,
    obs_feature_matrix_path: str,
    obs_segmentations_path: str,
    coordinate_system: str,
    description: Optional[str],
    output_file: Optional[str],
    data_http_url: str,
    obs_embedding_paths: Optional[Sequence[str]] = None,
    obs_embedding_names: Optional[Sequence[str]] = None,
):
    """Generate a Vitessce configuration JSON for a SpatialData store."""
    print(f"Loading SpatialData.zip from store: {sdata_store}")
    sdata_store_p = sdata_store
    is_zip = False
    if output_file is None:
        output_stem = sdata_store_p.rstrip("/").split("/")[-1].split(".")[0]
    if sdata_store.endswith(".zip"):
        print("Detected zip store...")
        if not data_http_url.endswith(".zip"):
            data_http_url += ".zip"
        is_zip = True
        sdata_store = zarr.storage.ZipStore(sdata_store)
        output_file = f"{output_stem}_vitessce_config_zipstore.json"
    elif output_file is None:
        output_file = f"{output_stem}_vitessce_config.json"
    else:
        output_file = output_file

    print("Generating Vitessce configuration...")
    wrapper = SpatialDataWrapper(
        sdata_store=sdata_store,
        image_path=image_path,
        table_path=table_path,
        obs_feature_matrix_path=obs_feature_matrix_path,
        obs_segmentations_path=obs_segmentations_path,
        coordinate_system=coordinate_system,
        coordination_values={
            "obsType": "cell",
        },
        obs_set_paths=["tables/table/obs/region"],
        obs_set_names=["Annotation"],
        is_zip=is_zip,
        obs_embedding_paths=obs_embedding_paths,
        obs_embedding_names=obs_embedding_names,
    )

    vc = build_config(wrapper, description)
    vw = vc.widget()
    config = vw.config.export(
        to="files",
        base_url=data_http_url,
        out_dir="./configs",
    )

    # Strip the HMS S3 specific suffix from the exported URL.
    old_url = config["datasets"][0]["files"][0]["url"]
    config["datasets"][0]["files"][0]["url"] = "/".join(old_url.split("/")[:-3])

    print(f"Writing Vitessce configuration to file {output_file}...")
    with open(output_file, "w") as json_file:
        json.dump(config, json_file, indent=4)


if __name__ == "__main__":
    cli()
