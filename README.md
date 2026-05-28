# PaSTa / Image-ST (old name)

PaSTa is a Nextflow pipeline for image-based spatial transcriptomics analysis. It
supports end-to-end processing from microscopy images to decoded transcripts,
cell segmentation, SpatialData output, and optional Vitessce and OMERO export
artifacts.

The pipeline currently targets:

- in situ sequencing-like encoding
- MERFISH-like encoding
- RNAScope-like labelling or peak counting

Core analysis steps include image registration, tiled segmentation, spot or peak
calling with Spotiflow, peak profile extraction, transcript decoding with
PoSTcode, SpatialData conversion, and optional OMERO table import.

## Requirements

- Nextflow
- Docker or Singularity
- Access to the input images and codebook files referenced in the params file

Install Nextflow from the official documentation:
<https://www.nextflow.io/docs/latest/getstarted.html>

Container runtime installation guides:

- Docker: <https://docs.docker.com/get-docker/>
- Singularity/Apptainer: <https://apptainer.org/docs/>

## Quick Start

Clone the repository:

```bash
git clone https://github.com/cellgeni/Image-ST.git
cd Image-ST
```

Copy one of the parameter templates and edit it for your data:

```bash
cp params_templates/image_series_template.yaml params.yaml
```

Run the pipeline:

```bash
nextflow run . \
  -profile local,docker \
  --workflow decoding_image_series \
  -params-file params.yaml \
  -resume
```

For HPC runs, replace the profile list with the profiles that match your
environment, for example:

```bash
nextflow run . \
  -profile lsf,singularity \
  --workflow decoding_image_stack \
  -params-file params.yaml \
  -c run.config \
  -resume
```

## Workflows

Select the workflow with `--workflow`.

| Workflow | Use case | Main input parameter |
| --- | --- | --- |
| `decoding_image_series` | Separated per-cycle images, one file per imaging cycle, that need registration before decoding | `images` |
| `decoding_image_stack` | Stacked multi-cycle images, with all cycles in one image stack, for decoding | `image_stack` |
| `peak_counting_image_series` | Peak counting from separated image files | `images` |
| `peak_counting_image_stack` | Peak counting from stacked images | `image_stack` |

## Parameter Files

Use the templates in `params_templates/` as starting points:

- `image_series_template.yaml`: separated per-cycle image files
- `image_stack_template.yaml`: stacked multi-cycle image files
- `template_RNAscope.yaml`: RNAScope-like peak counting input

Typical parameters:

| Parameter | Used by | Purpose | Example |
| --- | --- | --- | --- |
| `workflow` | all runs | Selects the pipeline entry point | `decoding_image_series` |
| `images` | image-series workflows | Sample IDs and separated per-cycle image paths | `[[{ id: sample_1 }, ["cycle1.ome.tif", "cycle2.ome.tif"]]]` |
| `image_stack` | image-stack workflows | Sample IDs and stacked multi-cycle image paths | `[[{ id: sample_1 }, "multi_cycle_stack.ome.tif"]]` |
| `codebook` | decoding workflows | Codebook and readout reference files | `[[{ id: sample_1 }, "codebook.xlsx", "./dummy.txt"]]` |
| `chs_to_call_peaks` | peak calling | Channel indexes used for peak calling | `[1, 2]` |
| `segmentation_method` | segmentation | Segmentation module to use | `CELLPOSE` |
| `registration_method` | image-series decoding | Registration method for separated images | `microaligner` or `valis` |
| `n_cycle_int` | image-stack decoding | Number of imaging cycles in each stack | `[[{ id: sample_1 }, 7]]` |
| `cell_diameters` | segmentation | Expected cell diameter values for tiled segmentation | `[30]` |
| `http_base_url` | Vitessce config | Base URL used when generating Vitessce configs | `http://webatlas.cog.sanger.ac.uk/s3` |
| `importsegmentation` | optional OMERO import | OMERO image and table metadata for importing cells and transcripts | `[[{ id: sample_1 }, 12345, "omero.example.org", "spots", "cells", "./omero_zarr"]]` |

The `id` metadata must match across related inputs for the same sample, such as
`images`, `image_stack`, `codebook`, `n_cycle_int`, and `importsegmentation`.

Example image-series params:

```yaml
chs_to_call_peaks: [1]
registration_method: microaligner
segmentation_method: CELLPOSE

codebook:
  - [
      { id: sample_1 },
      "codebook.xlsx",
      "./dummy.txt",
    ]

images:
  - [
      { id: sample_1 },
      [
        "cycle1.ome.tif",
        "cycle2.ome.tif",
        "cycle3.ome.tif",
      ],
    ]
```

Example image-stack params:

```yaml
chs_to_call_peaks: [1]
segmentation_method: CELLPOSE

n_cycle_int:
  - [{ id: sample_1 }, 7]

codebook:
  - [
      { id: sample_1 },
      "codebook.xlsx",
      "./dummy.txt",
    ]

image_stack:
  - [
      { id: sample_1 },
      "multi_cycle_multi_channel.ome.tif",
    ]
```

## Runtime Configuration

Use a Nextflow config file such as `run.config` to override resources, process
arguments, queues, or output locations without editing the pipeline.

Example:

```groovy
process {
  withName: IMAGING_POSTCODE {
    memory = { 20.GB * task.attempt }
    ext.args = "--codebook_target_col Gene --codebook_code_col Code --coding_col_prefix 'Readout_*' --min_prob 0.95"
  }

  withName: TO_SPATIALDATA {
    memory = { 20.GB * task.attempt }
    ext.args = "--feature_col 'Name' --expansion_in_pixels 30 --save_label_img False"
  }
}
```

Then pass the config with `-c run.config`.

Available profiles are defined in `nextflow.config`:

- `local`: local executor
- `lsf`: Sanger LSF profile using `conf/sanger.config`
- `docker`: enable Docker containers
- `singularity`: enable Singularity containers
- `tiger` and `cub`: work-in-progress GPU queue presets; these profiles do not
  currently work

Profiles can be combined, for example `-profile lsf,singularity`.

## Outputs

Pipeline outputs are written under `params.out_dir`, which defaults to
`./output`. Reports are written under `params.report_dir`, which defaults to
`./reports`.

Main outputs include:

- decoded transcript tables from PoSTcode
- tiled segmentation outputs
- SpatialData datasets
- OMERO-compatible cell and transcript tables
- Vitessce configuration files
- provenance files from the `nf-prov` plugin:
  - `bco.json`
  - `ro-crate-metadata.json`

To inspect a SpatialData output in napari, install `napari-spatialdata` and run:

```bash
python -m napari_spatialdata view <path_to_spatialdata_dataset>
```

## Development And Testing

Modules and subworkflows are organised under:

- `modules/local/`
- `modules/sanger-cellgeni/`
- `modules/nf-core/`
- `subworkflows/local/`
- `subworkflows/sanger-cellgeni/`

Many modules include nf-test tests under their `tests/` directories. Run an
individual test with:

```bash
nf-test test <path/to/main.nf.test>
```

## Troubleshooting

### Singularity cache fills `$HOME`

Set a project-local cache directory before running:

```bash
singularity cache clean
export SINGULARITY_CACHEDIR=./singularity_image_dir
export NXF_SINGULARITY_CACHEDIR=./singularity_image_dir
```

### Process-specific parameters need changing

Override module arguments in a Nextflow config file with `ext.args`:

```groovy
process {
  withName: IMAGING_POSTCODE {
    ext.args = "--codebook_target_col Gene --codebook_code_col Code"
  }
}
```

### Pretrained model download fails

Pretrained models for deep-learning tools such as Spotiflow and CellPose should
already be included in the pipeline containers. If a process still tries to
download a model at runtime, or fails because a pretrained model is missing,
please report it as a container or pipeline issue.

## Citation

If you use this repository, please cite it using the metadata in
`CITATION.cff`.

## License

This project is distributed under the MIT license. See `LICENSE`.
