PaSTa is a nextflow-based end-to-end image analysis pipeline for decoding image-based spatial transcriptomics data. It performs imaging cycle registration, cell segmentation and transcripts peak decoding. It is currently supports analysis of three types of ST technology:

- in-situ sequencing-like encoding
- MERFISH-like encoding
- RNAScope-like labelling

Prerequisites:

1. Nextflow. Installation guide: https://www.nextflow.io/docs/latest/getstarted.html
2. Docker or Singularity. Installation guide: https://docs.docker.com/get-docker/ or https://sylabs.io/guides/3.7/user-guide/quick_start.html

## Demo run with GitPod (outdated, but the principal is the same)

Check this HackMD from I2K2024 workshop: https://hackmd.io/w4DeWEDxRlKwIPTDCc77XA

## Basic run

## 1. Clone the repository

```
git clone https://github.com/cellgeni/Image-ST.git
```

## 2. Prepare the run.config file (if needed).

update/overwrite the default running parameters and I/O of the job on HPC/Cloud environment. If local, it will use as much as it can in term of resources. \*

```
process {
        withName: POSTCODE {
                memory = {20.Gb * task.attempt}
                storeDir = "./output/PoSTcode_decoding_output"
        }

        withName: TO_SPATIALDATA {
                memory = {20.Gb * task.attempt}
                ext.args = "--feature_col 'Name' --expansion_in_pixels 30 --save_label_img False"
        }
}
```

## 3. Prepare the parameters file for the workflow (e.g. params.yaml).

Depending on whether your data is pre-registered. You will need two different types of config file:

### 3.1. Stitched, but not registered.

See template in params_templates/image_series_template.yaml

### 3.2. Stitched and registered sequentially (i.e. same cycle order as in the codebook).

See template in params_templates/image_stack_template.yaml

## 4. Run the pipeline

Depending on the config file before you should use different pipeline entries:

### 4.1 for stitched and non-registered run:

```
nextflow run Image-ST -profile <profiles> --workflow decoding_image_series -params-file ./image_series_template.yaml <-c run.config> -resume
```

### 4.2 for stitched and registered run:

```
nextflow run Image-ST -profile <profiles> --workflow decoding_image_stack -params-file ./image_stack_template.yaml <-c run.config> -resume
```

## 5. Check the output in the specified storeDir.

Spin up Napari with napari-spatialdata plugin installed (https://spatialdata.scverse.org/projects/napari/en/latest/notebooks/spatialdata.html)

Then use the following command to visualize the output

```
python -m napari_spatialdata view <path_to_dataset>
```

\*: You may leave the process block empty if you want to use the default parameters.

# FAQ

1. My HOME dir is full when running Singularity image conversion on HPC.

A quick solution is to manually specify singularity dir by setting:

```
singularity cache clean
export SINGULARITY_CACHEDIR=./singularity_image_dir
export NXF_SINGULARITY_CACHEDIR=./singularity_image_dir
```

2. How do I modify parameters to specific process/step?

By following nf-core standard, it is possible to add any parameters to the main script using ext.args=”--[key] [value]” in the run.config file.

An example is

withName: POSTCODE {
ext.args = "--codebook_targer_col L-probe --codebook_code_col code "
}

3. Cannot download pretrained model for the deep-learning tools (Spotiflow/CellPose)

> Exception: URL fetch failure on https://drive.switch.ch/index.php/s/6AoTEgpIAeQMRvX/download: None -- [Errno -3] Temporary failure in name resolution
> Or CellPose
> urllib.error.URLError: <urlopen error [Errno -3] Temporary failure in name resolution>

Mostly likely you've reached max download (?), wait a bit and try later OR manually download those models and update the configuration file.
