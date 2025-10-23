PaSTa is a nextflow-based end-to-end image analysis pipeline for decoding image-based spatial transcriptomics data. It performs imaging cycle registration, cell segmentation and transcripts peak decoding. It is currently supports analysis of three types of ST technology:

- in-situ sequencing-like encoding
- MERFISH-like encoding
- RNAScope-like labelling

Prerequisites:
1. Nextflow. Installation guide: https://www.nextflow.io/docs/latest/getstarted.html
2. Docker or Singularity. Installation guide: https://docs.docker.com/get-docker/ or https://sylabs.io/guides/3.7/user-guide/quick_start.html

Demo run with GitPod
---

Check this HackMD from I2K2024 workshop: https://hackmd.io/w4DeWEDxRlKwIPTDCc77XA

Basic run
---
## 1. Clone the repository
```
git clone https://github.com/cellgeni/Image-ST.git
```
## 2. Prepare the run.config file. 
which specifies the running parameters and I/O of the job on HPC/Cloud environment. If local, it will use as much as it can in term of resources. *
```
process {
        withName: CELLPOSE {
                ext.args = "--channels [0,0]"
                storeDir = "./output/naive_cellpose_segmentation/"
        }

        withName: POSTCODE {
                memory = {20.Gb * task.attempt}
                storeDir = "./output/PoSTcode_decoding_output"
        }

        withName: TO_SPATIALDATA {
                memory = {20.Gb * task.attempt}
                ext.args = "--feature_col 'Name' --expansion_in_pixels 30 --save_label_img False"
        }

        withName: MERGE_OUTLINES {
                storeDir = "./output/merged_cellpose_segmentation/"
        }

        withName: BIOINFOTONGLI_MICROALIGNER {
                memory = {50.Gb * task.attempt}
                storeDir = "./output/registered_stacks"
        }

        withName: BIOINFOTONGLI_TILEDSPOTIFLOW {
                memory = {30.Gb * task.attempt}
                storeDir = "./output/spotiflow_peaks/"
        }

        withName: BIOINFOTONGLI_MERGEPEAKS {
                memory = {50.Gb * task.attempt}
                storeDir = "./output/spotiflow_peaks/"
        }

        withName: BIOINFOTONGLI_CONCATENATEWKTS {
                memory = {50.Gb * task.attempt}
                storeDir = "./output/spotiflow_peaks/"
        }

        withName: EXTRACT_PEAK_PROFILE {
                memory = {50.Gb * task.attempt}
                storeDir = "./output/peak_profiles/"
        }
}
```
## 3. Prepare the parameters file for the workflow (e.g. iss.yaml).
Depending on whether your data is pre-registered. You will need two different types of config file:   
### 3.1. Stitched, but not registered.
```
images:
   - ['id': "test",
       [
         "cycle1.ome.tiff",
         "cycle2.ome.tiff",
         "cycle3.ome.tiff",
         "cycle4.ome.tiff",
         "cycle5.ome.tiff",
         "cycle6.ome.tiff",
       ]
     ]
cell_diameters: [30]
chs_to_call_peaks: [1,2] // channels to call peaks, can be multiple
codebook:
  - ['id': "test", "./codebook.csv", "./dummy.txt"] // has to match the meta in `images` variable
segmentation_method: "CELLPOSE" // or DEEPCELL or STARDIST or INSTANSEG

out_dir: "./output"
```    
### 3.2. Stitched and registered sequentially (i.e. same cycle order as in the codebook)
```
cell_diameters: [30]
chs_to_call_peaks: [27]
codebook:
	- [
	      id: A02,
	      "./codebook.csv",
	      "./dummmy.txt",
	]
image_stack:
	- [
	    id: A02,
	    "my-stitched-and-registered-hyper-stack.ome.tif",
	]
n_cycle_int:
  - [id: A02, 6] # crucial for the decoding
```     
## 4. Run the pipeline
Depending on the config file before you should use different pipeline entries:
### 4.1 for stitched and non-registered run:
```
nextflow run Image-ST/main.nf -profile lsf,singularity -entry RUN_DECODING_IMAGE_SERIES -params-file ./iss.yaml -c run.config -resume
```
### 4.2 for stitched and registered run:
```
nextflow run actions/Image-ST/main.nf -profile lsf,singularity -entry RUN_DECODING_IMAGE_STACK -params-file ./manifest.yaml -c run.config -resume
```

## 5. Check the output in the specified storeDir.

Spin up Napari with napari-spatialdata plugin installed (https://spatialdata.scverse.org/projects/napari/en/latest/notebooks/spatialdata.html)

Then use the following command to visualize the output
```
from napari_spatialdata import Interactive
import spatialdata as spd

data = spd.read_zarr([path-to-.sdata-folder])
Interactive(data)
```

*: You may leave the process block empty if you want to use the default parameters.

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
Or CellPose
urllib.error.URLError: <urlopen error [Errno -3] Temporary failure in name resolution>

Mostly likely you've reached max download (?), wait a bit and try later OR manually download those models and update the configuration file.
