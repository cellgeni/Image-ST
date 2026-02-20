#!/usr/bin/env/ nextflow

include { REGISTER_AS_SPATIALDATA } from '../subworkflows/local/registration'
include { MICRO_ALIGNER_REGISTRATION } from '../subworkflows/sanger-cellgeni/microaligner_registration/main'
include { TILED_SEGMENTATION } from '../subworkflows/sanger-cellgeni/tiled_segmentation/main'
include { TILED_SPOTIFLOW } from '../subworkflows/sanger-cellgeni/tiled_spotiflow/main'
include { IMAGING_EXTRACTPEAKPROFILE as EXTRACT_PEAK_PROFILE } from '../modules/sanger-cellgeni/imaging/extractpeakprofile/main'
include { IMAGING_POSTCODE as POSTCODE } from '../modules/sanger-cellgeni/imaging/postcode/main'
include { TO_SPATIALDATA } from '../modules/local/to_spatialdata'
include { SPATIAL_GENERATEVITESSCECONFIG } from '../modules/sanger-cellgeni/spatial/generatevitessceconfig/main'


workflow DECODE_PEAKS_FROM_IMAGE_SERIES {
    take:
    images
    segmentation_method
    chs_to_call_peaks
    coding_references

    main:
    n_cycle = images.map { it ->
        [it[0], it[1].size()]
    }
    MICRO_ALIGNER_REGISTRATION(images)
    EXTRACT_AND_DECODE(
        MICRO_ALIGNER_REGISTRATION.out.image,
        segmentation_method,
        chs_to_call_peaks,
        coding_references,
        n_cycle,
    )

    emit:
    spatialdata = EXTRACT_AND_DECODE.out.spatialdata // channel: [ val(meta), [ spatialdata ] ]
}

workflow SIMPLE_PEAK_COUNTING {
    take:
    image_stack

    main:
    TILED_SEGMENTATION(image_stack, params.segmentation_method)
    TILED_SPOTIFLOW(image_stack, params.chs_to_call_peaks)
    TO_SPATIALDATA(
        TILED_SPOTIFLOW.out.spots_csv.combine(TILED_SEGMENTATION.out.geojson, by: 0).combine(image_stack, by: 0)
    )

    emit:
    spatialdata = TO_SPATIALDATA.out.spatialdata // channel: [ val(meta), [ spatialdata ] ]
}

workflow REGISTER_AND_PEAK_COUNTING {
    take:
    images

    main:
    MICRO_ALIGNER_REGISTRATION(images)
    SIMPLE_PEAK_COUNTING(MICRO_ALIGNER_REGISTRATION.out.image)

    emit:
    spatialdata = SIMPLE_PEAK_COUNTING.out.spatialdata // channel: [ val(meta), [ spatialdata ] ]
}

workflow EXTRACT_AND_DECODE {
    take:
    image_stack
    segmentation_method
    chs_to_call_peaks
    coding_references
    n_cycle

    main:
    TILED_SEGMENTATION(image_stack, segmentation_method)
    TILED_SPOTIFLOW(image_stack, channel.from(chs_to_call_peaks))
    // Run the decoding
    EXTRACT_PEAK_PROFILE(image_stack.join(TILED_SPOTIFLOW.out.spots_csv))
    codebook = channel.from(coding_references)
        .map { meta, codebook, readouts ->
            [
                meta,
                file(codebook, checkIfExists: true, type: 'file'),
                file(readouts, checkIfExists: false, type: 'file'),
            ]
        }
    POSTCODE(EXTRACT_PEAK_PROFILE.out.peak_profile.join(codebook).join(n_cycle))
    // Contrsuct the spatial data object
    TO_SPATIALDATA(
        POSTCODE.out.decoded_peaks.combine(TILED_SEGMENTATION.out.geojson, by: 0).combine(image_stack, by: 0)
    )
    SPATIAL_GENERATEVITESSCECONFIG(
        TO_SPATIALDATA.out.spatialdata.map { meta, sdata ->
            def raw_name = "raw_image"
            def label_name = "cell_labels"
            def table_name = "table"
            def http_base_url = params.http_base_url + "/${sdata.name}" ?: "http://webatlas.cog.sanger.ac.uk/s3/${sdata.name}"
            return [meta, sdata, raw_name, label_name, table_name, http_base_url]
        }
    )

    emit:
    spatialdata = TO_SPATIALDATA.out.spatialdata // channel: [ val(meta), [ spatialdata ] ]
    vitessce_config = SPATIAL_GENERATEVITESSCECONFIG.out.vitessce_config // channel: [ val(meta), path(vitessce_config) ]
}
