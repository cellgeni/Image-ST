#!/usr/bin/env nextflow
//
// adopted from https://github.com/cellgeni/spatialxe/blob/dev/subworkflows/local/proseg_preset_proseg2baysor/main.nf
// Runs proseg for the and proseg2baysor to generate cell ploygons
//

include { PROSEG } from './../../modules/nf-core/proseg/proseg/main'
include { PROSEG_TO_BAYSOR } from './../../modules/nf-core/proseg/proseg_to_baysor/main'

workflow PROSEG_PRESET_PROSEG2BAYSOR {
    take:
    ch_transcripts_csv // channel: [ val(meta), [ "transcripts.csv" ] ]

    main:

    ch_versions = channel.empty()

    // run proseg with the bespoke format
    PROSEG(
        ch_transcripts_csv,
        null,
        ["parquet", "csv", "csv"],
    )
    ch_versions = ch_versions.mix(PROSEG.out.versions)

    // run proseg-to-baysor on the data generated with the proseg run
    // PROSEG_TO_BAYSOR(PROSEG.out.transcript_metadata, PROSEG.out.cell_polygons)
    PROSEG_TO_BAYSOR(PROSEG.out.sd_zarr)
    ch_versions = ch_versions.mix(PROSEG_TO_BAYSOR.out.versions)

    emit:
    proseg_sd_zarr = PROSEG.out.sd_zarr // channel: [ val(meta), [ "cell-polygons.geojson.gz" ] ]
    xr_polygons = PROSEG_TO_BAYSOR.out.baysor_cell_polygons // channel: [ val(meta), [ "xr-cell-polygons.geojson" ] ]
    xr_metadata = PROSEG_TO_BAYSOR.out.baysor_transcript_metadata // channel: [ [ "xr-transcript-metadata.csv" ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
