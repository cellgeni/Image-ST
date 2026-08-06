include { POSTCODE_DECODE      } from '../../../modules/sanger-cellgeni/postcode/decode/main'
include { POSTCODE_PREPROCESS  } from '../../../modules/sanger-cellgeni/postcode/preprocess/main'
include { POSTCODE_POSTPROCESS } from '../../../modules/sanger-cellgeni/postcode/postprocess/main'

workflow POSTCODE_DECODING {
    take:
    ch_input // channel: [ val(meta), path(profile), path(tabular_codebook), path(readout_file), val(R), path(spot_loc) ]

    main:

    def ch_versions = channel.empty()

    def ch_for_preprocess = ch_input
        .map { meta, profile, tabular_codebook, readout_file, R, _spot_loc -> [ meta, profile, tabular_codebook, readout_file, R ] }
    def ch_loc = ch_input
        .map { meta, _profile, _tabular_codebook, _readout_file, _R, spot_loc -> [ meta, spot_loc ] }

    POSTCODE_PREPROCESS(ch_for_preprocess)
    ch_versions = ch_versions.mix(POSTCODE_PREPROCESS.out.versions.first())

    POSTCODE_DECODE(POSTCODE_PREPROCESS.out.for_decoding.combine(ch_loc, by: 0))
    ch_versions = ch_versions.mix(POSTCODE_DECODE.out.versions.first())

    POSTCODE_POSTPROCESS(POSTCODE_DECODE.out.model_params_and_losses.combine(ch_loc, by: 0))
    ch_versions = ch_versions.mix(POSTCODE_POSTPROCESS.out.versions.first())

    emit:
    processed_profiles = POSTCODE_PREPROCESS.out.for_decoding // channel: [ val(meta), [ pixel/cell ] ]
    decoded_profiles   = POSTCODE_POSTPROCESS.out.decoded_profile // channel: [ val(meta), [ decoded_profile ] ]
    versions           = ch_versions // channel: [ versions.yml ]
}
