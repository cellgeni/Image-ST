include { POSTCODE_DECODE      } from '../../../modules/sanger-cellgeni/postcode/decode/main'
include { POSTCODE_PREPROCESS  } from '../../../modules/sanger-cellgeni/postcode/preprocess/main'
include { POSTCODE_POSTPROCESS } from '../../../modules/sanger-cellgeni/postcode/postprocess/main'

workflow POSTCODE_DECODING {
    take:
    ch_profile // channel: [ val(meta), [ profile ] ]
    ch_tabular_codebook // channel: [ val(meta), [ tabular_codebook ] ]
    ch_readout_file // channel: [ val(meta), [ readout_file ] ]
    ch_R // channel: [ val(meta), [ R ] ]
    ch_loc // channel: [ val(meta), [ spot_loc ] ]

    main:

    def ch_versions = channel.empty()

    def for_decoding = ch_profile
        .combine(
            ch_tabular_codebook,
            by: 0
        )
        .combine(
            ch_readout_file,
            by: 0
        )
        .combine(ch_R, by: 0)

    // for_decoding.view()

    POSTCODE_PREPROCESS(for_decoding)
    ch_versions = ch_versions.mix(POSTCODE_PREPROCESS.out.versions.first())
    // POSTCODE_PREPROCESS.out.for_decoding.combine(ch_loc, by: 0).view()

    POSTCODE_DECODE(POSTCODE_PREPROCESS.out.for_decoding.combine(ch_loc, by: 0))
    ch_versions = ch_versions.mix(POSTCODE_DECODE.out.versions.first())

    POSTCODE_POSTPROCESS(POSTCODE_DECODE.out.model_params_and_losses.combine(ch_loc, by: 0))
    ch_versions = ch_versions.mix(POSTCODE_POSTPROCESS.out.versions.first())

    emit:
    processed_profiles = POSTCODE_PREPROCESS.out.for_decoding // channel: [ val(meta), [ pixel/cell ] ]
    decoded_profiles   = POSTCODE_POSTPROCESS.out.decoded_profile // channel: [ val(meta), [ decoded_profile ] ]
    versions           = ch_versions // channel: [ versions.yml ]
}
