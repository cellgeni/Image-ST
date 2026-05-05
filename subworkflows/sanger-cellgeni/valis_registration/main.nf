include { VALIS_REGISTER } from '../../../modules/sanger-cellgeni/valis/register/main'
include { VALIS_WARP     } from '../../../modules/sanger-cellgeni/valis/warp/main'
include { VALIS_MERGE    } from '../../../modules/sanger-cellgeni/valis/merge/main'

workflow VALIS_REGISTRATION {
    take:
    // [ val(meta), path(reference_img), path(moving_imgs) ]
    // reference_img: fixed image; moving_imgs: one or more files aligned to reference
    ch_slides
    // [ val(meta), path(channel_names_json) ] — pass [ meta, [] ] when not used
    ch_channel_names

    main:
    ch_versions = Channel.empty()

    VALIS_REGISTER(ch_slides)
    ch_versions = ch_versions.mix(VALIS_REGISTER.out.versions.first())

    VALIS_WARP(VALIS_REGISTER.out.registrar)
    ch_versions = ch_versions.mix(VALIS_WARP.out.versions.first())

    // Join registrar with optional channel_names_json on meta
    ch_merge_input = VALIS_REGISTER.out.registrar
        .join(ch_channel_names, by: 0)

    VALIS_MERGE(ch_merge_input)
    ch_versions = ch_versions.mix(VALIS_MERGE.out.versions.first())

    emit:
    results_dir       = VALIS_REGISTER.out.results_dir       // channel: [ val(meta), path(results_dir) ]
    registrar         = VALIS_REGISTER.out.registrar         // channel: [ val(meta), path(*_registrar.pickle) ]
    transforms        = VALIS_REGISTER.out.transforms        // channel: [ val(meta), path(*_transforms.json) ]
    registered_slides = VALIS_WARP.out.registered_slides     // channel: [ val(meta), path(warp_dir) ]
    merged            = VALIS_MERGE.out.merged               // channel: [ val(meta), path(*_merged.ome.tiff) ]
    versions          = ch_versions                          // channel: [ versions ]
}
