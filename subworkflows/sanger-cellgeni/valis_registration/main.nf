include { VALIS_REGISTER } from '../../../modules/sanger-cellgeni/valis/register/main'
include { VALIS_WARP     } from '../../../modules/sanger-cellgeni/valis/warp/main'
include { VALIS_MERGE    } from '../../../modules/sanger-cellgeni/valis/merge/main'

workflow VALIS_REGISTRATION {
    take:
    ch_slides
    ch_channel_names

    main:
    ch_versions = channel.empty()

    // Split the image list: first element → reference, rest → moving
    ch_register_input = ch_slides.map { meta, images ->
        def imgs = images instanceof List ? images : [images]
        def ref = imgs[0]
        def moving = imgs.size() > 1 ? imgs[1..-1] : imgs
        tuple(meta, ref, moving)
    }

    VALIS_REGISTER(ch_register_input)
    ch_versions = ch_versions.mix(VALIS_REGISTER.out.versions.first())

    VALIS_WARP(VALIS_REGISTER.out.results_dir)
    ch_versions = ch_versions.mix(VALIS_WARP.out.versions.first())

    // MERGE needs the slide files from results_dir and must run after WARP
    // (VALIS updates the registrar pickle during warp).
    ch_merge_input = VALIS_WARP.out.registered_slides
        .join(VALIS_REGISTER.out.results_dir, by: 0)
        .join(ch_channel_names, by: 0)
        .map { meta, _slides_dir, results_dir, channel_names -> tuple(meta, results_dir, channel_names) }

    VALIS_MERGE(ch_merge_input)
    ch_versions = ch_versions.mix(VALIS_MERGE.out.versions.first())

    emit:
    results_dir       = VALIS_REGISTER.out.results_dir // channel: [ val(meta), path(results_dir) ]
    registrar         = VALIS_REGISTER.out.registrar // channel: [ val(meta), path(*_registrar.pickle) ]
    transforms        = VALIS_REGISTER.out.transforms // channel: [ val(meta), path(*_transforms.json) ]
    registered_slides = VALIS_WARP.out.registered_slides // channel: [ val(meta), path(warp_dir) ]
    merged            = VALIS_MERGE.out.merged // channel: [ val(meta), path(*_merged.ome.tiff) ]
    versions          = ch_versions // channel: [ versions ]
}
