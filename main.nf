include { DECODE_PEAKS_FROM_IMAGE_SERIES ; EXTRACT_AND_DECODE ; SIMPLE_PEAK_COUNTING ; REGISTER_AND_PEAK_COUNTING } from './workflows/main'

params.images = [
    [
        ['id': "run_id"],
        [
            "cycle1",
            "cycle2",
        ],
    ]
]

params.image_stack = [
    [["id": "image_stack_id"], "image_stack.tif"]
]

params.chs_to_call_peaks = [1, 2]

params.codebook = [["id": ''], "codebook.csv", "readouts.csv"]
params.segmentation_method = "CELLPOSE"
params.http_base_url = "http://webatlas.cog.sanger.ac.uk/s3/"


workflow RUN_PEAK_COUNTING_IMAGE_SERIES {
    SIMPLE_PEAK_COUNTING(channel.from(params.images))
}

workflow RUN_PEAK_COUNTING_IMAGE_STACK {
    REGISTER_AND_PEAK_COUNTING(channel.from(params.image_stack))
}

workflow RUN_DECODING_IMAGE_SERIES {
    DECODE_PEAKS_FROM_IMAGE_SERIES(
        channel.from(params.images),
        params.segmentation_method,
        params.chs_to_call_peaks,
        params.codebook,
    )
}

workflow RUN_DECODING_IMAGE_STACK {
    image_stack = channel.from(params.image_stack)
    n_cycle = params.n_cycle_int ? channel.from(params.n_cycle_int) : image_stack.map { it -> [it[0], params.n_cycle_int] }
    EXTRACT_AND_DECODE(
        image_stack,
        params.segmentation_method,
        params.chs_to_call_peaks,
        params.codebook,
        n_cycle,
    )
}

workflow {
    if (params.workflow == "peak_counting_image_series") {
        RUN_PEAK_COUNTING_IMAGE_SERIES()
    }
    else if (params.workflow == "peak_counting_image_stack") {
        RUN_PEAK_COUNTING_IMAGE_STACK()
    }
    else if (params.workflow == "decoding_image_series") {
        RUN_DECODING_IMAGE_SERIES()
    }
    else if (params.workflow == "decoding_image_stack") {
        RUN_DECODING_IMAGE_STACK()
    }
    else {
        println("Invalid workflow specified: ${params.workflow}")
    }
}
