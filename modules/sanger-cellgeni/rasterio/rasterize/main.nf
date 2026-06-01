process RASTERIO_RASTERIZE {
    tag "${meta.id}"
    label 'process_medium'

    // conda "${moduleDir}/environment.yml"
    container "quay.io/cellgeni/rasterio:1.5.0"

    input:
    tuple val(meta), path(geojson), path(reference_image)

    output:
    tuple val(meta), path("${prefix}_labels.tif"), emit: label_image
    tuple val("${task.process}"), val('rasterio'), eval('rio --version | sed "s/rio, version //"'), topic: versions, emit: versions_rasterio_rasterize

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    rio rasterize \\
        ${geojson} \\
        --like ${reference_image} \\
        ${args} \\
        -o ./${prefix}_labels.tif
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_labels.tif
    """
}
