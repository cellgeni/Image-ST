process IMAGING_MERGEOUTLINES {
    tag "${meta.id}"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/imagetileprocessor:0.1.9'
        : 'quay.io/cellgeni/imagetileprocessor:0.1.9'}"

    input:
    tuple val(meta), path(outlines)

    output:
    tuple val(meta), path("${prefix}.wkt"), emit: multipoly_wkts, optional: true
    tuple val(meta), path("${prefix}.geojson"), emit: multipoly_geojsons
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_merged"
    """
    merge-polygons \\
        --wkts ${outlines} \\
        --output_prefix "${prefix}" \\
        ${args} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mergeoutlines: \$(merge-polygons --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_merged"
    """
    touch ${prefix}.geojson

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mergeoutlines: \$(merge-polygons --version)
    END_VERSIONS
    """
}
