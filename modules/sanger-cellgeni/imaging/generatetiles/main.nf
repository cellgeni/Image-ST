process IMAGING_GENERATETILES {
    tag "${meta.id}"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? "quay.io/cellgeni/imagetileprocessor:0.1.16"
        : "quay.io/cellgeni/imagetileprocessor:0.1.16"}"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("${output_name}"), emit: tile_coords
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    output_name = "${prefix}_tile_coords.csv"
    """
    tile-2d-image run \\
        --image ${image} \\
        --output_name "${output_name}" \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        IMAGING: \$(tile-2d-image version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    output_name = "${prefix}_tile_coords.csv"
    """
    touch "${output_name}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        IMAGING: \$(tile-2d-image version)
    END_VERSIONS
    """
}
