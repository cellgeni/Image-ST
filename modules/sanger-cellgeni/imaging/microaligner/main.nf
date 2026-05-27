params.debug = false

process IMAGING_MICROALIGNER {
    tag "${meta.id}"
    label 'process_large'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? "quay.io/cellgeni/microaligner:v1.0.7"
        : "quay.io/cellgeni/microaligner:v1.0.7"}"

    input:
    tuple val(meta), path(config), path(images)
    val method

    output:
    tuple val(meta), path("${prefix}_${method}_reg_result_stack.tif"), emit: registered_image
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    microaligner --config ${config} \
        --NumberOfWorkers ${task.cpus} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(microaligner --version 2>&1) | sed 's/^.*microaligner //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_optflow_reg_result_stack.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(microaligner --version 2>&1) | sed 's/^.*microaligner //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
