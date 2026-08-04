process IMAGING_INSTANSEG {
    tag "${meta.id}"
    label 'process_medium'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/instanseg:v1-0.0.2'
        : 'quay.io/cellgeni/instanseg:v1-0.0.2'}"

    input:
    tuple val(meta), val(x_min), val(y_min), val(x_max), val(y_max), path(img)

    output:
    tuple val(meta), path("${output_name}"), emit: wkts
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${x_min}_${y_min}_${x_max}_${y_max}"
    output_name = "${prefix}_instanseg_outlines.wkt"
    """
    instanseg_helper.py run \\
        -image-path ${img} \\
        -x-min ${x_min} \\
        -y-min ${y_min} \\
        -x-max ${x_max} \\
        -y-max ${y_max} \\
        -output_name ${output_name} \\
        ${args} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        imaging: \$(instanseg_helper.py version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${x_min}_${y_min}_${x_max}_${y_max}"
    output_name = "${prefix}_instanseg_outlines.wkt"
    """
    touch ${output_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        imaging: \$(instanseg_helper.py version)
    END_VERSIONS
    """
}
