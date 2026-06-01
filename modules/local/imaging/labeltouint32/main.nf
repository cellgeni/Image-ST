process IMAGING_LABELTOUINT32 {
    tag "${meta.id}"
    label 'process_low'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/extract_peak_profile:0.1.0'
        : 'quay.io/cellgeni/extract_peak_profile:0.1.0'}"

    input:
    tuple val(meta), path(label_image)

    output:
    tuple val(meta), path("${prefix}_uint32_labels.npy"), emit: label_image
    tuple val("${task.process}"), val('numpy'), eval('python -c "import numpy; print(numpy.__version__)"'), topic: versions, emit: versions_labeltouint32

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    label_to_uint32.py \\
        --input ${label_image} \\
        --output ${prefix}_uint32_labels.npy \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_uint32_labels.npy
    """
}
