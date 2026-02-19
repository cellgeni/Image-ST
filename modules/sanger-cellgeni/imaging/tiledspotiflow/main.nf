process IMAGING_TILEDSPOTIFLOW {
    tag "${meta.id}"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? "quay.io/cellgeni/tiled_spotiflow:0.6.4"
        : "quay.io/cellgeni/tiled_spotiflow:0.6.4"}"

    input:
    tuple val(meta), val(x_min), val(y_min), val(x_max), val(y_max), path(image), val(ch_ind)

    output:
    tuple val(meta), val(ch_ind), path("${output_name}"), emit: peaks
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    output_name = "${meta.id}_ch_${ch_ind}_peaks_Y${y_min}_Y${y_max}_X${x_min}_X${x_max}.csv"
    """
    Spotiflow_call_peaks.py run \
        -image_path ${image} \
        -x_min ${x_min} \
        -y_min ${y_min} \
        -x_max ${x_max} \
        -y_max ${y_max} \
        -C ${ch_ind} \
        -output_name "${output_name}" \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        IMAGING: \$(Spotiflow_call_peaks.py version)
    END_VERSIONS
    """

    stub:
    output_name = "${meta.id}_ch_${ch_ind}_peaks_Y${y_min}_Y${y_max}_X${x_min}_X${x_max}.csv"
    """
    touch "${output_name}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        IMAGING: \$(Spotiflow_call_peaks.py version)
    END_VERSIONS
    """
}
