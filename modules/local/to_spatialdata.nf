process TO_SPATIALDATA {
    tag "${meta.id}"
    label 'process_medium'

    container "quay.io/bioinfotongli/spatialdata:0.2.2"
    publishDir params.out_dir + "/spatialdata", mode: 'copy'

    input:
    tuple val(meta), path(transcripts), path(cells), path(registered_image)

    output:
    tuple val(meta), path("${out_name}"), emit: spatialdata
    tuple val(meta), path("${prefix}_count_matrix.csv"), optional: true
    tuple val(meta), path("${prefix}_cell_props.csv"), optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    out_name = "${prefix}.sdata"
    """
    export NUMBA_CACHE_DIR=/tmp/numba_cache
    /opt/conda/bin/python ${workflow.projectDir}/bin/to_spatialdata.py run \\
        --transcripts ${transcripts} \\
        --cells ${cells} \\
        --out_name ${out_name} \\
        --registered_image ${registered_image} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(/opt/conda/bin/python ${workflow.projectDir}/bin/to_spatialdata.py version))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sdata

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(/opt/conda/bin/python ${workflow.projectDir}/bin/to_spatialdata.py version))
    END_VERSIONS
    """
}
