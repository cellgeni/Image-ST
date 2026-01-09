process SPATIAL_GENERATEVITESSCECONFIG {
    tag "${meta.id}"
    label 'process_single'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'vitessce_config:latest'
        : 'vitessce_config:latest'}"

    input:
    tuple val(meta), path(sdata), val(raw_name), val(label_name), val(table_name), val(http_base_url)

    output:
    tuple val(meta), path("${out_json}"), emit: vitessce_config
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_json = "${prefix}_vitessce_config.json"
    """
    spatialdata_vitessce_config.py \\
        --sdata-store ${sdata} \\
        --image-path images/${raw_name} \\
        --table-path tables/${table_name} \\
        --obs-feature-matrix-path tables/${table_name}/X \\
        --obs-segmentations-path labels/${label_name} \\
        --data_http_url ${http_base_url} \\
        --output-file ${out_json} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spatial: \$(spatialdata_vitessce_config.py --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_json = "${prefix}_vitessce_config.json"
    """
    echo ${args}
    
    touch ${out_json}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spatial: \$(spatialdata_vitessce_config.py --version)
    END_VERSIONS
    """
}
