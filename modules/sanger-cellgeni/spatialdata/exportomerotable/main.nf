process SPATIALDATA_EXPORTOMEROTABLE {
    tag "${meta.id}"
    label 'process_single'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/spatialdata-io:0.6.0'
        : 'quay.io/cellgeni/spatialdata-io:0.6.0'}"

    input:
    tuple val(meta), path(bundle)

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.csv"), emit: cells_csv
    tuple val(meta), path("${task.ext.prefix ?: meta.id}_transcripts.csv"), emit: transcripts_csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out_csv = "${prefix}.csv"
    def out_transcripts_csv = "${prefix}_transcripts.csv"
    """
    spatialdata_shapes_to_csv.py \\
        ${bundle} \\
        --out ${out_csv} \\
        --out-transcripts ${out_transcripts_csv} \\
        ${args}

    python - <<'PY' > versions.yml
    import platform
    try:
        import spatialdata
        spatialdata_version = getattr(spatialdata, "__version__", "unknown")
    except Exception:
        spatialdata_version = "unknown"

    print("${task.process}:")
    print(f"  python: {platform.python_version()}")
    print(f"  spatialdata: {spatialdata_version}")
    PY
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.csv
    touch ${prefix}_transcripts.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
        spatialdata: stub
    END_VERSIONS
    """
}
