process OMERO_IMPORTSEGMENTATION {
    tag "${meta.id}"
    label 'process_single'
    cpus 30

    // conda "${moduleDir}/environment.yml"
    container "quay.io/cellgeni/roi_convert_ngff:0.6.2"
    secret 'OMERO_USER'
    secret 'OMERO_PASSWORD'

    input:
    tuple val(meta), path(csv), val(image_id), val(host), val(table_name), val(roi_name), val(out_dir)

    output:
    tuple val(meta), path("*.importsegmentation.done.txt"), emit: done
    tuple val("${task.process}"), val('ROI_Converter_NGFF'), val('0.6.2'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ROI_CONVERTER_LOG=${prefix}.ROI_Converter_NGFF.log

    ROI_Converter_NGFF \\
        -r ${image_id} \\
        --server ${host} \\
        --user \$OMERO_USER \\
        --password \$OMERO_PASSWORD \\
        --table_name "${table_name}" \\
        -i ${csv} \\
        --max_procs ${task.cpus} \\
        --name "${roi_name}" \\
        --table \\
        -d "${out_dir}" \\
        ${args} \\
        2>&1 | tee \${ROI_CONVERTER_LOG}

    roi_emitted=\$(grep -Eio 'roi([ _-]?id)?[^0-9]*[0-9]+' \${ROI_CONVERTER_LOG} | grep -Eo '[0-9]+' | tail -n 1 || true)
    if [[ -z "\${roi_emitted}" ]]; then
        roi_emitted="unknown"
    fi
    printf "%s\\n" "\${roi_emitted}" > ${prefix}.importsegmentation.done.txt
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    printf "%s\\n" "unknown" > ${prefix}.importsegmentation.done.txt
    """
}
