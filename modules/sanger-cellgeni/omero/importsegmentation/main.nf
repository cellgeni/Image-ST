process OMERO_IMPORTSEGMENTATION {
    tag "${meta.id}"
    label 'process_single'
    cpus 30

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/roi_convert_ngff:0.6.2'
        : 'quay.io/cellgeni/roi_convert_ngff:0.6.2'}"

    input:
    tuple val(meta), path(csv), val(image_id), val(host), val(table_name), val(roi_name), val(out_dir)

    secret 'OMERO_USER'
    secret 'OMERO_PASSWORD'

    output:
    tuple val(meta), path("*.importsegmentation.done.txt"), emit: done
    path "versions.yml", emit: versions

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

    python - <<'PY' > versions.yml
    import platform
    import subprocess

    version = "unknown"
    for cmd in (["ROI_Converter_NGFF", "--version"], ["ROI_Converter_NGFF", "-h"]):
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
            text = (proc.stdout or proc.stderr).strip()
            if text:
                version = text.splitlines()[0]
                break
        except Exception:
            pass

    print("${task.process}:")
    print(f"  python: {platform.python_version()}")
    print(f"  ROI_Converter_NGFF: {version}")
    PY
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    printf "%s\\n" "${roi_id}" > ${prefix}.importsegmentation.done.txt

    cat <<-EOF > versions.yml
    ${task.process}:
      python: stub
      ROI_Converter_NGFF: stub
    EOF
    """
}
