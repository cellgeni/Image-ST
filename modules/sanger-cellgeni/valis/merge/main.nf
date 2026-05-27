process VALIS_MERGE {
    tag "${meta.id}"
    label 'process_medium'

    // conda "${moduleDir}/environment.yml"
    container "quay.io/cellgeni/valis:1.2.0"

    // Singularity: VALIS writes to bf_formats.txt inside the Python env (read-only in Singularity).
    // A writable copy is created in the work dir and bind-mounted to the expected container path.
    containerOptions {
        if (workflow.containerEngine == 'singularity') {
            def bf = new File("${task.workDir}/bf_formats.txt")
            if (!bf.exists()) bf.createNewFile()
            "--bind ${bf.absolutePath}:/env/lib/python3.10/site-packages/valis/data/bf_formats.txt"
        } else {
            ''
        }
    }

    input:
    // results_dir is the full VALIS_REGISTER output dir; staging it makes
    // input_slides/ accessible at the same relative path the registrar expects.
    tuple val(meta), path(results_dir), path(channel_names_json)

    output:
    tuple val(meta), path("${prefix}_merged.ome.tiff"), emit: merged
    tuple val("${task.process}"), val('valis'), eval('pip show valis | grep "^Version" | sed \'s/Version: //\''), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args               = task.ext.args ?: ''
    prefix                 = task.ext.prefix ?: "${meta.id}"
    def channel_names_arg  = (channel_names_json.name && channel_names_json.name != 'NO_FILE') ? "--channel-names-json ${channel_names_json}" : ''
    """
    valis-cli merge \\
        ${results_dir}/${prefix}/data/${prefix}_registrar.pickle \\
        ${prefix}_merged.ome.tiff \\
        ${channel_names_arg} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_merged.ome.tiff"
    chmod a+w "${prefix}_merged.ome.tiff"
    """
}
