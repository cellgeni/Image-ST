process VALIS_WARP {
    tag "${meta.id}"
    label 'process_high'

    // conda "${moduleDir}/environment.yml"
    container "quay.io/cellgeni/valis:1.2.0"

    // Singularity: VALIS writes to bf_formats.txt inside the Python env (read-only in Singularity).
    // A writable copy is created in the work dir and bind-mounted to the expected container path.
    containerOptions {
        if (workflow.containerEngine == 'singularity') {
            def bf = new File("${task.workDir}/bf_formats.txt")
            if (!bf.exists()) {
                bf.createNewFile()
            }
            "--bind ${bf.absolutePath}:/env/lib/python3.10/site-packages/valis/data/bf_formats.txt"
        }
        else {
            ''
        }
    }

    input:
    // results_dir is the full VALIS_REGISTER output dir; staging it makes
    // input_slides/ accessible at the same relative path the registrar expects.
    tuple val(meta), path(results_dir)

    output:
    tuple val(meta), path("${out_dir}"), emit: registered_slides
    tuple val("${task.process}"), val('valis'), eval('pip show valis | grep "^Version" | sed \'s/Version: //\''), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_dir = "${prefix}_valis_warp"
    """
    export VIPS_CONCURRENCY=1

    mkdir -p ${out_dir}
    valis-cli warp \\
        ${results_dir}/${prefix}/data/${prefix}_registrar.pickle \\
        ${out_dir} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_dir = "${prefix}_valis_warp"
    """
    mkdir -p ${out_dir}
    touch "${out_dir}/${prefix}_registered.ome.tiff"
    chmod -R a+w ${out_dir}
    """
}
