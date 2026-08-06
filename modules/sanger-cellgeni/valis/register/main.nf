process VALIS_REGISTER {
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
    // stageAs puts reference and moving into separate sub-dirs, avoiding name collisions
    // when the same file is passed for both (e.g. self-registration checks).
    tuple val(meta), path(reference_img, stageAs: 'reference/*'), path(moving_imgs, stageAs: 'moving/*')

    output:
    tuple val(meta), path("${out_dir}"), emit: results_dir
    tuple val(meta), path("${reg_dir}/data/*_registrar.pickle"), emit: registrar
    tuple val(meta), path("${reg_dir}/data/*_transforms.json"), emit: transforms
    tuple val("${task.process}"), val('valis'), eval('pip show valis | grep "^Version" | sed \'s/Version: //\''), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_dir = "${prefix}_valis_register"
    // VALIS stores its data under {dst_dir}/{name}/ — pass --name so the path is predictable.
    reg_dir = "${out_dir}/${prefix}"
    """
    # pyvips reads TIFF strips in parallel by default; sequential TIFFs require
    # single-threaded access to avoid "out of order read" errors.
    export VIPS_CONCURRENCY=1

    valis-cli register \\
        ${reference_img} \\
        ${moving_imgs} \\
        ${out_dir} \\
        --name ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_dir = "${prefix}_valis_register"
    reg_dir = "${out_dir}/${prefix}"
    """
    mkdir -p ${reg_dir}/data
    touch "${reg_dir}/data/${prefix}_registrar.pickle"
    echo '{"@context":"https://ngff.openmicroscopy.org/rfc/5","registration_name":"${prefix}","coordinateTransformations":[]}' \
        > "${reg_dir}/data/${prefix}_transforms.json"
    chmod -R a+w ${out_dir}
    """
}
