process POSTCODE_POSTPROCESS {
    tag "${meta.id}"
    label 'process_single'

    container "quay.io/cellgeni/postcode:0.2.1"

    input:
    tuple val(meta), path(model_params_and_losses_path), path(starfish_codebook_path), path(spot_loc_path)

    output:
    tuple val(meta), path("${decoded_profile}"), emit: decoded_profile
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    decoded_profile = "${prefix}_decoded.csv"
    """
    decode_postprocess.py run \\
        --model_file_path ${model_params_and_losses_path} \\
        --spot_locations_path ${spot_loc_path} \\
        --starfish_codebook_path ${starfish_codebook_path} \\
        -decoded_spots_df_path ${decoded_profile} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        imaging: \$(decode_postprocess.py version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    decoded_profile = "${prefix}_decoded.csv"
    """
    echo ${args}
    
    touch ${decoded_profile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        imaging: \$(decode_postprocess.py version)
    END_VERSIONS
    """
}
