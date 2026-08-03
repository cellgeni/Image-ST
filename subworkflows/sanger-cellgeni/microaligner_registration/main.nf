#!/usr/bin/env/ nextflow

include { IMAGING_MICROALIGNER                            } from '../../../modules/sanger-cellgeni/imaging/microaligner/main'

include { IMAGING_MICROALIGNER as MICROALIGNER_FEATREG ; IMAGING_MICROALIGNER as MICROALIGNER_OPTFLOWREG } from '../../../modules/sanger-cellgeni/imaging/microaligner/main'


workflow MICRO_ALIGNER_REGISTRATION {
    take:
    images

    main:

    ch_versions = channel.empty()
    GENERATE_FEAT_REG_YAML(images)
    GENERATE_OPTFLOW_REG_YAML(images)
    MICROALIGNER_FEATREG(GENERATE_FEAT_REG_YAML.out.combine(images, by: 0), "feature")
    ch_versions = ch_versions.mix(MICROALIGNER_FEATREG.out.versions.first())
    MICROALIGNER_OPTFLOWREG(GENERATE_OPTFLOW_REG_YAML.out.combine(MICROALIGNER_FEATREG.out.registered_image, by: 0), "optflow")
    ch_versions = ch_versions.mix(MICROALIGNER_OPTFLOWREG.out.versions.first())

    emit:
    image    = MICROALIGNER_OPTFLOWREG.out.registered_image // channel: [ val(meta), [ image ] ]
    versions = ch_versions // channel: [ versions.yml ]
}


process GENERATE_FEAT_REG_YAML {
    tag "${meta.id}"

    publishDir params.out_dir + "/registration_configs"

    input:
    tuple val(meta), path(images)

    output:
    tuple val(meta), path("${meta.id}_feat_reg.yaml")

    script:
    def cycles_str = ""
    images.eachWithIndex { img, i ->
        cycles_str += "        Cycle ${i + 1}: ./${img}\n    "
    }
    """
    echo "# Input
    # If your input is not structured as a multichannel stack, you can specify
    # individual channels per cycle like described bellow.
    # If file contains more than one page, they will be treated as z-planes.

    Input:
        InputImagePaths:
    ${cycles_str}
        ReferenceCycle: ${params.reference_cycle}
        ReferenceChannel: ${params.reference_channel}

    # Output
    # Images will be saved to a directory
    # And because SaveOutputToCycleStack is true in this example,
    # and OutputPrefix is empty,
    # the file names will have the following patterns:
    #   feature_reg_result_stack.tif
    #   optflow_reg_result_stack.tif

    Output:
        OutputDir: "./"
        OutputPrefix: ${meta.id}_
        SaveOutputToCycleStack: true

    # Registration parameters
    # Will do only linear feature based registration,
    # will not do non-linear optical flow based registration.
    # For information on all the registration parameters check the file
    # config_with_all_info.yaml

    RegistrationParameters:
        FeatureReg:
            NumberPyramidLevels: 3
            NumberIterationsPerLevel: 3
            TileSize: 2000
            Overlap: 100
            NumberOfWorkers: 0
            UseFullResImage: false
            UseDOG: true" >> ./${meta.id}_feat_reg.yaml
    """
}


process GENERATE_OPTFLOW_REG_YAML {
    tag "${meta.id}"

    publishDir params.out_dir + "/registration_configs"

    input:
    tuple val(meta), path(images)

    output:
    tuple val(meta), path("${meta.id}_optflow_reg.yaml")

    script:
    """
    echo "# Input
    # If your input image is a stack that contains channels from all cycles

    Input:
        InputImagePaths:
            CycleStack: ${meta.id}_feature_reg_result_stack.tif
        ReferenceCycle: ${params.reference_cycle}
        ReferenceChannel: ${params.reference_channel}

    # Output
    # Images will be saved to a directory
    #   "/path/to/out/registered_imgs"
    # And because SaveOutputToCycleStack is false in this example,
    # the file names will have the following patterns:
    #   experiment_002_feature_reg_result_cyc001.tif
    #   experiment_002_feature_reg_result_cyc00N.tif
    #   experiment_002_optflow_reg_result_cyc001.tif
    #   experiment_002_optflow_reg_result_cyc00N.tif

    Output:
        OutputDir: "./"
        OutputPrefix: "${meta.id}"_
        SaveOutputToCycleStack: true

    # Registration parameters
    # Will do only non-linear optical flow based registration,
    # will skip the feature based registration.
    # For information on all the registration parameters check the file
    # config_with_all_info.yaml

    RegistrationParameters:
        OptFlowReg:
            NumberPyramidLevels: 3
            NumberIterationsPerLevel: 3
            TileSize: 1000
            Overlap: 100
            NumberOfWorkers: 0
            UseFullResImage: true
            UseDOG: false" >> ./${meta.id}_optflow_reg.yaml
    """
}
