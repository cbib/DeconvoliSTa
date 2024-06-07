nextflow.enable.dsl=2

// Deconvolution methods
include { runRCTD } from './rctd/run_method.nf'
include { buildCell2locationModel; fitCell2locationModel} from './cell2location/run_method.nf'
// Helper functions
include { convertBetweenRDSandH5AD as convert_sc ; convertBetweenRDSandH5AD as convert_sp } from '../helper_processes'
include { formatTSVFile as formatStereoscope; formatTSVFile as formatC2L; formatTSVFile as formatDestVI; formatTSVFile as formatDSTG;
          formatTSVFile as formatTangram; formatTSVFile as formatSTRIDE } from '../helper_processes'
include { createDummyFile } from '../helper_processes'

workflow runMethods {
    take:
        sc_input_ch
        sp_input_ch
        sc_input_type
        sp_input_type

    main:

        output_ch = Channel.empty() // collect output channels


                buildCell2locationModel(sc_input_conv)

                // Repeat model output for each spatial file
                buildCell2locationModel.out.combine(sp_input_pair)
                .multiMap { model_sc_file, sp_file_h5ad, sp_file_rds ->
                            model: model_sc_file
                            sp_input: tuple sp_file_h5ad, sp_file_rds }
                .set{ c2l_combined_ch }

                fitCell2locationModel(c2l_combined_ch.sp_input,
                                      c2l_combined_ch.model)
                formatC2L(fitCell2locationModel.out)
                output_ch = output_ch.mix(formatC2L.out)
    emit:
        output_ch
}

workflow {

    if (!(params.mode ==~ /run_dataset/)){
        throw new Exception("Error: can only run this with the 'run_dataset' mode")
    }

    // RUN ON YOUR OWN DATA
    println("Running the pipeline on the provided data...")
    
    // Print inputs (the timing isn't right with with view(), so do this instead)
    // Although view() has the advantage that it gives the absolute path (sc_input_ch.view())
    if (params.verbose) {
        println("Single-cell reference dataset:")
        println(file(params.sc_input))

        println("\nSpatial dataset(s):") 
        // With glob pattern, there will be multiple files
        params.sp_input =~ /\*/ ? file(params.sp_input).each{println "$it"} : println (file(params.sp_input))
    }

    sc_input_ch = Channel.fromPath(params.sc_input) // Can only have 1 file
    sp_input_ch = Channel.fromPath(params.sp_input) // Can have one or more files

    runMethods(sc_input_ch, sp_input_ch)
}
