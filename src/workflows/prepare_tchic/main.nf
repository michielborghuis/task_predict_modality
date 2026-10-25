workflow auto {
  findStates(params, meta.config)
    | meta.workflow.run(auto: [publish: "state"])
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | prepare_tchic_input.run(
      fromState: { id, state ->
        [
          input_raw: state.input_raw,
          input_filtered: state.input_filtered,
          input_annotations: state.input_annotations,
          input_split: state.input_split,
          mark: state.mark,
          annotation_sample_prefix: state.annotation_sample_prefix,
          output_mod1: "tchic_raw_mod1.h5ad",
          output_mod2: "tchic_raw_mod2.h5ad"
        ]
      },
      toState: { id, output, state ->
        state + [raw_mod1: output.output_mod1, raw_mod2: output.output_mod2]
      }
    )
    | log_cp.run(
      key: "normalize_h3k4me3",
      fromState: { id, state ->
        [
          input: state.raw_mod1,
          output: "tchic_normalized_mod1.h5ad",
          n_cp: 10000,
          layer_output: "normalized",
          obs_size_factors: "size_factors",
          normalization_id: "log_cp10k"
        ]
      },
      toState: { id, output, state -> state + [normalized_mod1: output.output] }
    )
    | log_cp.run(
      key: "normalize_gex",
      fromState: { id, state ->
        [
          input: state.raw_mod2,
          output: "tchic_normalized_mod2.h5ad",
          n_cp: 10000,
          layer_output: "normalized",
          obs_size_factors: "size_factors",
          normalization_id: "log_cp10k"
        ]
      },
      toState: { id, output, state -> state + [normalized_mod2: output.output] }
    )
    | process_dataset.run(
      fromState: { id, state ->
        [
          input_mod1: state.normalized_mod1,
          input_mod2: state.normalized_mod2,
          output_train_mod1: "train_mod1.h5ad",
          output_train_mod2: "train_mod2.h5ad",
          output_test_mod1: "test_mod1.h5ad",
          output_test_mod2: "test_mod2.h5ad",
          dataset_id: "tchic_mgast_h3k4me3_genelevel_dev",
          swap: false,
          seed: state.seed,
          max_test_cells: state.max_test_cells,
          preserve_obs_columns: "batch,size_factors,sample,day,cell_type,celltype_general,lineage",
          preserve_var_columns: "gene_ids,hvg,hvg_score,feature_id,feature_name,Chromosome,Symbol,Biotype"
        ]
      },
      toState: [
        "output_train_mod1",
        "output_train_mod2",
        "output_test_mod1",
        "output_test_mod2"
      ]
    )
    | setState([
      "output_train_mod1",
      "output_train_mod2",
      "output_test_mod1",
      "output_test_mod2",
      "_meta"
    ])

  emit:
  output_ch
}
