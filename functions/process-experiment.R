process_experiment <-
  function(experiment, prepped_fishery, cloud_dir) {
    ex <- readRDS(glue::glue("{cloud_dir}/experiment_{experiment}.rds"))

    print(experiment)

    sampled_years <- prepped_fishery$sampled_years

    pex <-
      process_scrooge(fit = ex,
                      sampled_years = sampled_years,
                      to_tidy = "f_t")


    performance <-
      judge_performance(observed = prepped_fishery$simed_fishery,
                        predicted = pex)

    out <- performance$comparison_summary

    return(out)

  }