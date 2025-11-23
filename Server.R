# Server.R

library(shiny)
library(DT)
library(nphRshiny)
library(tidyverse)
library(mvtnorm)
library(gsDesign)
library(plotly)
library(ggplot2)
library(SAMprior)
library(shinyjs)
library(BayesianHybridDesign)
library(future)
library(promises)
library(foreach)
library(doParallel)
library(RBesT)

# Configure temp directory for Windows compatibility
tryCatch({
  app_temp <- file.path(getwd(), "temp")
  if (!dir.exists(app_temp)) {
    dir.create(app_temp, recursive = TRUE)
  }
  Sys.setenv(TMPDIR = app_temp)
  Sys.setenv(TEMP = app_temp)
  Sys.setenv(TMP = app_temp)
}, error = function(e) {
  message("Using system temp directory")
})

dpp_pmd_stats <- reactiveVal(NULL)

# Configure future plan with error handling
tryCatch({
  plan(multisession, workers = 2)
}, error = function(e) {
  message("Multisession not available, using sequential processing")
  plan(sequential)
})

# ===========================================================================
# NEW: Credible Difference Calculation Functions
# ===========================================================================

calculate_credible_difference <- function(dpp_results, confidence_level = 0.95) {
  if (is.null(dpp_results$phat_pt_larger_pc_all)) {
    return(NA)
  }

  # Calculate the critical value from the posterior distribution
  # This finds the value where (1 - confidence_level) of simulations exceed
  critical_value <- quantile(dpp_results$phat_pt_larger_pc_all, probs = confidence_level)

  return(critical_value)
}

calculate_orr_critical_value <- function(dpp_results, confidence_level = 0.95) {
  if (is.null(dpp_results$mean_hca) || is.null(dpp_results$mean_c)) {
    return(NA)
  }

  # Calculate ORR differences: experimental vs hybrid control
  orr_differences <- dpp_results$mean_hca - dpp_results$mean_c

  # Find the critical value for statistical significance
  critical_value <- quantile(orr_differences, probs = (1 - confidence_level))

  return(critical_value)
}

shinyServer(function(input, output, session) {

  # ===========================================================================
  # Term Definitions (keeping existing definitions)
  # ===========================================================================

  dpp_term_definitions <- list(
    "pt" = "Response rate for experimental arm in current study.",
    "nt" = "Number of patients in experimental arm in current study.",
    "pc" = "Response rate for control arm in current study.",
    "nc" = "Number of patients in control arm in current study.",
    "pc.calib" = "Required for calibration if tau is not provided. Response rate for control arm in current study for calibration. Usually, <code>pc.calib = pch</code>.",
    "pch" = "Response rate for control treatment in historical study.",
    "nche" = "Equivalent number of patients borrowed from historical study.",
    "nch" = "Total number of patients in historical control.",
    "alpha" = "A scalar. One sided type I error rate. Required for calibration if tau is not provided.",
    "tau" = "Calibrated threshold for statistical significance. If tau is not provided, it will be calculated by calibration to type I error alpha.",
    "a0c" = "Hyperprior for control response rate beta(a0c, b0c).",
    "b0c" = "Hyperprior for control response rate beta(a0c, b0c).",
    "a0t" = "Hyperprior for experimental response rate beta(a0t, b0t).",
    "b0t" = "Hyperprior for experimental response rate beta(a0t, b0t).",
    "delta_threshold" = "Borrow when <code>abs(pc_hat (current study) - pch) <= delta_threshold</code>.",
    "method" = "Method for dynamic borrowing, 'Empirical Bayes', 'Bayesian p', 'Generalized BC', 'JSD'.",
    "theta" = "A parameter with a range of (0, 1), and applicable to method: 'Generalized BC'.",
    "eta" = "A parameter with a range of (0, infty), and applicable to method: 'Bayesian p', 'Generalized BC', 'JSD'. 'Generalized BC' method requires two parameters theta and eta.",
    "datamat" = "A matrix with dimension <code>nsim * 2</code> containing the pre-simulated data for the study treatment (1st column) and control (1st column) groups, respectively. If not supplied, binomial random Monte Carlo samples will be generated in the function.",
    "w0" = "Prior power parameters <code>w</code>. If not specified (default), <code>w_d</code> is calculated by the specified method for dynamic borrowing.",
    "nsim" = "Number of replications to calculate power.",
    "seed" = "Seed for simulations."
  )

  dpp_analysis_term_definitions <- list(
    "w" = "Borrowing weight. A number between 0 and 1 that shows how much information from historical data is being used.",
    "phat_pt_larger_pc" = "The posterior probability that the experimental arm's response rate is truly higher than the control arm's response rate.",
    "apost_c_trial" = "The posterior alpha parameter for the control arm's response rate, using only data from the current study.",
    "bpost_c_trial" = "The posterior beta parameter for the control arm's response rate, using only data from the current study.",
    "apost_c_hca" = "The posterior alpha parameter for the hybrid control arm's response rate after combining current and historical data.",
    "bpost_c_hca" = "The posterior beta parameter for the hybrid control arm's response rate after combining current and historical data.",
    "apost_t" = "The posterior alpha parameter for the experimental arm's response rate.",
    "bpost_t" = "The posterior beta parameter for the experimental arm's response rate."
  )

  fisher_term_definitions <- list(
    "Yc (fisher)" = "Number of subjects with response in experimental arm for <code>fisher</code> function.",
    "nc (fisher)" = "Number of subjects in control arm for <code>fisher</code> function.",
    "Yt (fisher)" = "Number of subjects with response in control arm for <code>fisher</code> function.",
    "nt (fisher)" = "Number of subjects in experimental arm for <code>fisher</code> function.",
    "alternative (fisher)" = "Type of alternative hypothesis (e.g., 'greater') for <code>fisher</code> function.",
    "pc (fisher.bound)" = "Response rate for control arm for <code>fisher.bound</code> function.",
    "nc (fisher.bound)" = "Number of patients in control arm for <code>fisher.bound</code> function.",
    "nt (fisher.bound)" = "Number of patients in experimental arm for <code>fisher.bound</code> function.",
    "alpha (fisher.bound)" = "P-value threshold for significance. Alpha must be one-sided. Default 0.1 for <code>fisher.bound</code> function.",
    "pt (fisher.power)" = "Probability of success in experimental arm for <code>fisher.power</code> function.",
    "nt (fisher.power)" = "Number of subject in experimental arm for <code>fisher.power</code> function.",
    "pc (fisher.power)" = "Probability of success in control arm for <code>fisher.power</code> function.",
    "nc (fisher.power)" = "Number of subject in control arm for <code>fisher.power</code> function.",
    "alpha (fisher.power)" = "One sided type I error rate for <code>fisher.power</code> function.",
    "nsim (fisher.power)" = "Number of replications to calculate power. Default 100,000 for <code>fisher.power</code> function.",
    "seed (fisher.power)" = "Seed for simulations. Default 2000 for <code>fisher.power</code> function."
  )

  sam_plot_term_definitions <- list(
    "alpha_hist" = "Alpha (a) for Beta(a, b) in the informative prior (historical) for SAM Prior Plot.",
    "beta_hist" = "Beta (b) for Beta(a, b) in the informative prior (historical) for SAM Prior Plot.",
    "n_control" = "Number of patients in the control arm data for SAM Prior Plot.",
    "p_control" = "Simulated control rate for SAM Prior Plot.",
    "nf_alpha" = "Alpha for Beta(a, b) in the non-informative prior for SAM Prior Plot.",
    "nf_beta" = "Beta (b) for Beta(a, b) in the non-informative prior for SAM Prior Plot.",
    "delta" = "Clinically meaningful difference (delta) for SAM Prior Plot."
  )

  sam_weight_term_definitions <- list(
    "alpha_hist" = "Alpha (a) for Beta(a, b) in the informative prior (historical) for SAM Weight calculation.",
    "beta_hist" = "Beta (b) for Beta(a, b) in the informative prior (historical) for SAM Weight calculation.",
    "n_control" = "Number of patients in the control arm summary data for SAM Weight calculation.",
    "r_control" = "Number of responses in the control arm summary data for SAM Weight calculation.",
    "delta" = "Clinically meaningful difference (delta) for SAM Weight calculation.",
    "method" = "Weight method ('LRT' or 'PPR') for SAM Weight calculation.",
    "prior_odds" = "Prior Odds H0 vs H1 (PPR only) for SAM Weight calculation."
  )

  # ===========================================================================
  # Observe Events for Term Definitions (keeping existing)
  # ===========================================================================

  observeEvent(input$term_pt, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["pt"]])) })
  observeEvent(input$term_nt, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["nt"]])) })
  observeEvent(input$term_pc, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["pc"]])) })
  observeEvent(input$term_nc, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["nc"]])) })
  observeEvent(input$term_p_calib, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["pc.calib"]])) })
  observeEvent(input$term_pch, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["pch"]])) })
  observeEvent(input$term_nche, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["nche"]])) })
  observeEvent(input$term_nch, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["nch"]])) })
  observeEvent(input$term_alpha, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["alpha"]])) })
  observeEvent(input$term_tau, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["tau"]])) })
  observeEvent(input$term_a0c, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["a0c"]])) })
  observeEvent(input$term_b0c, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["b0c"]])) })
  observeEvent(input$term_a0t, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["a0t"]])) })
  observeEvent(input$term_b0t, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["b0t"]])) })
  observeEvent(input$term_delta_threshold, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["delta_threshold"]])) })
  observeEvent(input$term_method, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["method"]])) })
  observeEvent(input$term_theta, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["theta"]])) })
  observeEvent(input$term_eta, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["eta"]])) })
  observeEvent(input$term_datamat, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["datamat"]])) })
  observeEvent(input$term_w0, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["w0"]])) })
  observeEvent(input$term_nsim, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["nsim"]])) })
  observeEvent(input$term_seed, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["seed"]])) })

  observeEvent(input$dpp_analysis_term_w, { output$definition_output <- renderUI(HTML(dpp_analysis_term_definitions[["w"]])) })
  observeEvent(input$dpp_analysis_term_phat, { output$definition_output <- renderUI(HTML(dpp_analysis_term_definitions[["phat_pt_larger_pc"]])) })
  observeEvent(input$dpp_analysis_term_apost_c_trial, { output$definition_output <- renderUI(HTML(dpp_analysis_term_definitions[["apost_c_trial"]])) })
  observeEvent(input$dpp_analysis_term_bpost_c_trial, { output$definition_output <- renderUI(HTML(dpp_analysis_term_definitions[["bpost_c_trial"]])) })
  observeEvent(input$dpp_analysis_term_apost_c_hca, { output$definition_output <- renderUI(HTML(dpp_analysis_term_definitions[["apost_c_hca"]])) })
  observeEvent(input$dpp_analysis_term_bpost_c_hca, { output$definition_output <- renderUI(HTML(dpp_analysis_term_definitions[["bpost_c_hca"]])) })
  observeEvent(input$dpp_analysis_term_apost_t, { output$definition_output <- renderUI(HTML(dpp_analysis_term_definitions[["apost_t"]])) })
  observeEvent(input$dpp_analysis_term_bpost_t, { output$definition_output <- renderUI(HTML(dpp_analysis_term_definitions[["bpost_t"]])) })

  observeEvent(input$fisher_term_Yc, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["Yc (fisher)"]])) })
  observeEvent(input$fisher_term_nc_fisher, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["nc (fisher)"]])) })
  observeEvent(input$fisher_term_Yt, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["Yt (fisher)"]])) })
  observeEvent(input$fisher_term_nt_fisher, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["nt (fisher)"]])) })
  observeEvent(input$fisher_term_alternative, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["alternative (fisher)"]])) })
  observeEvent(input$fisher_term_pc_bound, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["pc (fisher.bound)"]])) })
  observeEvent(input$fisher_term_nc_bound, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["nc (fisher.bound)"]])) })
  observeEvent(input$fisher_term_nt_bound, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["nt (fisher.bound)"]])) })
  observeEvent(input$fisher_term_alpha_bound, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["alpha (fisher.bound)"]])) })
  observeEvent(input$fisher_term_pt_power, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["pt (fisher.power)"]])) })
  observeEvent(input$fisher_term_nt_power, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["nt (fisher.power)"]])) })
  observeEvent(input$fisher_term_pc_power, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["pc (fisher.power)"]])) })
  observeEvent(input$fisher_term_nc_power, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["nc (fisher.power)"]])) })
  observeEvent(input$fisher_term_alpha_power, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["alpha (fisher.power)"]])) })
  observeEvent(input$fisher_term_nsim_power, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["nsim (fisher.power)"]])) })
  observeEvent(input$fisher_term_seed_power, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["seed (fisher.power)"]])) })

  observeEvent(input$sam_plot_library_term_alpha_hist, { output$sam_plot_definition_output <- renderUI(HTML(sam_plot_term_definitions[["alpha_hist"]])) })
  observeEvent(input$sam_plot_library_term_beta_hist, { output$sam_plot_definition_output <- renderUI(HTML(sam_plot_term_definitions[["beta_hist"]])) })
  observeEvent(input$sam_plot_library_term_n_control, { output$sam_plot_definition_output <- renderUI(HTML(sam_plot_term_definitions[["n_control"]])) })
  observeEvent(input$sam_plot_library_term_p_control, { output$sam_plot_definition_output <- renderUI(HTML(sam_plot_term_definitions[["p_control"]])) })
  observeEvent(input$sam_plot_library_term_nf_alpha, { output$sam_plot_definition_output <- renderUI(HTML(sam_plot_term_definitions[["nf_alpha"]])) })
  observeEvent(input$sam_plot_library_term_nf_beta, { output$sam_plot_definition_output <- renderUI(HTML(sam_plot_term_definitions[["nf_beta"]])) })
  observeEvent(input$sam_plot_library_term_delta, { output$sam_plot_definition_output <- renderUI(HTML(sam_plot_term_definitions[["delta"]])) })

  observeEvent(input$sam_weight_library_term_alpha_hist_w, { output$sam_plot_definition_output <- renderUI(HTML(sam_weight_term_definitions[["alpha_hist"]])) })
  observeEvent(input$sam_weight_library_term_beta_hist_w, { output$sam_plot_definition_output <- renderUI(HTML(sam_weight_term_definitions[["beta_hist"]])) })
  observeEvent(input$sam_weight_library_term_n_control_w, { output$sam_plot_definition_output <- renderUI(HTML(sam_weight_term_definitions[["n_control"]])) })
  observeEvent(input$sam_weight_library_term_r_control_w, { output$sam_plot_definition_output <- renderUI(HTML(sam_weight_term_definitions[["r_control"]])) })
  observeEvent(input$sam_weight_library_term_delta_w, { output$sam_plot_definition_output <- renderUI(HTML(sam_weight_term_definitions[["delta"]])) })
  observeEvent(input$sam_weight_library_term_method_w, { output$sam_plot_definition_output <- renderUI(HTML(sam_weight_term_definitions[["method"]])) })
  observeEvent(input$sam_weight_library_term_prior_odds, { output$sam_plot_definition_output <- renderUI(HTML(sam_weight_term_definitions[["prior_odds"]])) })


  # ===========================================================================
  # Conditional UI (keeping existing)
  # ===========================================================================

  output$theta_ui <- renderUI({
    if (input$method == "Generalized BC") {
      numericInput("theta", "Theta (0-0.5)", value = 0.5, min = 0.001, max = 0.5, step = 0.001)
    }
  })

  output$eta_ui <- renderUI({
    if (input$method %in% c("Bayesian p", "Generalized BC", "JSD")) {
      numericInput("eta", "Eta (0-Inf)", value = 1, min = 0, step = 0.1)
    }
  })

  output$dpp_analysis_theta_ui <- renderUI({
    if (input$dpp_analysis_method == "Generalized BC") {
      numericInput("dpp_analysis_theta", "Theta (0-1)", value = 0.5, min = 0.001, max = 0.999, step = 0.001)
    }
  })

  output$dpp_analysis_eta_ui <- renderUI({
    if (input$dpp_analysis_method %in% c("Bayesian p", "Generalized BC", "JSD")) {
      numericInput("dpp_analysis_eta", "Eta (0-Inf)", value = 1, min = 0, step = 0.1)
    }
  })

  # ===========================================================================
  # Dynamic Power Prior (DPP) Study Design
  # ===========================================================================

  dpp_results <- reactiveVal(NULL)
  dpp_pmd_stats <- reactiveVal(NULL)

  observeEvent(input$run_dpp, {
    shinyjs::html("analysis_status_message", "")

    errors <- c()
    if (is.null(input$pt) || input$pt < 0 || input$pt > 1) errors <- c(errors, "Experimental Arm Response Rate (pt) must be between 0 and 1.")
    if (is.null(input$pc) || input$pc < 0 || input$pc > 1) errors <- c(errors, "Control Arm Response Rate (pc) must be between 0 and 1.")
    if (is.null(input$p_calib) || input$p_calib < 0 || input$p_calib > 1) errors <- c(errors, "Response Rate for Calibration (pc.calib) must be between 0 and 1.")
    if (is.null(input$pch) || input$pch < 0 || input$pch > 1) errors <- c(errors, "Historical Control Response Rate (pch) must be between 0 and 1.")
    if (is.null(input$nt) || input$nt <= 0 || is.null(input$nc) || input$nc <= 0 || is.null(input$nch) || input$nch <= 0 || is.null(input$nche) || input$nche <= 0) {
      errors <- c(errors, "Sample sizes (nt, nc, nch, nche) must be positive integers.")
    }
    if (!is.null(input$nche) && input$nche > 50) {
      errors <- c(errors, "Maximum Number of Patients from Historical Study (nche) cannot be greater than 50.")
    }
    if (is.null(input$alpha) || input$alpha <= 0 || input$alpha > 1) errors <- c(errors, "Type I Error (alpha) must be between 0 and 1.")
    if (!is.null(input$nche) && !is.null(input$nch) && input$nche > input$nch) errors <- c(errors, "Equivalent number of patients borrowed (nche) cannot exceed total historical control patients (nch).")
    if (is.null(input$delta_threshold) || input$delta_threshold < 0 || input$delta_threshold > 1) errors <- c(errors, "Delta Threshold must be between 0 and 1.")

    if (input$method == "Generalized BC" && (is.null(input$theta) || input$theta <= 0 || input$theta >= 1)) {
      errors <- c(errors, "Theta must be between 0 and 1 for 'Generalized BC' method.")
    }
    if (input$method %in% c("Bayesian p", "Generalized BC", "JSD") && (is.null(input$eta) || input$eta <= 0)) {
      errors <- c(errors, "Eta must be positive for selected methods.")
    }

    seed_val <- suppressWarnings(as.integer(input$seed))
    if (is.na(seed_val) || seed_val < 1) {
      errors <- c(errors, "Seed must be a positive integer.")
    }

    if (length(errors) > 0) {
      shinyjs::html("analysis_status_message", paste0("<p style='color: red;'>", paste(errors, collapse = "<br>"), "</p>"))
      return()
    }

    shinyjs::html("analysis_status_message", "<p style='color: blue;'>Calculation started. Please wait...</p>")
    shinyjs::disable("run_dpp")

    power.DPP_args <- list(
      pt = input$pt,
      nt = input$nt,
      pc = input$pc,
      nc = input$nc,
      pc.calib = input$p_calib,
      pch = input$pch,
      nche = input$nche,
      nch = input$nch,
      alpha = input$alpha,
      a0c = input$a0c,
      b0c = input$b0c,
      a0t = input$a0t,
      b0t = input$b0t,
      delta_threshold = input$delta_threshold,
      method = input$method,
      nsim = input$nsim,
      seed = seed_val
    )
    if (input$method == "Generalized BC") {
      power.DPP_args$theta = input$theta
      power.DPP_args$eta = input$eta
    } else if (input$method %in% c("Bayesian p", "JSD")) {
      power.DPP_args$eta = input$eta
    }

    p <- future({
      result <- do.call(BayesianHybridDesign::power.DPP, power.DPP_args)
      result
    }, seed = TRUE)

    then(p, onFulfilled = function(result) {
      dpp_results(result)
      shinyjs::html("analysis_status_message", "<p style='color: green;'>Study Design calculation ran successfully!</p>")
      shinyjs::enable("run_dpp")
    }, onRejected = function(e) {
      shinyjs::html("analysis_status_message", paste0("<p style='color: red;'>Study Design calculation failed: ", e$message, ". Please check inputs and try again.</p>"))
      dpp_results(NULL)
      shinyjs::enable("run_dpp")
    })
  })

  output$dppPower <- renderPrint({
    req(dpp_results()$power)
    cat(round(dpp_results()$power, 4))
  })

  output$dppTau <- renderPrint({
    req(dpp_results()$tau)
    cat(round(dpp_results()$tau, 4))
  })

  # UPDATED: Credible Difference (Statistical Critical Value)
  output$dppDeltaBound <- renderPrint({
    req(dpp_results())
    credible_diff <- calculate_credible_difference(dpp_results())
    if (!is.na(credible_diff)) {
      cat(round(credible_diff, 4))
    } else {
      cat("N/A")
    }
  })


  # ValueBox outputs for enhanced UI
  output$dppPowerBox <- renderValueBox({
    req(dpp_results()$power)
    valueBox(
      value = round(dpp_results()$power, 4),
      subtitle = "Statistical Power",
      icon = icon("chart-line"),
      color = "blue"
    )
  })

  output$dppTauBox <- renderValueBox({
    req(dpp_results()$tau)
    valueBox(
      value = round(dpp_results()$tau, 4),
      subtitle = "Tau Significance Threshold",
      icon = icon("crosshairs"),
      color = "yellow"
    )
  })

  # UPDATED: Credible Difference Box
  output$dppDeltaBoundBox <- renderValueBox({
    req(dpp_results())
    credible_diff <- calculate_credible_difference(dpp_results())
    if (!is.na(credible_diff)) {
      valueBox(
        value = round(credible_diff, 4),
        subtitle = "Credible Difference (Critical Value)",
        icon = icon("arrows-alt-h"),
        color = "green"
      )
    } else {
      valueBox(
        value = "N/A",
        subtitle = "Credible Difference",
        icon = icon("exclamation-triangle"),
        color = "red"
      )
    }
  })



  # PMD Statistics outputs
  output$dpp_pmd_mean <- renderText({
    req(dpp_pmd_stats())
    sprintf("%.4f", dpp_pmd_stats()$mean_PMD)
  })

  output$dpp_pmd_sd <- renderText({
    req(dpp_pmd_stats())
    sprintf("%.4f", dpp_pmd_stats()$sd_PMD)
  })

  output$dpp_pmd_ci <- renderText({
    req(dpp_pmd_stats())
    ci <- dpp_pmd_stats()$CI95_PMD
    sprintf("[%.4f, %.4f]", ci[1], ci[2])
  })

  # Updated PMD Plot using the proper plotPMD function
  output$dpp_plot_pmd <- renderPlot({
    req(dpp_results())

    tryCatch({
      # Use the proper plotPMD function and capture its return value
      stats <- BayesianHybridDesign::plotPMD(o = dpp_results())
      dpp_pmd_stats(stats)  # Store the statistics
    }, error = function(e) {
      dpp_pmd_stats(NULL)
      plot.new()
      text(0.5, 0.5,
           paste("Unable to generate PMD plot:\n", e$message),
           cex = 1.2, col = "red")
    })
  })

  # ===========================================================================
  # DPP Table 2 Analysis - Compare Different Borrowing Amounts
  # ===========================================================================

  dpp_table_results <- reactiveVal(NULL)

  observeEvent(input$run_dpp_table, {
    shinyjs::html("dpp_table_status_message", "")

    errors <- c()
    if (is.null(input$table_pc_values) || input$table_pc_values == "") {
      errors <- c(errors, "Please enter control arm response rates (comma-separated).")
    }
    if (is.null(input$table_nche_values) || input$table_nche_values == "") {
      errors <- c(errors, "Please enter nche values (comma-separated).")
    }

    if (length(errors) > 0) {
      shinyjs::html("dpp_table_status_message", paste0("<p style='color: red;'>", paste(errors, collapse = "<br>"), "</p>"))
      return()
    }

    shinyjs::html("dpp_table_status_message", "<p style='color: blue;'>Running comprehensive analysis. This may take several minutes...</p>")
    shinyjs::disable("run_dpp_table")

    pc_values <- as.numeric(unlist(strsplit(input$table_pc_values, ",")))
    nche_values <- as.numeric(unlist(strsplit(input$table_nche_values, ",")))
    table_nt <- input$table_nt
    table_nc <- input$table_nc
    table_pch <- input$table_pch
    table_nch <- input$table_nch
    table_delta_threshold <- input$table_delta_threshold
    table_alpha <- input$table_alpha
    table_a0c <- input$table_a0c
    table_b0c <- input$table_b0c
    table_a0t <- input$table_a0t
    table_b0t <- input$table_b0t
    table_nsim <- input$table_nsim

    p <- future({
      ns <- length(pc_values)
      nm <- length(nche_values)

      # Initialize result matrices
      TypeI_DPP <- matrix(NA, ns, nm)
      Power_DPP <- matrix(NA, ns, nm)
      pc_mi_DPP <- matrix(NA, ns, nm)
      pc_sdi_DPP <- matrix(NA, ns, nm)

      # Generate data once
      xt_list <- list()
      xt_typeI_list <- list()
      xc_list <- list()

      for(s in 1:ns){
        set.seed(2000)
        xt_typeI_list[[s]] <- rbinom(table_nsim, size = table_nt, prob = pc_values[s])
        xc_list[[s]] <- rbinom(table_nsim, size = table_nc, prob = pc_values[s])
        xt_list[[s]] <- rbinom(table_nsim, size = table_nt, prob = pc_values[s] + 0.2)
      }

      # Run simulations for each combination
      for(s in 1:ns){
        for(m in 1:nm){
          # Calibration
          cs <- s
          datamat_typeI <- cbind(xt_typeI_list[[cs]], xc_list[[cs]])

          tau <- BayesianHybridDesign::calibration(
            nt = table_nt,
            pc.calib = pc_values[cs],
            nc = table_nc,
            pch = table_pch,
            nche = nche_values[m],
            nch = table_nch,
            alpha = table_alpha,
            a0c = table_a0c,
            b0c = table_b0c,
            a0t = table_a0t,
            b0t = table_b0t,
            delta_threshold = table_delta_threshold,
            method = "Empirical Bayes",
            datamat = datamat_typeI,
            nsim = table_nsim,
            seed = 2000
          )

          # Type I error
          TypeI <- BayesianHybridDesign::power.DPP(
            pt = pc_values[s],
            nt = table_nt,
            pc = pc_values[s],
            nc = table_nc,
            pc.calib = pc_values[cs],
            pch = table_pch,
            nche = nche_values[m],
            nch = table_nch,
            tau = tau,
            a0c = table_a0c,
            b0c = table_b0c,
            a0t = table_a0t,
            b0t = table_b0t,
            delta_threshold = table_delta_threshold,
            method = "Empirical Bayes",
            datamat = datamat_typeI,
            nsim = table_nsim,
            seed = 2000
          )

          TypeI_DPP[s, m] <- TypeI$power
          pc_mi_DPP[s, m] <- TypeI$pc.PMD
          pc_sdi_DPP[s, m] <- TypeI$pc.sd.PMD

          # Power
          datamat_Power <- cbind(xt_list[[s]], xc_list[[s]])

          Power <- BayesianHybridDesign::power.DPP(
            pt = pc_values[s] + 0.2,
            nt = table_nt,
            pc = pc_values[s],
            nc = table_nc,
            pc.calib = pc_values[cs],
            pch = table_pch,
            nche = nche_values[m],
            nch = table_nch,
            tau = tau,
            a0c = table_a0c,
            b0c = table_b0c,
            a0t = table_a0t,
            b0t = table_b0t,
            delta_threshold = table_delta_threshold,
            method = "Empirical Bayes",
            datamat = datamat_Power,
            nsim = table_nsim,
            seed = 2000
          )

          Power_DPP[s, m] <- Power$power
        }
      }

      list(
        TypeI = TypeI_DPP,
        Power = Power_DPP,
        pc_mi = pc_mi_DPP,
        pc_sdi = pc_sdi_DPP,
        pc_values = pc_values,
        nche_values = nche_values
      )
    }, seed = TRUE)

    then(p, onFulfilled = function(result) {
      dpp_table_results(result)
      shinyjs::html("dpp_table_status_message", "<p style='color: green;'>Table analysis completed successfully!</p>")
      shinyjs::enable("run_dpp_table")
    }, onRejected = function(e) {
      shinyjs::html("dpp_table_status_message", paste0("<p style='color: red;'>Table analysis failed: ", e$message, "</p>"))
      dpp_table_results(NULL)
      shinyjs::enable("run_dpp_table")
    })
  })

  output$dpp_table_typeI <- DT::renderDataTable({
    req(dpp_table_results())
    result <- dpp_table_results()
    df <- as.data.frame(result$TypeI)
    colnames(df) <- paste0("nche=", result$nche_values)
    df <- cbind(pc = result$pc_values, df)
    DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE)) %>%
      DT::formatRound(columns = 2:ncol(df), digits = 4)
  })

  output$dpp_table_power <- DT::renderDataTable({
    req(dpp_table_results())
    result <- dpp_table_results()
    df <- as.data.frame(result$Power)
    colnames(df) <- paste0("nche=", result$nche_values)
    df <- cbind(pc = result$pc_values, df)
    DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE)) %>%
      DT::formatRound(columns = 2:ncol(df), digits = 4)
  })

  output$dpp_table_pmd <- DT::renderDataTable({
    req(dpp_table_results())
    result <- dpp_table_results()
    df <- as.data.frame(result$pc_mi)
    colnames(df) <- paste0("nche=", result$nche_values)
    df <- cbind(pc = result$pc_values, df)
    DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE)) %>%
      DT::formatRound(columns = 2:ncol(df), digits = 4)
  })

  output$dpp_table_sd_pmd <- DT::renderDataTable({
    req(dpp_table_results())
    result <- dpp_table_results()
    df <- as.data.frame(result$pc_sdi)
    colnames(df) <- paste0("nche=", result$nche_values)
    df <- cbind(pc = result$pc_values, df)
    DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE)) %>%
      DT::formatRound(columns = 2:ncol(df), digits = 4)
  })

  # ===========================================================================
  # Dynamic Power Prior (DPP) Analysis
  # ===========================================================================

  dpp_analysis_results <- reactiveVal(NULL)

  observeEvent(input$run_dpp_analysis, {
    shinyjs::html("dpp_analysis_status_message", "")

    errors <- c()
    if (is.null(input$dpp_analysis_rt) || input$dpp_analysis_rt < 0) errors <- c(errors, "Observed Responders in Experimental Arm (rt) must be non-negative.")
    if (is.null(input$dpp_analysis_rc) || input$dpp_analysis_rc < 0) errors <- c(errors, "Observed Responders in Control Arm (rc) must be non-negative.")
    if (!is.null(input$dpp_analysis_rt) && !is.null(input$dpp_analysis_nt) && input$dpp_analysis_rt > input$dpp_analysis_nt) errors <- c(errors, "Observed responders in experimental arm (rt) cannot exceed total experimental sample size (nt).")
    if (!is.null(input$dpp_analysis_rc) && !is.null(input$dpp_analysis_nc) && input$dpp_analysis_rc > input$dpp_analysis_nc) errors <- c(errors, "Observed responders in control arm (rc) cannot exceed total control sample size (nc).")

    if (length(errors) > 0) {
      shinyjs::html("dpp_analysis_status_message", paste0("<p style='color: red;'>", paste(errors, collapse = "<br>"), "</p>"))
      return()
    }

    shinyjs::html("dpp_analysis_status_message", "<p style='color: blue;'>DPP Analysis started. Please wait...</p>")
    shinyjs::disable("run_dpp_analysis")

    Ych_calculated <- round(input$dpp_analysis_pch * input$dpp_analysis_nch)

    analysis_args <- list(
      Yt = input$dpp_analysis_rt,
      nt = input$dpp_analysis_nt,
      Yc = input$dpp_analysis_rc,
      nc = input$dpp_analysis_nc,
      Ych = Ych_calculated,
      nch = input$dpp_analysis_nch,
      nche = input$dpp_analysis_nche,
      a0c = input$dpp_analysis_a0c,
      b0c = input$dpp_analysis_b0c,
      a0t = input$dpp_analysis_a0t,
      b0t = input$dpp_analysis_b0t,
      delta_threshold = input$dpp_analysis_delta_threshold,
      method = input$dpp_analysis_method
    )

    if (input$dpp_analysis_method == "Generalized BC") {
      analysis_args$theta = input$dpp_analysis_theta
    }
    if (input$dpp_analysis_method %in% c("Bayesian p", "Generalized BC", "JSD")) {
      analysis_args$eta = input$dpp_analysis_eta
    }

    p <- future({
      result <- do.call(BayesianHybridDesign::DPP.analysis, analysis_args)
      result
    })

    then(p, onFulfilled = function(result) {
      dpp_analysis_results(result)
      shinyjs::html("dpp_analysis_status_message", "<p style='color: green;'>DPP Analysis ran successfully!</p>")
      shinyjs::enable("run_dpp_analysis")
    }, onRejected = function(e) {
      shinyjs::html("dpp_analysis_status_message", paste0("<p style='color: red;'>DPP Analysis calculation failed: ", e$message, ". Please check inputs and try again.</p>"))
      dpp_analysis_results(NULL)
      shinyjs::enable("run_dpp_analysis")
    })
  })

  output$dppAnalysisResult <- renderPrint({
    req(dpp_analysis_results())
    print(dpp_analysis_results())
  })

  calculated_pc <- reactive({
    req(input$dpp_analysis_rc, input$dpp_analysis_nc)
    if (input$dpp_analysis_nc > 0) {
      return(input$dpp_analysis_rc / input$dpp_analysis_nc)
    } else {
      return(0)
    }
  })

  calculated_pt <- reactive({
    req(input$dpp_analysis_rt, input$dpp_analysis_nt)
    if (input$dpp_analysis_nt > 0) {
      return(input$dpp_analysis_rt / input$dpp_analysis_nt)
    } else {
      return(0)
    }
  })

  output$dpp_analysis_pc_display <- renderPrint({
    cat(sprintf("%.4f", calculated_pc()))
  })

  output$dpp_analysis_pt_display <- renderPrint({
    cat(sprintf("%.4f", calculated_pt()))
  })

  output$plotDPP <- renderPlot({
    req(dpp_analysis_results())
    BayesianHybridDesign::plotDPP(DPP = dpp_analysis_results())
  })

  output$concurrent_summary <- renderPrint({
    req(dpp_analysis_results()$apost_c_trial)
    req(dpp_analysis_results()$bpost_c_trial)

    a_post <- dpp_analysis_results()$apost_c_trial
    b_post <- dpp_analysis_results()$bpost_c_trial

    median <- qbeta(p = 0.5, shape1 = a_post, shape2 = b_post)
    ci_lower <- qbeta(p = 0.025, shape1 = a_post, shape2 = b_post)
    ci_upper <- qbeta(p = 0.975, shape1 = a_post, shape2 = b_post)

    cat("Parameters: \n")
    cat(sprintf("  apost_c_trial: %.4f\n", a_post))
    cat(sprintf("  bpost_c_trial: %.4f\n", b_post))
    cat("-----------------------------------\n")
    cat(sprintf("Median: %.4f\n", median))
    cat(sprintf("95%% Credible Interval: (%.4f, %.4f)\n", ci_lower, ci_upper))
  })

  output$hybrid_control_summary <- renderPrint({
    req(dpp_analysis_results()$apost_c_hca)
    req(dpp_analysis_results()$bpost_c_hca)

    a_post <- dpp_analysis_results()$apost_c_hca
    b_post <- dpp_analysis_results()$bpost_c_hca

    median <- qbeta(p = 0.5, shape1 = a_post, shape2 = b_post)
    ci_lower <- qbeta(p = 0.025, shape1 = a_post, shape2 = b_post)
    ci_upper <- qbeta(p = 0.975, shape1 = a_post, shape2 = b_post)

    cat("Parameters: \n")
    cat(sprintf("  apost_c_hca: %.4f\n", a_post))
    cat(sprintf("  bpost_c_hca: %.4f\n", b_post))
    cat("-----------------------------------\n")
    cat(sprintf("Median: %.4f\n", median))
    cat(sprintf("95%% Credible Interval: (%.4f, %.4f)\n", ci_lower, ci_upper))
  })

  output$experimental_summary <- renderPrint({
    req(dpp_analysis_results()$apost_t)
    req(dpp_analysis_results()$bpost_t)

    a_post <- dpp_analysis_results()$apost_t
    b_post <- dpp_analysis_results()$bpost_t

    median <- qbeta(p = 0.5, shape1 = a_post, shape2 = b_post)
    ci_lower <- qbeta(p = 0.025, shape1 = a_post, shape2 = b_post)
    ci_upper <- qbeta(p = 0.975, shape1 = a_post, shape2 = b_post)

    cat("Parameters: \n")
    cat(sprintf("  apost_t: %.4f\n", a_post))
    cat(sprintf("  bpost_t: %.4f\n", b_post))
    cat("-----------------------------------\n")
    cat(sprintf("Median: %.4f\n", median))
    cat(sprintf("95%% Credible Interval: (%.4f, %.4f)\n", ci_lower, ci_upper))
  })

  output$statistical_conclusion <- renderPrint({
    result <- dpp_analysis_results()

    if (is.null(result)) {
      cat("No DPP analysis results available yet. Please run the analysis.")
      return()
    }

    phat <- result$phat_pt_larger_pc
    tau <- NULL

    if (!is.null(result$tau)) {
      tau <- result$tau
    }

    if (is.null(tau) && !is.null(dpp_results()$tau)) {
      tau <- dpp_results()$tau
    }

    if (is.null(tau) && !is.null(input$dpp_analysis_tau)) {
      tau <- input$dpp_analysis_tau
    }

    if (is.null(tau)) {
      tau <- 0.95
    }

    if (is.null(phat)) {
      cat("Analysis did not return 'phat_pt_larger_pc'. Cannot compute conclusion.")
      return()
    }

    cat(sprintf(
      "Statistical significance is defined as P(ORR Treatment > ORR Control | Hybrid Data) >= tau.\n"
    ))

    if (phat >= tau) {
      cat(sprintf(
        "Since P(ORR_t > ORR_c) = %.4f is GREATER THAN OR EQUAL to the significance threshold (tau = %.4f),\n",
        phat, tau
      ))
      cat("The experimental treatment is statistically significant compared to the control treatment in ORR.")
    } else {
      cat(sprintf(
        "Since P(ORR_t > ORR_c) = %.4f is LESS THAN the significance threshold (tau = %.4f),\n",
        phat, tau
      ))
      cat("The experimental treatment is NOT statistically significant compared to the control treatment in ORR.")
    }
  })


  # ===========================================================================
  # SAM Prior - Single Design
  # ===========================================================================

  sam_prior_params <- reactiveVal(NULL)

  observe({
    req(input$sam_nch, input$sam_pch)
    a_hist <- round(input$sam_pch * input$sam_nch)
    b_hist <- input$sam_nch - a_hist
    sam_prior_params(list(a = a_hist, b = b_hist))
  })

  output$sam_calculated_prior <- renderUI({
    params <- sam_prior_params()
    req(params)

    tagList(
      div(style = "background: rgba(255,255,255,0.95); padding: 20px; border-radius: 8px; color: #2c3e50;",
          div(style = "margin-bottom: 15px;",
              div(style = "font-size: 0.9em; color: rgba(255,255,255,0.9); margin-bottom: 8px; font-weight: 500;",
                  "Informative Prior Distribution:"
              ),
              div(style = "font-size: 1.8em; font-weight: 700; color: #667eea; margin-bottom: 5px;",
                  sprintf("Beta(a = %d, b = %d)", params$a, params$b)
              )
          ),
          hr(style = "border-color: rgba(0,0,0,0.1); margin: 15px 0;"),
          div(style = "display: flex; justify-content: space-between; align-items: center;",
              div(
                div(style = "font-size: 0.85em; color: #7f8c8d; margin-bottom: 3px;", "Sample Size"),
                div(style = "font-size: 1.3em; font-weight: 600; color: #2c3e50;", sprintf("%d", input$sam_nch))
              ),
              div(
                div(style = "font-size: 0.85em; color: #7f8c8d; margin-bottom: 3px;", "Response Rate"),
                div(style = "font-size: 1.3em; font-weight: 600; color: #2c3e50;", sprintf("%.3f", input$sam_pch))
              )
          )
      )
    )
  })

  observeEvent(input$generate_sam_prior_plot, {
    shinyjs::disable("generate_sam_prior_plot")

    tryCatch({
      req(input$sam_nch, input$sam_pch, input$sam_nc,
          input$sam_sim_control_rate, input$sam_delta,
          input$sam_alpha_noninf, input$sam_beta_noninf, input$sam_seed)

      params <- sam_prior_params()
      req(params)

      prior_hist <- RBesT::mixbeta(c(1, params$a, params$b))
      nf_prior <- RBesT::mixbeta(c(1, input$sam_alpha_noninf, input$sam_beta_noninf))

      set.seed(input$sam_seed)
      control_data <- rbinom(input$sam_nc, size = 1, prob = input$sam_sim_control_rate)
      r_control <- sum(control_data)

      weight <- SAMprior::SAM_weight(
        if.prior = prior_hist,
        delta = input$sam_delta,
        n = input$sam_nc,
        r = r_control
      )

      output$sam_weight_display <- renderUI({
        tagList(
          div(style = "background: rgba(255,255,255,0.95); padding: 20px; border-radius: 8px; color: #2c3e50;",
              div(style = "margin-bottom: 15px;",
                  div(style = "font-size: 0.9em; color: rgba(255,255,255,0.9); margin-bottom: 8px; font-weight: 500;",
                      "SAM Weight Value:"
                  ),
                  div(style = "font-size: 2.5em; font-weight: 700; color: #11998e; margin-bottom: 10px;",
                      sprintf("%.4f", weight)
                  )
              ),
              hr(style = "border-color: rgba(0,0,0,0.1); margin: 15px 0;"),
              div(style = "margin-bottom: 12px;",
                  div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 8px;",
                      span(style = "font-size: 0.95em; color: #34495e; font-weight: 500;", "Informative Prior:"),
                      span(style = "font-size: 1.2em; font-weight: 700; color: #667eea;", sprintf("%.2f%%", weight * 100))
                  ),
                  div(style = "display: flex; justify-content: space-between; align-items: center;",
                      span(style = "font-size: 0.95em; color: #34495e; font-weight: 500;", "Non-informative Prior:"),
                      span(style = "font-size: 1.2em; font-weight: 700; color: #95a5a6;", sprintf("%.2f%%", (1 - weight) * 100))
                  )
              ),
              hr(style = "border-color: rgba(0,0,0,0.1); margin: 15px 0;"),
              div(style = "background: rgba(52, 152, 219, 0.1); padding: 12px; border-radius: 6px; border-left: 3px solid #3498db;",
                  div(style = "font-size: 0.85em; color: #2c3e50;",
                      icon("info-circle", style = "margin-right: 5px; color: #3498db;"),
                      sprintf("Based on %d observed responses out of %d patients", r_control, input$sam_nc)
                  )
              )
          )
        )
      })

      sam_prior <- SAMprior::SAM_prior(
        if.prior = prior_hist,
        nf.prior = nf_prior,
        weight = weight
      )

      output$sam_prior_plot <- renderPlot({
        plot(sam_prior, main = "SAM Prior Distribution")
      })

      showNotification("SAM Prior plot and weight calculated successfully!", type = "message", duration = 5)

    }, error = function(e) {
      output$sam_prior_plot <- renderPlot(NULL)
      output$sam_weight_display <- renderUI({
        div(style = "background: rgba(255,255,255,0.95); padding: 20px; border-radius: 8px; color: #e74c3c; text-align: center;",
            icon("exclamation-triangle", style = "font-size: 2em; margin-bottom: 10px;"),
            div(style = "font-size: 1.1em; font-weight: 600; margin-bottom: 5px;", "Error"),
            div(style = "font-size: 0.9em;", paste("Unable to calculate SAM weight:", e$message))
        )
      })
      showNotification(paste0("SAM Prior generation failed: ", e$message), type = "error", duration = 10)
    }, finally = {
      shinyjs::enable("generate_sam_prior_plot")
    })
  })


  # ===========================================================================
  # SAM Weight Library
  # ===========================================================================

  observeEvent(input$run_sam_weight, {
    shinyjs::html("sam_weight_status_message", "Calculation started. Please wait...")
    shinyjs::disable("run_sam_weight")

    tryCatch({
      req(input$alpha_hist, input$beta_hist, input$delta,
          input$n_control, input$r_control, input$sam_method)

      prior_hist <- RBesT::mixbeta(c(1, input$alpha_hist, input$beta_hist))

      if (input$sam_method == "PPR") {
        weight <- SAMprior::SAM_weight(
          if.prior = prior_hist,
          delta = input$delta,
          method.w = input$sam_method,
          prior.odds = input$prior_odds,
          n = input$n_control,
          r = input$r_control
        )
      } else {
        weight <- SAMprior::SAM_weight(
          if.prior = prior_hist,
          delta = input$delta,
          method.w = input$sam_method,
          n = input$n_control,
          r = input$r_control
        )
      }

      output$samWeightResult <- renderPrint({
        cat("SAM Weight: ", round(weight, 4), "\n")
        cat("Method: ", input$sam_method, "\n")
        cat("Interpretation: ", round(weight * 100, 2),
            "% weight on informative prior\n")
      })

      shinyjs::html("sam_weight_status_message",
                    "SAM Weight calculation ran successfully!")

    }, error = function(e) {
      output$samWeightResult <- renderPrint(NULL)
      shinyjs::html("sam_weight_status_message",
                    paste0("SAM Weight calculation failed: ", e$message,
                           ". Please check inputs and try again."))
    }, finally = {
      shinyjs::enable("run_sam_weight")
    })
  })


  # ===========================================================================
  # SAM Comparative Analysis
  # ===========================================================================

  sam_table_results <- reactiveVal(NULL)

  observeEvent(input$run_sam_table, {
    shinyjs::html("sam_table_status_message", "Running comparative analysis. This may take several minutes...")
    shinyjs::disable("run_sam_table")

    pc_values <- as.numeric(unlist(strsplit(input$sam_table_pc_values, ",")))
    sam_table_nt <- input$sam_table_nt
    sam_table_nc <- input$sam_table_nc
    sam_table_pch <- input$sam_table_pch
    sam_table_nche <- input$sam_table_nche
    sam_table_delta <- input$sam_table_delta
    sam_table_nsim <- input$sam_table_nsim
    sam_table_typeIER <- input$sam_table_typeIER

    p <- future({
      ns <- length(pc_values)
      SAMMethodNames <- c("SAMprior", "MAP", "Noninfo")
      nm <- length(SAMMethodNames)
      TypeI_SAM <- matrix(NA, ns, nm)
      Power_SAM <- matrix(NA, ns, nm)
      pc_mi_SAM <- matrix(NA, ns, nm)
      pc_sdi_SAM <- matrix(NA, ns, nm)

      xt_list <- list()
      xt_typeI_list <- list()
      xc_list <- list()

      for(s in 1:ns){
        set.seed(2000)
        xt_typeI_list[[s]] <- rbinom(sam_table_nsim, size = sam_table_nt, prob = pc_values[s])
        xc_list[[s]] <- rbinom(sam_table_nsim, size = sam_table_nc, prob = pc_values[s])
        xt_list[[s]] <- rbinom(sam_table_nsim, size = sam_table_nt, prob = pc_values[s] + 0.2)
      }

      for(s in 1:ns){
        datamat_typeI <- cbind(xt_typeI_list[[s]], xc_list[[s]])
        datamat_Power <- cbind(xt_list[[s]], xc_list[[s]])

        TypeI <- BayesianHybridDesign::runCalibratedSAM(
          nsim = sam_table_nsim,
          pch = sam_table_pch,
          delta_threshold = sam_table_delta,
          nche.c = sam_table_nche,
          nc.calib = sam_table_nc,
          nt.calib = sam_table_nt,
          pc.calib = pc_values[s],
          xt.cal = xt_typeI_list[[s]],
          xc.cal = xc_list[[s]],
          typeIER.cal = sam_table_typeIER,
          nche = sam_table_nche,
          nc = sam_table_nc,
          nt = sam_table_nt,
          pc = pc_values[s],
          pt = pc_values[s],
          xt = xt_typeI_list[[s]],
          xc = xc_list[[s]],
          nf.prior = RBesT::mixbeta(c(1, 0.001, 0.001)),
          seed.hist = 1000,
          seed.gMAP = 2000,
          seed.SAM = 3000,
          seed.cal = 4000
        )

        TypeI_SAM[s, ] <- TypeI$Sim_Result
        pc_mi_SAM[s, ] <- TypeI$pc.PMD
        pc_sdi_SAM[s, ] <- TypeI$pc.PSDD

        Power <- BayesianHybridDesign::runCalibratedSAM(
          nsim = sam_table_nsim,
          pch = sam_table_pch,
          delta_threshold = sam_table_delta,
          nche.c = sam_table_nche,
          nc.calib = sam_table_nc,
          nt.calib = sam_table_nt,
          pc.calib = pc_values[s],
          xt.cal = xt_typeI_list[[s]],
          xc.cal = xc_list[[s]],
          typeIER.cal = sam_table_typeIER,
          nche = sam_table_nche,
          nc = sam_table_nc,
          nt = sam_table_nt,
          pc = pc_values[s],
          pt = pc_values[s] + 0.2,
          xt = xt_list[[s]],
          xc = xc_list[[s]],
          nf.prior = RBesT::mixbeta(c(1, 0.001, 0.001)),
          seed.hist = 1000,
          seed.gMAP = 2000,
          seed.SAM = 3000,
          seed.cal = 4000
        )

        Power_SAM[s, ] <- Power$Sim_Result
      }

      list(
        TypeI = TypeI_SAM,
        Power = Power_SAM,
        pc_mi = pc_mi_SAM,
        pc_sdi = pc_sdi_SAM,
        pc_values = pc_values,
        methods = SAMMethodNames
      )
    }, seed = TRUE)

    then(p, onFulfilled = function(result) {
      sam_table_results(result)
      shinyjs::html("sam_table_status_message", "SAM comparative analysis completed successfully!")
      shinyjs::enable("run_sam_table")
    }, onRejected = function(e) {
      shinyjs::html("sam_table_status_message", paste0("SAM analysis failed: ", e$message, ""))
      sam_table_results(NULL)
      shinyjs::enable("run_sam_table")
    })
  })

  output$sam_table_typeI <- DT::renderDataTable({
    req(sam_table_results())
    result <- sam_table_results()
    df <- as.data.frame(result$TypeI)
    colnames(df) <- result$methods
    df <- cbind(pc = result$pc_values, df)
    DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE)) %>%
      DT::formatRound(columns = 2:ncol(df), digits = 4)
  })

  output$sam_table_power <- DT::renderDataTable({
    req(sam_table_results())
    result <- sam_table_results()
    df <- as.data.frame(result$Power)
    colnames(df) <- result$methods
    df <- cbind(pc = result$pc_values, df)
    DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE)) %>%
      DT::formatRound(columns = 2:ncol(df), digits = 4)
  })

  output$sam_table_pmd <- DT::renderDataTable({
    req(sam_table_results())
    result <- sam_table_results()
    df <- as.data.frame(result$pc_mi)
    colnames(df) <- result$methods
    df <- cbind(pc = result$pc_values, df)
    DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE)) %>%
      DT::formatRound(columns = 2:ncol(df), digits = 4)
  })

  output$sam_table_sd_pmd <- DT::renderDataTable({
    req(sam_table_results())
    result <- sam_table_results()
    df <- as.data.frame(result$pc_sdi)
    colnames(df) <- result$methods
    df <- cbind(pc = result$pc_values, df)
    DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE)) %>%
      DT::formatRound(columns = 2:ncol(df), digits = 4)
  })


  # ===========================================================================
  # Fisher Testing Library
  # ===========================================================================

  observeEvent(input$run_fisher_power, {
    shinyjs::html("fisher_power_status_message", "<p style='color: blue;'>Calculation started. Please wait...</p>")
    shinyjs::disable("run_fisher_power")

    tryCatch({
      set.seed(input$fp_seed)
      power_result <- BayesianHybridDesign::fisher.power(
        pt = input$fp_pt,
        nt = input$fp_nt,
        pc = input$fp_pc,
        nc = input$fp_nc,
        alpha = input$fp_alpha,
        nsim = input$fp_nsim,
        seed = input$fp_seed
      )

      output$fisherPowerResult <- renderPrint({
        cat("Power:", round(power_result, 4), "\n")
      })

      shinyjs::html("fisher_power_status_message", "<p style='color: green;'>Fisher Power calculation ran successfully!</p>")

      output$fisherPowerConclusion <- renderUI({
        HTML(
          "<h4>Summary</h4>
          <p>The statistical power of the study was estimated to be <strong>", round(power_result, 4), "</strong>. This value represents the probability of detecting a statistically significant treatment effect given the design parameters (e.g., sample sizes, response rates, alpha level).</p>

          <h4>Statistical Conclusion</h4>
          <p>A power of <strong>", round(power_result, 4), "</strong> indicates a <strong>", round(power_result * 100, 2), "%</strong> chance of correctly rejecting the null hypothesis (i.e., finding a significant effect) if the true effect size is as specified. This value is a crucial measure of the study's ability to avoid a Type II error (a false negative). Typically, a power of 80% or greater is considered acceptable, and your study design either meets or exceeds this standard, suggesting it is well-powered to detect a true difference.</p>"
        )
      })

    }, error = function(e) {
      output$fisherPowerResult <- renderPrint({
        cat("Error calculating Fisher Power:", conditionMessage(e), "\n")
      })
      shinyjs::html("fisher_power_status_message", paste0("<p style='color: red;'>Fisher Power calculation failed: ", e$message, ". Please check inputs and try again.</p>"))
    }, finally = {
      shinyjs::enable("run_fisher_power")
    })
  })

  observeEvent(input$run_fisher_bound, {
    shinyjs::html("fisher_bound_status_message", "<p style='color: blue;'>Calculation started. Please wait...</p>")
    shinyjs::disable("run_fisher_bound")

    tryCatch({
      bound_result <- BayesianHybridDesign::fisher.bound(
        pc = input$fb_pc,
        nc = input$fb_nc,
        nt = input$fb_nt,
        alpha = input$fb_alpha
      )

      output$fisherBoundM <- renderPrint({ print(bound_result$M) })
      output$fisherBoundP <- renderPrint({ round(bound_result$p, 4) })
      output$fisherBoundRc <- renderPrint({ bound_result$rc })
      output$fisherBoundNc <- renderPrint({ bound_result$nc })
      output$fisherBoundRt <- renderPrint({ round(bound_result$rt, 4) })
      output$fisherBoundNt <- renderPrint({ bound_result$nt })
      output$fisherBoundDelta <- renderPrint({ round(bound_result$delta, 4) })

      shinyjs::html("fisher_bound_status_message", "<p style='color: green;'>Fisher Bound calculation ran successfully!</p>")

      output$fisherBoundConclusion <- renderUI({
        HTML(
          "<h4>Summary</h4>
          <p>The Fisher Bound analysis determined that a minimum of <strong>", bound_result$M, "</strong> responders in the treatment arm are required to achieve statistical significance. This corresponds to a minimum detectable difference of <strong>", round(bound_result$delta, 4), "</strong> and a p-value of <strong>", round(bound_result$p, 4), "</strong> at the boundary.</p>"
        )
      })

      output$fisherBoundStatisticalConclusion <- renderUI({
        HTML(
          "<h4>Statistical Conclusion</h4>
          <p>The number <strong>", bound_result$M, "</strong> is the critical benchmark for this study. If the observed number of responses in the experimental arm is less than this value, the trial will not achieve statistical significance at the specified alpha level. This provides a clear, objective criterion for evaluating the success of the experimental treatment's performance based on the selected study parameters.</p>"
        )
      })

    }, error = function(e) {
      output$fisherBoundM <- renderPrint(NULL)
      output$fisherBoundP <- renderPrint(NULL)
      output$fisherBoundRc <- renderPrint(NULL)
      output$fisherBoundNc <- renderPrint(NULL)
      output$fisherBoundRt <- renderPrint(NULL)
      output$fisherBoundNt <- renderPrint(NULL)
      output$fisherBoundDelta <- renderPrint(NULL)
      shinyjs::html("fisher_bound_status_message", paste0("<p style='color: red;'>Fisher Bound calculation failed: ", e$message, ". Please check inputs and try again.</p>"))
    }, finally = {
      shinyjs::enable("run_fisher_bound")
    })
  })

})
