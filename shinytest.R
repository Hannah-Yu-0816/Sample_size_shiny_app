library(shiny)
library(shinydashboard)
library(causalTree)
library(purrr)
library(grf)
library(dplyr)
library(tidyverse)
library(parallel)
library(iterators)
library(doParallel)
library(foreach)
library(ggplot2)
library(plotly)
library(DT)

# ======================== Function ========================

generate_data <- function(seed = 42, n = 1000, k = 1, d = 2, x_prob = NULL, beta = NULL, gamma = NULL, noise_sd = 0.1) {
  
  # Input parameter
  # n: sample size
  # k: number of effect modifiers
  # d: number of total covariates (the first k are effect modifiers)
  # x_prob: the probability that each covariate equal to 1 (parameters for the Bernoulli distribution)
  # beta: Coefficients controlling the effect of covariates on the baseline risk (including intercept)
  # gamma: Coefficients controlling the effect of effect modifiers on the treatment effect (including main effect)
  # noise_sd: Standard deviation of the Gaussian noise term
  # seed: Random seed for reproducibility
  
  # Set random seed for reproducibility
  set.seed(seed)
  
  # Input check and eror print 
  
  # set default value for x_prob and error check
  if (is.null(x_prob)) x_prob <- rep(0.5, d)
  
  if (!is.numeric(x_prob)) stop("Error: The 'x_prob' parameter must be a numeric vector.")
  
  if (length(x_prob) != d) stop("Error: The length of 'x_prob' must equal the total number of covariates.")
  
  # set default value for beta and error check
  if (is.null(beta)) beta <- rep(1, d+1)
  
  if (!is.numeric(beta)) stop("Error: The 'beta' parameter must be a numeric vector.")
  
  if (length(beta) != d+1) stop("Error: The length of 'beta' must equal the total number of covariates + 1 (main baseline risk).")
  
  # set default value for gamma and error check
  # if (is.null(gamma)) gamma <- rep(1, k+1)
  # 
  # if (!is.numeric(gamma)) stop("Error: The 'gamma' parameter must be a numeric vector.")
  # 
  # if (length(gamma) != k+1) stop("Error: The length of 'gamma' must equal the number of effect modifiers + 1 (main effect).")
  
  # error check
  if (d < k) stop("Error: The total number of the covariates should be larger than the number of effect modifiers.")
  
  if(k > 5) stop("Error: The maximum number of the effect modifier is 5.")
  
  # Generate n observations of d binary covariates
  X <- sapply(x_prob, function(p) rbinom(n, size = 1, prob = p))
  colnames(X) <- paste0("X", 1:d)
  
  # Generate phenotype      
  phenotype <- as.integer(apply(X[, 1:k, drop = FALSE], 1, prod) == 1) 
  
  # Random treatment assignment A (Bernoulli(0.5))
  A <- rbinom(n, 1, 0.5)
  
  # Random noise (epsilon), normally distributed
  epsilon <- rnorm(n, mean = 0, sd = noise_sd)
  
  # Compute potential outcomes
  Y1_linear <- beta[1] + X %*% beta[2:(d+1)] + gamma[1] + apply(X[, 1:k, drop = FALSE], 1, prod) * gamma[2] + epsilon
  
  # Y1_linear <- beta[1] + X %*% beta[2:(d+1)] + gamma[1] + X[, 1:k, drop = FALSE] %*% gamma[2:(k+1)] + epsilon
  
  Y0_linear <- if (d-k > 0) beta[1] + X %*% beta[2:(d+1)] + epsilon else epsilon
  
  p1 <- plogis(Y1_linear)
  p0 <- plogis(Y0_linear)
  
  # Observed outcome Y depends on treatment assignment A
  Y1 <- rbinom(n, 1, p1)
  Y0 <- rbinom(n, 1, p0)
  #Y <- A * Y1 + (1-A) *Y0
  Y <- ifelse(A==1, Y1, Y0)
  
  # Return the dataset as a data frame
  data <- data.frame(
    ID = 1:n,
    A = A,
    Y = Y,
    #Y1 = Y1,
    #Y0 = Y0,
    CATE = p1-p0
  )
  
  data <- cbind(data, as.data.frame(X), as.data.frame(phenotype))
  
  return(data)
}


check_trees_purity <- function(forest, data, num_trees, mode = c("x","phenotype")) {
  mode <- match.arg(mode)
  
  success <- 0
  unsuccess <- 0
  partial_success <- 0
  
  for(i in seq_len(num_trees)) {
    
    leaf_list <- forest$`_leaf_samples`[[i]]
    leaf_nodes <- lapply(leaf_list, function(x) length(x) > 0)
    leaf_node_indices <- leaf_list[unlist(leaf_nodes)]
    
    if (mode == "x"){
      
      x1_props <- sapply(leaf_node_indices, function(leaf) mean(data$X1[leaf + 1]))
      
      if (all(x1_props %in% c(0,1))){
        success <- success + 1
        #next
      } else if (any(x1_props == 1) && any(x1_props > 0 & x1_props < 1)) {
        partial_success <- partial_success + 1
      } else{
        unsuccess <- unsuccess + 1
      }
      
    } else if(mode == "phenotype"){
      
      # number of samples with phenotype == 1
      all_leaf_samples <- unlist(leaf_node_indices) + 1
      n_current_tree_all_pheno <- sum(data[all_leaf_samples, ]$phenotype)
      
      # if no phenotype==1 sample, skip this tree
      if (n_current_tree_all_pheno == 0) {
        unsuccess <- unsuccess + 1
        next
      }
      
      is_full_phenotype <- vapply(leaf_node_indices, function(leaf) {
        sample_indices <- leaf + 1
        all(data[sample_indices,]$phenotype == 1)
      }, logical(1))
      
      phenotype_covered <- sum(vapply(leaf_node_indices[is_full_phenotype], length, integer(1)))
      prop_phenotype <- phenotype_covered / n_current_tree_all_pheno
      
      if (prop_phenotype == 1) {
        success <- success + 1
      } else if (prop_phenotype > 0) {
        partial_success <- partial_success + 1
      } else {
        unsuccess <- unsuccess + 1
      }
      
    }
    
  }
  
  return(list(
    success_ratio = success / num_trees,
    partial_success_ratio = partial_success / num_trees,
    unsuccess_ratio = unsuccess / num_trees
  ))
}


forest_built <- function(data, 
                         num_trees, 
                         k, 
                         sample.fraction = 0.5, 
                         honesty = TRUE, 
                         honesty.fraction = 0.5){
  
  # Train Causal Forest
  if(honesty == FALSE){
    causal_forest_model <- causal_forest(data %>% select(starts_with("X")), 
                                         data$Y, 
                                         data$A,
                                         num.trees = num_trees,
                                         sample.fraction = sample.fraction,
                                         honesty = honesty,
                                         seed = 123)
  }else{
    causal_forest_model <- causal_forest(data %>% select(starts_with("X")), 
                                         data$Y, 
                                         data$A,
                                         num.trees = num_trees,
                                         sample.fraction = sample.fraction,
                                         honesty = honesty,
                                         honesty.fraction = honesty.fraction,
                                         seed = 123)
  }
  
  # predicted cate
  tau.hat.est <- causal_forest_model$predictions
  
  # heterogeneity captured
  pheno_data_index <- data$phenotype == 1
  
  truth_cate <- mean(data$CATE[pheno_data_index]) - mean(data$CATE[!pheno_data_index])
  
  cate_phenotype_1 <- mean(tau.hat.est[pheno_data_index])
  cate_phenotype_0 <- mean(tau.hat.est[!pheno_data_index])
  heterogeneity_captured <- ((cate_phenotype_1 - cate_phenotype_0) / truth_cate)
  
  # only one effect modifier
  if(k==1){
    tmp_result <- check_trees_purity(causal_forest_model, data, num_trees, mode="x")
  }else{ # more than 1 effect modifiers
    tmp_result <- check_trees_purity(causal_forest_model, data, num_trees, mode="phenotype")
  }
  
  success_ratio <- tmp_result[[1]]
  partial_success_ratio <- tmp_result[[2]]
  unsuccess_ratio <- tmp_result[[3]]
  
  return(list(heterogeneity_captured = heterogeneity_captured,
              success_ratio = success_ratio,
              partial_success_ratio = partial_success_ratio,
              unsuccess_ratio = unsuccess_ratio))
  
}


run_simulation <- function(sim_no = 500, 
                           sample_param = NULL, 
                           param = NULL,
                           k = 1, 
                           d = 2, 
                           x_prob, 
                           beta, 
                           gamma, 
                           num_trees = 2000){
  
  is_tune_mode <- !is.null(param)
  is_sample_mode <- !is_tune_mode
  
  if (is_tune_mode && length(sample_param) != 1) {
    stop("In tuning mode, sample size must be a single fixed number.")
  }
  
  if (is_sample_mode) {
    param_grid <- data.frame(sample_size = sample_param)
  } else {
    param_grid <- expand.grid(param)
  }
  
  n_param_grid <- nrow(param_grid)
  
  heterogeneity <- numeric(n_param_grid)
  success <- numeric(n_param_grid)
  partial_success <- numeric(n_param_grid)
  unsuccess <- numeric(n_param_grid)
  
  # Detect cores and set up parallel backend
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  for(i in seq_len(n_param_grid)){
    start_time <- Sys.time()
    set.seed(1)
    seeds <- sample(1:(sim_no*2), sim_no, replace = FALSE)
    
    # Parallel loop over simulations
    results <- foreach(j = 1:sim_no,
                       .combine = rbind,
                       .packages = c("grf", "dplyr"),
                       .export = c("generate_data", 
                                   "forest_built",
                                   "check_trees_purity")) %dopar% {
                                     
                                     sim_data <- generate_data(
                                       seed = seeds[j], 
                                       n = if (is_sample_mode) param_grid$sample_size[i] else sample_param, 
                                       k = k, 
                                       d = d, 
                                       x_prob = x_prob, 
                                       beta = beta, 
                                       gamma = gamma
                                     )
                                     
                                     # Build forest
                                     result <- forest_built(
                                       data = sim_data, 
                                       num_trees = num_trees, 
                                       k = k,
                                       sample.fraction = if (!is_sample_mode) param_grid[i, 1] else 0.5,
                                       honesty.fraction = if (!is_sample_mode) param_grid[i, 2] else 0.5
                                     )
                                     
                                     # Combine results as a row
                                     c(result$heterogeneity_captured, 
                                       result$success_ratio, 
                                       result$partial_success_ratio, 
                                       result$unsuccess_ratio)
                                   }
    
    
    # Compute the average/median from the parallel results
    heterogeneity[i] <- median(results[, 1])
    success[i] <- median(results[, 2])
    partial_success[i] <- median(results[, 3])
    unsuccess[i] <- median(results[,4])
    
    end_time <- Sys.time()
    if (is_sample_mode) {
      cat("Sample size", param_grid$sample_size[i], ": running time", end_time - start_time, "\n")
    } else {
      cat("Tuning param", paste(names(param_grid), unlist(param_grid[i, ]), sep = "=", collapse = ", "), ": running time", end_time - start_time, "\n")
    }
  }
  
  stopCluster(cl)
  
  return_result <- data.frame(
    heterogeneity_captured = heterogeneity,
    success_ratio = success,
    partial_success_ratio = partial_success,
    unsuccess_ratio = unsuccess
  )
  
  return_result <- cbind(param_grid, return_result)
  return(return_result)
}


# ======================== UI ==============================
ui <- dashboardPage(
  dashboardHeader(title = "Causal Forest Simulation"),
  
  dashboardSidebar(
    width = 400,
    sidebarMenu(
      id = "current_tab",
      menuItem("Simulation Controls", tabName = "controls", icon = icon("sliders-h")),
      menuItem("Parameter Tuning", tabName = "tuning_param", icon = icon("cogs")),
      
      conditionalPanel(
        condition = "input.current_tab == 'controls'",
        
        numericInput("k", "Number of Effect Modifiers:", value = 1, min = 0),
        uiOutput("k_warning"),
        numericInput("m", "Number of Other Risk Factors:", value = 1, min = 0),
        
        # conditionalPanel(
        #   condition = "input.k > 1",
        #   selectInput("model_type", 
        #               label = "Choose Treatment Effect Structure:", 
        #               choices = c("Linear" = "linear", "Phenotype" = "pheno"),
        #               selected = "pheno")
        # ),
        
        numericInput("beta0", 
                     label = HTML("Baseline Risk Intercept \\( \\beta_0 \\):"), 
                     value = 0),
        numericInput("gamma0", 
                     label = HTML("Baseline Treatment Effect Intercept \\( \\gamma_0 \\):"), 
                     value = 0),
        #numericInput("gamma", "Effect Size of Phenotype:", value = 1),
        
        conditionalPanel(
          condition = "input.k == 1",
          numericInput("gamma", 
                       label = HTML("Treatment Effect Size \\( \\gamma \\):"), 
                       value = 1, step = 0.01)
        ),
        conditionalPanel(
          # condition = "input.k > 1 && input.model_type == 'pheno'",
          condition = "input.k > 1",
          numericInput("gamma", 
                       label = HTML("Effect Size of Phenotype \\( \\gamma \\):"), 
                       value = 1, step = 0.01)
        ),
        
        uiOutput("covariate_inputs"),
        
        numericInput("sim_no", "Simulation Iterations:", value = 500, min = 10),
        numericInput("num_trees", "Number of Trees in the Forest:", value = 1000, min = 10),
        textInput("sample_param", "Sample Sizes:", value = "100,200,400,600,800,1000"),
        
        actionButton("run_simulation", "Run Simulation", icon = icon("play"))
      ),
      
      conditionalPanel(
        condition = "input.current_tab == 'tuning_param'",
        numericInput("overall_sample_param", "Overall Sample Size:", value = 600, min = 100),
        sliderInput("sample_fraction", "Sample Fraction:", 
                    min = 0.1, max = 0.5, value = c(0.2,0.4), step = 0.1),
        sliderInput("honesty_fraction", "Honesty Fraction:", 
                    min = 0.3, max = 0.8, value = c(0.4,0.7), step = 0.1),
        actionButton("run_tuning_param", "Run Parameter Tuning Simulation", icon = icon("play"))
      )
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .well {
          background-color: black !important;  
          color: white !important;             
          border: 1px solid #444 !important;   
        }
      "))
    ),
    tabItems(
      tabItem(tabName = "controls",
              fluidRow(
                box(
                  title = "Simulation Data Generation Formula and Example Data",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  
                  # LaTeX-rendered formula
                  uiOutput("formula_text"),
                  br(),
                  
                  # Preview of example data
                  DTOutput("example_data")
                ),
                
                box(
                  title = "Simulation Progress", 
                  status = "primary", 
                  solidHeader = TRUE,
                  width = 12,
                  verbatimTextOutput("progress_status")
                ),
                
                box(
                  title = "Simulation Results", 
                  status = "primary", 
                  solidHeader = TRUE,
                  width = 12,
                  plotlyOutput("result_plot", height = 600)
                )
              )
      ),
      tabItem(tabName = "tuning_param",
              fluidRow(
                box(
                  title = "Simulation Progress", 
                  status = "primary", 
                  solidHeader = TRUE,
                  width = 12,
                  verbatimTextOutput("tuning_progress_status")
                ),
                box(
                  title = "Tuning Results: Heterogeneity Captured", 
                  status = "primary", 
                  solidHeader = TRUE,
                  width = 12,
                  plotOutput("tuning_plot", height = 500)
                ),
                box(
                  title = "Tuning Results: Success Ratio of Pure Effect Modifier Leaf", 
                  status = "primary", 
                  solidHeader = TRUE,
                  width = 12,
                  plotOutput("tuning_plot2", height = 500)
                )
              )
      )
    )
  )
)

# ======================== SERVER ===========================
server <- function(input, output, session) {
  
  # ======================== Panel 1 - Run Simulation ===========================
  
  # Reactive value to store simulation progress logs
  progress_text <- reactiveVal("")
  
  # Render progress logs in the Simulation Progress panel
  output$progress_status <- renderText({
    progress_text()
  })
  
  # Set the limitation of number of effect modifiers
  output$k_warning <- renderUI({
    req(input$k)
    if (input$k > 5) {
      HTML("<span style='color: red;'The maximum number of effect modifiers is 5.</span>")
    } else {NULL}
  })
  # observeEvent(input$k, {
  #   if (!is.null(input$k) && input$k > 5) {
  #     updateNumericInput(session, "k", value = 5)
  #     # showNotification("The maximum number of effect modifiers is 5.", type = "warning")
  #   }
  # })
  
  # covariate inputs
  output$covariate_inputs <- renderUI({
    req(input$m)
    req(input$k)
    
    ui_list <- lapply(1:(input$m + input$k), function(i) {
      
      box_title <- if (i <= input$k) {
        paste0("Effect Modifier X", i)
      } else {
        paste0("Other Risk Factor X", i)
      }
      
      withMathJax(
        wellPanel(
          strong(box_title),
          # checkboxInput(paste0("effmod_", i), "Effect Modifier", value = FALSE),
          # conditionalPanel(
          #   condition = sprintf("input.k > 1 && input.effmod_%d == true", i),
          #   checkboxInput(paste0("pheno_", i), "Phenotype", value = FALSE)
          # ),
          numericInput(paste0("beta_", i), 
                       label = paste0("Baseline Covariate Effect \\( \\beta_{", i, "} \\):"), 
                       value = 0, step = 0.01),
          # conditionalPanel(
          #   # condition = sprintf("input.k > 1 && input.model_type == 'linear' && %d <= input.k", i),
          #   numericInput(paste0("gamma_", i), "Treatment Effect Size:", value = 0)
          # ),
          numericInput(paste0("xprob_", i), "Probability:", value = 0.5, min = 0, max = 1, step = 0.01)
        )
      )
      
    })
    
    do.call(tagList, ui_list)
  })
  
  # run simulation button
  observeEvent(input$run_simulation, {
    
    # Extract all user inputs
    d_val <- input$m + input$k
    k_val <- input$k
    sim_no_val <- input$sim_no
    num_trees_val <- input$num_trees
    # model_type_val <- input$model_type
    
    beta0 <- input$beta0
    beta_vals <- sapply(1:d_val, function(i) input[[paste0("beta_", i)]])
    beta_vals <- c(beta0, beta_vals)
    
    gamma0 <- input$gamma0
    gamma <- input$gamma
    gamma_vals <- c(gamma0, gamma)
    # if (k_val == 1 | model_type_val == "pheno"){
    #   gamma <- input$gamma
    #   gamma_vals <- c(gamma0, gamma)
    # }else{
    #   gamma_vals <- sapply(1:k_val, function(i) input[[paste0("gamma_", i)]])
    #   gamma_vals <- c(gamma0, gamma_vals)
    # }
    
    xprob_vals <- sapply(1:d_val, function(i) input[[paste0("xprob_", i)]])
    sample_param_vals <- as.numeric(unlist(strsplit(input$sample_param, ",")))
    # pheno_vals <- sapply(1:k_val, function(i) {
    #   if (!is.null(input[[paste0("pheno_", i)]]) && input[[paste0("pheno_", i)]]) {
    #     return(paste0("X", i))
    #   } else {
    #     return(NULL)
    #   }
    # })
    # pheno_vals <- pheno_vals[!sapply(pheno_vals, is.null)]
    
    # Input validation
    # if (d_val < k_val) {
    #   showNotification("d must be >= k.", type = "error")
    #   return(NULL)
    # }
    
    # Initialize progress text
    progress_text("Starting simulation...\n")
    
    withProgress(message = 'Running Simulation...', value = 0, {
      
      is_sample_mode <- TRUE
      
      # if (is_tune_mode && length(sample_param) != 1) {
      #   stop("In tuning mode, sample size must be a single fixed number.")
      # }
      
      if (is_sample_mode) {
        param_grid <- data.frame(sample_size = sample_param_vals)
      } else {
        param_grid <- expand.grid(param)
      }
      
      n_param_grid <- nrow(param_grid)
      
      heterogeneity <- numeric(n_param_grid)
      success <- numeric(n_param_grid)
      partial_success <- numeric(n_param_grid)
      unsuccess <- numeric(n_param_grid)
      
      # Detect cores and set up parallel backend
      num_cores <- detectCores() - 1
      cl <- makeCluster(num_cores)
      registerDoParallel(cl)
      
      time_logs <- c()
      
      # Loop over each sample size
      for (i in seq_len(n_param_grid)) {
        
        sample_size_now <- param_grid$sample_size[i]
        
        start_time <- Sys.time()
        set.seed(1)
        seeds <- sample(1:(sim_no_val*2), sim_no_val, replace = FALSE)
        
        # Parallel simulation loop for each repetition
        results <- foreach(j = 1:sim_no_val,
                           .combine = rbind,
                           .packages = c("grf", "dplyr"),
                           .export = c("generate_data", 
                                       "forest_built",
                                       "check_trees_purity")) %dopar% {
                             
                             sim_data <- generate_data(seed = seeds[j], n = sample_size_now, k = k_val, d = d_val,
                                                       x_prob = xprob_vals, beta = beta_vals, gamma = gamma_vals)

                             result <- forest_built(sim_data, num_trees = num_trees_val, k = k_val)
                             
                             # Combine results as a row
                             c(result$heterogeneity_captured, 
                               result$success_ratio, 
                               result$partial_success_ratio, 
                               result$unsuccess_ratio)
                           }
        
        # Compute the average/median from the parallel results
        heterogeneity[i] <- median(results[, 1])
        success[i] <- median(results[, 2])
        partial_success[i] <- median(results[, 3])
        unsuccess[i] <- median(results[,4])
        
        end_time <- Sys.time()
        run_time <- round(difftime(end_time, start_time, units = "secs"), 4)
        
        log_msg <- paste0("Sample size ", sample_size_now, " : running time ", run_time, " sec")
        time_logs <- c(time_logs, log_msg)
        
        # Update Simulation Progress panel with logs
        progress_text(paste(time_logs, collapse = "\n"))
        
        # Update progress bar in the bottom-right corner
        incProgress(1 / n_param_grid, detail = log_msg)
      }
      
      stopCluster(cl)
      
      return_result <- data.frame(
        heterogeneity_captured = heterogeneity,
        success_ratio = success,
        partial_success_ratio = partial_success,
        unsuccess_ratio = unsuccess
      )
      
      return_result <- cbind(param_grid, return_result)
      
      # Final simulation status update
      progress_text(paste(c(time_logs, "Simulation completed!"), collapse = "\n"))
      
      # Render simulation results in a Plotly plot
      output$result_plot <- renderPlotly({
        
        plot_ly(data = return_result) %>%
          
          # Heterogeneity Captured
          add_trace(
            x = ~sample_size,
            y = ~heterogeneity_captured * 100,
            type = 'scatter',
            mode = 'lines+markers',
            name = 'Heterogeneity Captured',
            line = list(color = 'black', width = 2),
            marker = list(symbol = 'circle', size = 8, color = 'black'),
            hovertemplate = paste(
              'Sample Size: %{x}<br>',
              'Heterogeneity Captured: %{y:.2f}%<extra></extra>'
            )
          ) %>%
          
          # Success X
          add_trace(
            x = ~sample_size,
            y = ~success_ratio * 100,
            type = 'scatter',
            mode = 'lines+markers',
            name = 'Success X Tree',
            line = list(color = 'blue', width = 2),
            marker = list(symbol = 'circle', size = 8, color = 'blue'),
            hovertemplate = paste(
              'Sample Size: %{x}<br>',
              'Success X: %{y:.2f}%<extra></extra>'
            )
          ) %>%
          
          # Partial Success X
          add_trace(
            x = ~sample_size,
            y = ~partial_success_ratio * 100,
            type = 'scatter',
            mode = 'lines+markers',
            name = 'Partial Success X Tree',
            line = list(color = 'lightblue', width = 2),
            marker = list(symbol = 'circle', size = 8, color = 'lightblue'),
            hovertemplate = paste(
              'Sample Size: %{x}<br>',
              'Partial Success X: %{y:.2f}%<extra></extra>'
            )
          ) %>%
          
          # customize layout
          layout(
            title = '',
            xaxis = list(title = 'Sample Size', tickvals = return_result$sample_size,
                         ticktext = ifelse(return_result$sample_size == 0, "", return_result$sample_size)),
            yaxis = list(title = 'Percentage (%)', range = c(0, 105)),
            hovermode = 'closest'
          )
      })
    })
  })
  
  # show formula text
  output$formula_text <- renderUI({
    req(input$k, input$m)
    
    k <- input$k
    m <- input$m
    d <- k + m
    
    # beta term
    beta_terms <- paste0("\\beta_", 1:d, " x_{i", 1:d, "}")
    beta_expr <- paste(beta_terms, collapse = " + ")
    
    # gamma term
    gamma_expr <- if (k > 1) {
      paste0("\\gamma ", paste0("x_{i", 1:k, "}", collapse = ""))
      #paste0("\\gamma", paste0("x_{i", 1:k, "}"))
    } else if (k == 1) {
      paste0("\\gamma x_{i1}")
    } else {
      ""
    }
    
    # formula
    PY1 <- paste0(
      "P(Y_i(1)) = \\beta_0 + (", beta_expr, ")",
      if (k > 0) paste0(" + A_i (\\gamma_0 + ", gamma_expr, ")") else "",
      " + \\varepsilon"
    )
    
    Y1 <- paste0(
      "Y_i(1) \\mid X \\sim \\text{Bernoulli} \\left(\\text{logit}^{-1} \\left(\\beta_0 + ", 
      beta_expr, 
      if (k > 0) paste0(" + \\gamma_0 + ", gamma_expr) else "",
      " + \\varepsilon \\right) \\right)"
    ) 
    
    PY0 <- paste0("P(Y_i(0)) = \\beta_0 + ", beta_expr, " + \\varepsilon")
    
    Y0 <- paste0(
      "Y_i(0) \\mid X \\sim \\text{Bernoulli} \\left(\\text{logit}^{-1} \\left(\\beta_0 + ", 
      beta_expr, 
      " + \\varepsilon \\right) \\right)"
    )
    
    # Final tau
    tau_formula <- paste0(
      "\\tau(x_i) = \\gamma_0 + ", gamma_expr
    )
    
    # withMathJax(helpText(paste0("$$", latex_formula, "$$")))
    
    withMathJax(HTML(paste0(
      "<b>1. Data Structure</b><br>",
      "- Features: \\( X_i = (x_{i1}, \\ldots, x_{id}) \\in \\{0,1\\}^d \\)<br>",
      "- Response: \\( Y_i \\in \\{0,1\\} \\)<br>",
      "- First \\( k \\) covariates \\( x_1, \\ldots, x_k \\) are effect modifiers: covariates contribute to the treatment effect<br>", 
      "- \\( x_{k+1}, \\ldots, x_{d} \\) are other risk factors: covariates act as a distractor to the identification of variables that are relevant for the treatment effect<br>",
      
      "<b>2. Treatment Assignment</b><br>",
      "- Treatment indicator \\( A_i \\in \\{0,1\\} \\)<br>",
      "- Random assignment<br>",
      "- Error \\( \\varepsilon \\sim N(0, \\sigma^2) \\)<br>",
      
      "<b>3. Outcome Generation</b><br>",
      "<u>Treated units \\( (A_i = 1) \\):</u><br>",
      "$$", PY1, "$$",
      "$$", Y1, "$$<br>",
      
      "<u>Control units \\( (A_i = 0) \\):</u><br>",
      "$$", PY0, "$$",
      "$$", Y0, "$$<br>",
      
      "<b>4. Individual Causal Effect</b><br>",
      "$$\\tau(x_i) = E[Y_i(1) - Y_i(0) \\mid X_i = x]$$<br>",
      
      "<b>5. Ground-truth for \\( \\tau(x_i) \\)</b><br>",
      "$$", tau_formula, "$$<br>",
      
      "<hr style='border-top: 2px solid #bbb; margin-top: 20px; margin-bottom: 20px;'>",
      
      "<b>Example Simulation Data</b>"
    )))
    
  })
  
  
  # example data
  output$example_data <- renderDT({
    req(input$k, input$m)
    
    d <- input$k + input$m
    k <- input$k
    xprob_vals <- sapply(1:d, function(i) {
      val <- input[[paste0("xprob_", i)]]
      if (is.null(val)) 0.5 else as.numeric(val)
    })
    
    sim_data <- generate_data(
      n = 100,
      k = k,
      d = d,
      x_prob = xprob_vals,
      beta = c(input$beta0, sapply(1:d, function(i) input[[paste0("beta_", i)]])),
      # gamma = if (k == 1 || input$model_type == "pheno") {
      #   c(input$gamma0, input$gamma)
      # } else {
      #   c(input$gamma0, sapply(1:k, function(i) input[[paste0("gamma_", i)]]))
      # }
      gamma = c(input$gamma0, input$gamma)
    )
    
    sim_data$CATE <- round(sim_data$CATE, 3)
    sim_data <- sim_data[, !names(sim_data) %in% "ID"]
    
    if (k == 1 && "phenotype" %in% colnames(sim_data)) {
      sim_data <- sim_data[, !names(sim_data) %in% "phenotype"]
    }
    
    sim_data
  })
  
  # ======================== Panel 2 - Tuning Parameters ===========================
  
  # Reactive value to store tuning simulation results
  tuning_results <- reactiveVal(NULL)
  
  # Observe tuning param run button
  observeEvent(input$run_tuning_param, {
    # Save all input values locally to avoid referencing `input` inside parallel workers
    overall_sample_param_val <- input$overall_sample_param
    sample_fraction_range <- input$sample_fraction
    honesty_fraction_range <- input$honesty_fraction
    
    d_val <- input$m + input$k
    k_val <- input$k
    
    beta0 <- input$beta0
    beta_vals <- sapply(1:d_val, function(i) input[[paste0("beta_", i)]])
    beta_vals <- c(beta0, beta_vals)
    
    gamma0 <- input$gamma0
    gamma_val <- input$gamma
    gamma_vals <- c(gamma0, gamma_val)
    
    xprob_vals <- sapply(1:d_val, function(i) input[[paste0("xprob_", i)]])
    
    num_trees_val <- input$num_trees
    sim_no_val <- input$sim_no
    
    # Create the parameter grid to search
    param_list <- list(
      sample_fractions = seq(sample_fraction_range[1], sample_fraction_range[2], by = 0.1),
      honesty_fractions = seq(honesty_fraction_range[1], honesty_fraction_range[2], by = 0.1)
    )
    param_grid <- expand.grid(param_list)
    n_param_grid <- nrow(param_grid)
    
    heterogeneity <- numeric(n_param_grid)
    success <- numeric(n_param_grid)
    partial_success <- numeric(n_param_grid)
    unsuccess <- numeric(n_param_grid)
    
    # Initialize progress logs
    tuning_progress_logs <- reactiveVal("")
    tuning_progress_logs("Starting tuning simulation...\n")
    
    output$tuning_progress_status <- renderText({
      tuning_progress_logs()
    })
    
    withProgress(message = 'Running Tuning Simulation...', value = 0, {
      # Set up parallel backend
      num_cores <- detectCores() - 1
      cl <- makeCluster(num_cores)
      registerDoParallel(cl)
      
      # Loop through each parameter combination
      for (i in seq_len(n_param_grid)) {
        sample_fraction_now <- param_grid$sample_fractions[i]
        honesty_fraction_now <- param_grid$honesty_fractions[i]
        
        set.seed(1)
        seeds <- sample(1:(sim_no_val * 2), sim_no_val, replace = FALSE)
        
        start_time <- Sys.time()
        
        # Parallel simulation over seeds
        results <- foreach(j = 1:sim_no_val,
                           .combine = rbind,
                           .packages = c("grf", "dplyr"),
                           .export = c("generate_data", 
                                       "forest_built",
                                       "check_trees_purity")) %dopar% {
                                         
                                         sim_data <- generate_data(
                                           seed = seeds[j],
                                           n = overall_sample_param_val,
                                           k = k_val,
                                           d = d_val,
                                           x_prob = xprob_vals,
                                           beta = beta_vals,
                                           gamma = gamma_vals
                                         )
                                         
                                         result <- forest_built(
                                           data = sim_data,
                                           num_trees = num_trees_val,
                                           k = k_val,
                                           sample.fraction = sample_fraction_now,
                                           honesty.fraction = honesty_fraction_now
                                         )
                                         
                                         c(result$heterogeneity_captured, 
                                           result$success_ratio, 
                                           result$partial_success_ratio, 
                                           result$unsuccess_ratio)
                                       }
        
        # Store median results for each parameter set
        heterogeneity[i] <- median(results[, 1])
        success[i] <- median(results[, 2])
        partial_success[i] <- median(results[, 3])
        unsuccess[i] <- median(results[, 4])
        
        end_time <- Sys.time()
        run_time <- round(difftime(end_time, start_time, units = "secs"), 2)
        
        new_log <- paste0("Sample fraction=", sample_fraction_now, 
                          ", Honesty fraction=", honesty_fraction_now, 
                          " : running time ", run_time, " sec")
        
        new_log2 <- paste0(
          "Sample fraction=", sample_fraction_now, 
          ", Honesty fraction=", honesty_fraction_now, 
          " : running time ", run_time, " sec",
          ", Heterogeneity Captured=", round(heterogeneity[i],2),
          ", Success Ratio of Pure Effect Modifier Leaf=", round(success[i],2)
        )
        
        # Update simulation progress logs
        current_logs <- tuning_progress_logs()
        tuning_progress_logs(paste(current_logs, new_log2, sep = "\n"))
        
        # Update loading bar progress
        incProgress(1 / n_param_grid, detail = new_log)
      }
      
      stopCluster(cl)
      
      # Save final tuning simulation results
      tuning_result_final <- data.frame(
        sample_fractions = param_grid$sample_fractions,
        honesty_fractions = param_grid$honesty_fractions,
        heterogeneity_captured = heterogeneity,
        success_ratio = success,
        partial_success_ratio = partial_success,
        unsuccess_ratio = unsuccess
      )
      
      tuning_results(tuning_result_final)
      
      final_logs <- tuning_progress_logs()
      tuning_progress_logs(paste(final_logs, "\nTuning simulation completed!"))
    })
  })
  
  
  # Heatmap: Heterogeneity captured
  output$tuning_plot <- renderPlot({
    req(tuning_results())
    
    tuning_df <- tuning_results()
    
    ggplot(tuning_df, aes(x = sample_fractions, 
                          y = honesty_fractions, 
                          fill = heterogeneity_captured)) +
      geom_tile(color = "white") +  
      scale_fill_viridis_c(option = "viridis", direction = 1) +
      scale_x_continuous(breaks = unique(tuning_df$sample_fractions)) +
      scale_y_continuous(breaks = unique(tuning_df$honesty_fractions)) +
      labs(
        title = "Heatmap of Heterogeneity Captured",
        x = "Sample Fractions",
        y = "Honesty Fractions",
        fill = "Heterogeneity"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(vjust = 1, hjust = 1)  
      )
  })
  
  # Heatmap: Success ratio of pure effect modifier leaf
  output$tuning_plot2 <- renderPlot({
    req(tuning_results())
    
    tuning_df <- tuning_results()
    
    ggplot(tuning_df, aes(x = sample_fractions, 
                          y = honesty_fractions, 
                          fill = success_ratio)) +
      geom_tile(color = "white") +  
      scale_fill_viridis_c(option = "viridis", direction = 1) +
      scale_x_continuous(breaks = unique(tuning_df$sample_fractions)) +
      scale_y_continuous(breaks = unique(tuning_df$honesty_fractions)) +
      labs(
        title = "Heatmap of Success Ratio of Pure Effect Modifier Leaf",
        x = "Sample Fractions",
        y = "Honesty Fractions",
        fill = "Success Ratio"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(vjust = 1, hjust = 1)  
      )
  })

  
}

# ======================== RUN APP ===========================
shinyApp(ui, server)