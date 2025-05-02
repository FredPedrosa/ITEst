library(shiny)
library(ggplot2)
library(nls.multstart)

EstIT_HLM <- function(scores) {
  # Load necessary packages
  library(ggplot2)
  
  # Error function for harmonic model optimization
  harmonic_model_error <- function(params, scores) {
    freq <- params[1]  # Frequency
    f <- params[2]     # Linear Trend
    a <- params[3]     # Intercept
    b <- params[4]     # Sine Coefficient
    c <- params[5]     # Cosine Coefficient
    
    subject <- data.frame(Score = scores, Session = seq_along(scores))
    
    # Harmonic model with optimized frequency and linear trend
    subject$Predicted <- a + b * sin(freq * subject$Session) + c * cos(freq * subject$Session) + f * subject$Session
    
    # Error to be minimized (squared residuals)
    sum((subject$Score - subject$Predicted)^2)
  }
  
  # Initial parameters for the harmonic model
  start_params <- c(freq = 1, f = 0, a = mean(scores), b = 1, c = 1)
  
  # Optimization of harmonic model parameters
  optim_result <- optim(par = start_params, fn = harmonic_model_error, scores = scores)
  best_freq <- optim_result$par[1]
  best_f <- optim_result$par[2]
  best_a <- optim_result$par[3]
  best_b <- optim_result$par[4]
  best_c <- optim_result$par[5]
  
  # Final fit of the harmonic model with optimized parameters
  subject <- data.frame(Score = scores, Session = seq_along(scores))
  subject$Predicted_Harmonic <- best_a + best_b * sin(best_freq * subject$Session) + 
    best_c * cos(best_freq * subject$Session) + best_f * subject$Session
  
  # Fit of a simple linear regression model
  linear_model <- lm(Score ~ Session, data = subject)
  subject$Predicted_Linear <- predict(linear_model, newdata = subject)
  
  # Calculate R² for each model
  R2_harmonic <- 1 - sum((subject$Score - subject$Predicted_Harmonic)^2) / sum((subject$Score - mean(subject$Score))^2)
  R2_linear <- summary(linear_model)$r.squared
  
  # Weighted combination of harmonic and linear models by R²
  subject$Predicted_Combined <- (R2_harmonic * subject$Predicted_Harmonic + R2_linear * subject$Predicted_Linear) / 
    (R2_harmonic + R2_linear)
  
  # Calculate residuals
  subject$Residuals <- subject$Score - subject$Predicted_Combined
  
  # Calculate R² for the combined model
  R2_combined <- 1 - sum((subject$Score - subject$Predicted_Combined)^2) / sum((subject$Score - mean(subject$Score))^2)
  
  # Residuals vs. Fitted Values Plot
  plot_residuals <- ggplot(subject, aes(x = Predicted_Combined, y = Residuals)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
    theme_classic() +
    labs(title = "Residuals vs. Fitted Values",
         x = "Fitted Values",
         y = "Residuals")
  
  # Q-Q Plot to check normality of residuals
  plot_qq <- ggplot(subject, aes(sample = Residuals)) +
    stat_qq() +
    stat_qq_line(color = "blue") +
    theme_classic() +
    labs(title = "Residuals Q-Q plot")
  
  # Shapiro-Wilk test for normality of residuals
  shapiro_result <- shapiro.test(subject$Residuals)
  
  # Function to create the harmonic linear prediction plot
  plot_harmonic_linear <- function(model_result, limits = c(0, 68)) {
    # Load the ggplot2 package
    library(ggplot2)
    
    # Extract the adjusted data from the model result
    subject <- model_result$adjusted_data
    
    # Create a dataframe with session values and observed and predicted scores
    df <- data.frame(
      Session = subject$Session,              # Sessions
      Score = subject$Score,                  # Observed scores
      Predicted = subject$Predicted_Combined  # Predicted scores (combined by the model)
    )
    
    # Build the plot with observed and predicted scores
    p <- ggplot(df, aes(x = Session)) +
      geom_point(aes(y = Score, color = "Raw Data"), size = 2) +  # Points for observed data
      geom_line(aes(y = Predicted, color = "Combined Model"), linewidth = 1) +  # Line for predicted values
      scale_color_manual(values = c("Raw Data" = "black", "Combined Model" = "blue"),
                         name = "Legend",
                         labels = c("HL Model", "Raw Data" )) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = length(scores))) +  # Adjust x-axis breaks based on number of sessions
      labs(subtitle = "Estimated by HLM Model",
           x = "Sessions", y = "Score") +  # Axis labels
      ylim(limits) +  # Y-axis limits
      theme_classic()  # Classic theme for the plot
    
    # Return the plot
    return(p)
  }
  
  # Create the prediction plot using the internal function
  plot_model <- plot_harmonic_linear(list(adjusted_data = subject))
  
  # Final model results
  cat("\n========== Harmonic Linear Model Results ==========\n")
  cat("R²: ", round(R2_combined, 4), "\n")
  cat("Residual Normality (p-value): ", round(shapiro_result$p.value, 4), "\n")
  
  # Determine trend based on the linear regression coefficient
  trend <- ifelse(coef(linear_model)[2] > 0, "Growth", "Decline")
  
  # Output to indicate the individual's trend
  cat("Individual's Trend: ", trend, "\n\n")
  
  cat("Optimized Parameters of Harmonic Model:\n")
  cat("Intercept (a): ", round(best_a, 4), "\n")
  cat("Sine Coefficient (b): ", round(best_b, 4), "\n")
  cat("Cosine Coefficient (c): ", round(best_c, 4), "\n")
  cat("Linear Trend (f): ", round(best_f, 4), "\n")
  cat("Frequency (freq): ", round(best_freq, 4), "\n")
  cat("R² harmonic: ", round(R2_harmonic, 4), "\n\n")
  
  cat("Parameters of Linear Model:\n")
  cat("Linear Intercept: ", round(coef(linear_model)[1], 4), "\n")
  cat("Linear Coefficient (Slope): ", round(coef(linear_model)[2], 4), "\n")
  cat("R² linear: ", round(R2_linear, 4), "\n\n")
  
  # Return results 
  return(list(
    R2 = R2_combined,
    shapiro_p_value = shapiro_result$p.value,  # Store Shapiro-Wilk p-value in the output
    adjusted_data = subject,
    plot = plot_model,
    harmonic_parameters = list(
      intercept = best_a,
      sine_coef = best_b,
      cosine_coef = best_c,
      linear_trend = best_f,
      frequency = best_freq,
      R2_har = R2_harmonic
    ),
    linear_parameters = list(
      intercept = coef(linear_model)[1],
      slope_coef = coef(linear_model)[2],
      R2_lin = R2_linear
    )
  ))
}




EstIT_4pl <- function(scores, tol = c(.3, .3), hill = c(-100000, 100000), inflec = NULL, max.iter = 5000, limit = c(0, 68)) {
  
  # Load necessary packages
  if (!requireNamespace("nls.multstart", quietly = TRUE)) {
    install.packages("nls.multstart")
  }
  library(nls.multstart)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }
  library(ggplot2)
  
  # Logistic 4-parameter model function
  M.4pl <- function(Session, lower.asymp, upper.asymp, inflec, hill) {
    lower.asymp + ((upper.asymp - lower.asymp) / (1 + (Session / inflec)^-hill))
  }
  
  # Set inflection points
  if (is.null(inflec)) inflec <- c(2, length(scores) - 1)
  
  suj <- data.frame(Score = scores, Session = 1:length(scores))
  
  # Adjust tolerance if all values are equal
  if (all(scores == scores[1])) {
    tol <- c(.2, .1)
  }
  
  list_m <- vector("list", 100)
  i <- 1
  j <- 1
  k <- 101
  
  # Optimization loop
  while (i < k) {
    modelo_p <- nls_multstart(
      Score ~ M.4pl(Session, lower.asymp, upper.asymp, inflec, hill),
      data = suj,
      lower = c(lower.asymp = min(scores) - tol[1], upper.asymp = max(scores) - tol[2], inflec = inflec[1], hill = hill[1]),
      upper = c(lower.asymp = min(scores) + tol[2], upper.asymp = max(scores) + tol[1], inflec = inflec[2], hill = hill[2]),
      start_lower = c(lower.asymp = min(scores) - tol[1], upper.asymp = max(scores) - tol[2], inflec = inflec[1], hill = hill[1]),
      start_upper = c(lower.asymp = min(scores) + tol[2], upper.asymp = max(scores) + tol[1], inflec = inflec[2], hill = hill[2]),
      iter = 5000, supp_errors = "Y"
    )
    
    if (!is.null(modelo_p)) {
      list_m[[i]] <- modelo_p
      i <- i + 1
    }
    
    # Execution message
    cat("\014")  # Clears the console
    cat("\n =========== Running 4PL model for scores ==============\n")
    cat("Repetitions performed: ", j, "\n")
    cat(sprintf("\r Valid repetitions: %d out of %d\n", i, k - 1))
    
    j <- j + 1
    if (j > 100 & j < 200) k <- 50
    if (j >= 200) k <- 30
    if (j == max.iter) i <- k
  }
  
  # Filter valid models
  valid_models <- list_m[!sapply(list_m, is.null)]
  
  # Model results
  if (length(valid_models) > 0) {
    results <- data.frame(do.call(rbind, lapply(valid_models, coef)), 
                          resid = unlist(lapply(valid_models, function(x) sum(resid(x)^2))))
    results$R2 <- 1 - results$resid / sum((scores - mean(scores))^2)
    
    bm_4pl <- which.min(results$resid)
    
    # Add the fitted data to the dataframe
    suj$Predicted <- predict(valid_models[[bm_4pl]], newdata = suj)
    
    # Create the plot
    p <- ggplot(data = suj, aes(x = Session)) +
      geom_point(aes(y = Score, color = "Raw data")) +  # Raw data
      geom_line(aes(y = Predicted, color = "4-PL model"), linewidth = 1, linetype = "dashed") +  # 4-PL model
      scale_x_continuous(breaks = scales::pretty_breaks(n = length(scores))) +
      scale_color_manual(name = "Legend", values = c("Raw data" = "black", "4-PL model" = "blue")) +
      labs(subtitle = "Estimated by 4-PL Model", x = "Session", y = "Scores") + ylim(limit) +
      theme_classic()
    
    # Display the plot
    print(p)
    
    # Print results in the console
    cat("\n ================== 4PL model results ==================\n")
    cat("Best model's R²: ", round(results$R2[bm_4pl], 4), "\n")
    cat("Best model's parameters:\n")
    print(results[bm_4pl, 1:6])
    
    # Return the results and the plot
    return(list(
      results = results,                    # Model results
      best_model = results[bm_4pl,],        # Best model parameters
      plot = p,                             # Final plot
      best_model_object = valid_models[[bm_4pl]]  # Best model object
    ))
  } else {
    stop("The model didn't converge.")
  }
}

interpret_4pl <- function(output_4pl) {
  
  # Arredondar os valores
  lower_asymp <- round(output_4pl$best_model["lower.asymp"], 2)
  upper_asymp <- round(output_4pl$best_model["upper.asymp"], 2)
  inflec <- round(output_4pl$best_model["inflec"], 2)
  hill <- round(output_4pl$best_model["hill"], 2)
  R2 <- round(output_4pl$best_model["R2"], 2) * 100
  
  # Determinar a direção do desenvolvimento com base em 'hill'
  development_direction <- if (hill > 0) {
    "a positive growth trajectory."
  } else {
    "a decline in the trajectory."
  }
  
  # Construção da interpretação
  interpretation <- paste0(
    "The 4-parameter logistic model (4-pl) indicates ", development_direction, "\n\n",
    "The 4-pl model (Gomes & Farias apud Araujo & Farias, 2024) explains approximately ", R2, "% of 
    the variance in the data. The lower asymptote value is ", lower_asymp, ", representing the 
    lowest score estimated by the model, and the upper asymptote value is ", upper_asymp, ", representing 
    the highest score estimated by the model. The inflection point occurs at session ",inflec, ", which 
    marks the moment where the steepest change in scores happens. The slope of the curve, ", hill, ", 
    reflects the rate of change around this point.", "\n\n",
    "Reference:", "\n",
    "Araújo, J. de, & Farias, H. B. (2024). Avaliando a trajetória do processo psicológico 
    do indivíduo por meio de modelos. I Congrsso Brasileiro de Psicometria e Análise de Dados, 
    Porto Alegre. https://www.researchgate.net/publication/381741254_Avaliando_a_trajetoria_
    do_processo_psicologico_do_individuo_por_meio_de_modelos", "\n\n",
    "How to cite:", "\n",
    "Pedrosa, F.G. (2024). _Estimation of Individual Trajectory: 4-Parameter Logistic Model 
    and Harmonic Linear Modelr_.[Software]. https://fredpedrosa.shinyapps.io/estit/"
  )
  
  # Exibir a interpretação com quebra de linha
  cat(interpretation, "\n")
}


interpret_HLM <- function(output_HLM) {
  # Acessar os parâmetros corretos da saída do modelo
  intercept <- round(output_HLM$harmonic_parameters$intercept, 2)
  sine_coef <- round(output_HLM$harmonic_parameters$sine_coef, 2)
  cosine_coef <- round(output_HLM$harmonic_parameters$cosine_coef, 2)
  linear_trend <- round(output_HLM$harmonic_parameters$linear_trend, 2)
  frequency <- round(output_HLM$harmonic_parameters$frequency, 2)
  harmonic_R2 <- round(output_HLM$harmonic_parameters$R2_har, 2) * 100
  
  linear_intercept <- round(output_HLM$linear_parameters$intercept, 2)
  slope <- round(output_HLM$linear_parameters$slope_coef, 2)
  linear_R2 <- round(output_HLM$linear_parameters$R2_lin, 2) * 100
  
  combined_R2 <- round(output_HLM$R2, 2) * 100
  
  # Verificar se o valor-p do teste de normalidade dos resíduos está disponível
  normality_pvalue <- if (!is.null(output_HLM$shapiro_p_value)) {
    round(output_HLM$shapiro_p_value, 4)
  } else {
    NA  # Valor padrão caso o teste não esteja disponível
  }
  
  # Verificar se a frequência é negativa para indicar padrão não-harmônico
  if (frequency < 0) {
    cat("The data have non-harmonic pattern. This model should be disregarded.\n")
    return()  # Sai da função sem fazer mais nada
  }
  
  # Verificação de harmonia
  harmony_check <- "indicating the data follow a harmonic pattern, but with a linear trend as well."
  
  # Verificação da normalidade dos resíduos
  normality_check <- if (is.na(normality_pvalue)) {
    "normality test result is unavailable."
  } else if (normality_pvalue < 0.05) {
    "the model should be interpreted with caution as it may have issues with its assumptions."
  } else {
    "the model does not present problems with its assumptions."
  }
  
  # Determinar a direção da tendência geral
  trend_direction <- if (linear_trend > 0) {
    "the individual's trajectory indicates growth."
  } else {
    "the individual's trajectory indicates decline."
  }
  
  # Construção da interpretação final
  interpretation <- paste0(
    "The data have a harmonic pattern and ", trend_direction, "\n\n", 
    "The Harmonic Linear Model (Pedrosa, 2024) explained ", combined_R2, "% of the variation in the data, and the residual 
    normality test indicated that ", normality_check, "The sine 
    coefficient is ", sine_coef, ", and the cosine coefficient is ", cosine_coef, ", indicating a cyclical pattern. 
    The linear trend is ", linear_trend, ", and the frequency is ", frequency, ", 
    ", harmony_check, "\n\n",
    "Reference:", "\n",
    "Pedrosa, F.G.(2024). Harmonic Linear Model [R].
    https://github.com/FredPedrosa/HarmonicLinearModel", "\n\n",
    "How to cite:", "\n",
    "Pedrosa, F.G. (2024). _Estimation of Individual Trajectory: 4-Parameter Logistic Model 
    and Harmonic Linear Modelr_.[Software]. https://fredpedrosa.shinyapps.io/estit/"
  )
  
  # Exibir a interpretação com espaçamento adequado
  cat(interpretation, "\n")
}




# Define UI for application
ui <- fluidPage(
  titlePanel("Estimation of Individual Trajectory: 4-Parameter Logistic Model 
  and Harmonic Linear Model"),
  
  # Adding the subtitle
  tags$h4("Note: As the 4-PL model generates multiple fits during analysis to select the best model, the process may take a few minutes. 
          To ensure an accurate estimate of the individual trajectory, it is recommended to use measurement instruments that provide strong evidence of validity in an intrapersonal context."),
  
  sidebarLayout(
    sidebarPanel(
      textInput("scores", "Enter the scores (comma-separated):", value = "35, 36, 35, 31, 32, 35, 41, 47, 43"),
      numericInput("min_limit", "Minimum Score Limit:", value = 0),
      numericInput("max_limit", "Maximum Score Limit:", value = 68),
      
      # Adding model selection
      checkboxGroupInput("model_choice", "Choose the model(s) to run:", 
                         choices = list("Harmonic Linear Model (HLM)" = "HLM", 
                                        "4-Parameter Logistic Model (4PL)" = "4PL"),
                         selected = c("HLM", "4PL")),
      
      actionButton("run", "Run Analysis")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Harmonic Linear Model (HLM)",
                 plotOutput("plot_hlm"),
                 div(style = "white-space: pre-wrap;", verbatimTextOutput("interpret_hlm"))
        ),
        tabPanel("4-Parameter Logistic Model (4PL)",
                 plotOutput("plot_4pl"),
                 div(style = "white-space: pre-wrap;", verbatimTextOutput("interpret_4pl"))
        )
      )
    )
  )
)

# Define server logic required to draw the plots
server <- function(input, output) {
  
  observeEvent(input$run, {
    # Read the scores entered as a string and convert them into a numeric vector
    scores <- as.numeric(unlist(strsplit(input$scores, ",")))
    
    # Check if all values are numeric
    if (any(is.na(scores))) {
      output$interpret_hlm <- renderPrint({
        "Error: Some entered values are not valid numbers. Please input numeric scores separated by commas."
      })
      output$interpret_4pl <- renderPrint({
        "Error: Some entered values are not valid numbers. Please input numeric scores separated by commas."
      })
      return()
    }
    
    min_limit <- input$min_limit
    max_limit <- input$max_limit
    
    # Run HLM if selected
    if ("HLM" %in% input$model_choice) {
      result_HLM <- EstIT_HLM(scores)
      
      output$plot_hlm <- renderPlot({
        result_HLM$plot + ylim(min_limit, max_limit)
      })
      
      output$interpret_hlm <- renderPrint({
        cat("========== Harmonic Linear Model Results ==========\n")
        cat("R²: ", round(result_HLM$R2,2), "\n")
        cat("Residual normality (p-value): ", round(result_HLM$shapiro_p_value,2), "\n\n")
        cat("Optimized HML Parameters:\n")
        cat("Intercept: ", round(result_HLM$harmonic_parameters$intercept,2), "\n")
        cat("Sine coefficient: ", round(result_HLM$harmonic_parameters$sine_coef,2), "\n")
        cat("Cosine coefficient: ", round(result_HLM$harmonic_parameters$cosine_coef,2), "\n")
        cat("Linear trend: ", round(result_HLM$harmonic_parameters$linear_trend,2), "\n")
        cat("Frequency: ", round(result_HLM$harmonic_parameters$frequency,2), "\n")
        cat("Harmonic R²: ", round(result_HLM$harmonic_parameters$R2_har,2), "\n\n")
        cat("Linear Model Parameters:\n")
        cat("Linear Intercept: ", round(result_HLM$linear_parameters$intercept,2), "\n")
        cat("Slope coefficient: ", round(result_HLM$linear_parameters$slope_coef,2), "\n")
        cat("Linear R²: ", round(result_HLM$linear_parameters$R2_lin,2), "\n\n")
        
        interpret_HLM(result_HLM)
      })
    }
    
    # Run 4PL if selected
    if ("4PL" %in% input$model_choice) {
      result_4pl <- EstIT_4pl(scores)
      
      output$plot_4pl <- renderPlot({
        result_4pl$plot + ylim(min_limit, max_limit)
      })
      
      output$interpret_4pl <- renderPrint({
        cat("================== 4PL Model Results ==================\n")
        cat("R² of the best model: ", round(result_4pl$best_model$R2,2), "\n")
        cat("Best model parameters:\n")
        cat("Lower asymptote: ", round(result_4pl$best_model$lower.asymp,2), "\n")
        cat("Upper asymptote: ", round(result_4pl$best_model$upper.asymp,2), "\n")
        cat("Inflection point: ", round(result_4pl$best_model$inflec,2), "\n")
        cat("Curve slope: ", round(result_4pl$best_model$hill,2), "\n")
        cat("Residuals: ", round(result_4pl$best_model$resid,2), "\n\n")
        
        interpret_4pl(result_4pl)
      })
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

