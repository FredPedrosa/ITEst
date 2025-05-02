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
    labs(title = "Resíduos vs. Valores Ajustados",
         x = "Valores Ajustados",
         y = "Resíduos")
  
  # Q-Q Plot to check normality of residuals
  plot_qq <- ggplot(subject, aes(sample = Residuals)) +
    stat_qq() +
    stat_qq_line(color = "blue") +
    theme_classic() +
    labs(title = "Gráfico Q-Q dos Resíduos")
  
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
      geom_point(aes(y = Score, color = "Dados Brutos"), size = 2) +  # Points for observed data
      geom_line(aes(y = Predicted, color = "Modelo Combinado"), linewidth = 1) +  # Line for predicted values
      scale_color_manual(values = c("Dados Brutos" = "black", "Modelo Combinado" = "blue"),
                         name = "Legenda",
                         labels = c("Dados brutos", "HLM")) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = length(scores))) +  # Adjust x-axis breaks based on number of sessions
      labs(subtitle = "Estimado pelo Modelo HLM",
           x = "Sessões", y = "Pontuação") +  # Axis labels
      ylim(limits) +  # Y-axis limits
      theme_classic()  # Classic theme for the plot
    
    # Return the plot
    return(p)
  }
  
  # Create the prediction plot using the internal function
  plot_model <- plot_harmonic_linear(list(adjusted_data = subject))
  
  # Final model results
  cat("\n========== Resultados do Modelo Harmônico Linear ==========\n")
  cat("R²: ", round(R2_combined, 4), "\n")
  cat("Normalidade dos Resíduos (valor-p): ", round(shapiro_result$p.value, 4), "\n")
  
  # Determine trend based on the linear regression coefficient
  trend <- ifelse(coef(linear_model)[2] > 0, "Crescimento", "Declínio")
  
  # Output to indicate the individual's trend
  cat("Tendência do Indivíduo: ", trend, "\n\n")
  
  cat("Parâmetros Otimizados do Modelo Harmônico:\n")
  cat("Intercepto (a): ", round(best_a, 4), "\n")
  cat("Coeficiente do Seno (b): ", round(best_b, 4), "\n")
  cat("Coeficiente do Cosseno (c): ", round(best_c, 4), "\n")
  cat("Tendência Linear (f): ", round(best_f, 4), "\n")
  cat("Frequência (freq): ", round(best_freq, 4), "\n")
  cat("R² harmônico: ", round(R2_harmonic, 4), "\n\n")
  
  cat("Parâmetros do Modelo Linear:\n")
  cat("Intercepto Linear: ", round(coef(linear_model)[1], 4), "\n")
  cat("Coeficiente Linear (Inclinação): ", round(coef(linear_model)[2], 4), "\n")
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
  
  # Carregar pacotes necessários
  if (!requireNamespace("nls.multstart", quietly = TRUE)) {
    install.packages("nls.multstart")
  }
  library(nls.multstart)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }
  library(ggplot2)
  
  # Função do modelo logístico de 4 parâmetros
  M.4pl <- function(Session, lower.asymp, upper.asymp, inflec, hill) {
    lower.asymp + ((upper.asymp - lower.asymp) / (1 + (Session / inflec)^-hill))
  }
  
  # Definir pontos de inflexão
  if (is.null(inflec)) inflec <- c(2, length(scores) - 1)
  
  suj <- data.frame(Score = scores, Session = 1:length(scores))
  
  # Ajustar tolerância se todos os valores forem iguais
  if (all(scores == scores[1])) {
    tol <- c(.2, .1)
  }
  
  list_m <- vector("list", 100)
  i <- 1
  j <- 1
  k <- 101
  
  # Loop de otimização
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
    
    # Mensagem de execução
    cat("\014")  # Limpa o console
    cat("\n =========== Rodando o modelo 4PL para as pontuações ==============\n")
    cat("Repetições realizadas: ", j, "\n")
    cat(sprintf("\r Repetições válidas: %d de %d\n", i, k - 1))
    
    j <- j + 1
    if (j > 100 & j < 200) k <- 50
    if (j >= 200) k <- 30
    if (j == max.iter) i <- k
  }
  
  # Filtrar modelos válidos
  valid_models <- list_m[!sapply(list_m, is.null)]
  
  # Resultados do modelo
  if (length(valid_models) > 0) {
    results <- data.frame(do.call(rbind, lapply(valid_models, coef)), 
                          resid = unlist(lapply(valid_models, function(x) sum(resid(x)^2))))
    results$R2 <- 1 - results$resid / sum((scores - mean(scores))^2)
    
    bm_4pl <- which.min(results$resid)
    
    # Adicionar os dados ajustados ao dataframe
    suj$Predicted <- predict(valid_models[[bm_4pl]], newdata = suj)
    
    # Criar o gráfico
    p <- ggplot(data = suj, aes(x = Session)) +
      geom_point(aes(y = Score, color = "Dados brutos")) +  # Dados brutos
      geom_line(aes(y = Predicted, color = "Modelo 4-PL"), linewidth = 1, linetype = "dashed") +  # Modelo 4-PL
      scale_x_continuous(breaks = scales::pretty_breaks(n = length(scores))) +
      scale_color_manual(name = "Legenda", values = c("Dados brutos" = "black", "Modelo 4-PL" = "blue")) +
      labs(subtitle = "Estimado pelo Modelo 4-PL", x = "Sessão", y = "Pontuações") + ylim(limit) +
      theme_classic()
    
    # Exibir o gráfico
    print(p)
    
    # Imprimir resultados no console
    cat("\n ================== Resultados do modelo 4PL ==================\n")
    cat("R² do melhor modelo: ", round(results$R2[bm_4pl], 4), "\n")
    cat("Parâmetros do melhor modelo:\n")
    print(results[bm_4pl, 1:6])
    
    # Retornar os resultados e o gráfico
    return(list(
      results = results,                    # Resultados do modelo
      best_model = results[bm_4pl,],        # Parâmetros do melhor modelo
      plot = p,                             # Gráfico final
      best_model_object = valid_models[[bm_4pl]]  # Objeto do melhor modelo
    ))
  } else {
    stop("O modelo não convergiu.")
  }
}



interpret_4pl <- function(output_4pl) {
  
  # Arredondar os valores dos parâmetros do modelo
  lower_asymp <- round(output_4pl$best_model["lower.asymp"], 2)
  upper_asymp <- round(output_4pl$best_model["upper.asymp"], 2)
  inflec <- round(output_4pl$best_model["inflec"], 2)
  hill <- round(output_4pl$best_model["hill"], 2)
  R2 <- round(output_4pl$best_model["R2"], 2) * 100
  
  # Determinar a direção do desenvolvimento com base no parâmetro 'hill'
  development_direction <- if (hill > 0) {
    "uma trajetória de crescimento positivo."
  } else {
    "um declínio na trajetória."
  }
  
  # Construção da interpretação final
  interpretation <- paste0(
    "O modelo logístico de 4 parâmetros (4PL) indica ", development_direction, "\n\n",
    "O modelo 4PL (Gomes & Farias apud Araujo & Farias, 2024) explica aproximadamente ", R2, "% da 
    variância nos dados. O valor da assíntota inferior é ", lower_asymp, ", representando a 
    pontuação mais baixa estimada pelo modelo, e o valor da assíntota superior é ", upper_asymp, ", representando 
    a pontuação mais alta estimada pelo modelo. O ponto de inflexão ocorre na sessão ", inflec, ", que 
    marca o momento de maior mudança nas pontuações. A inclinação da curva, ", hill, ", 
    reflete a taxa de mudança ao redor desse ponto.", "\n\n",
    "Referência:", "\n",
    "Araújo, J. de, & Farias, H. B. (2024). Avaliando a trajetória do processo psicológico 
    do indivíduo por meio de modelos. I Congresso Brasileiro de Psicometria e Análise de Dados, 
    Porto Alegre. https://www.researchgate.net/publication/381741254_Avaliando_a_trajetoria_
    do_processo_psicologico_do_individuo_por_meio_de_modelos", "\n\n",
    "Como citar:", "\n",
    "Pedrosa, F.G. (2024). _Estimação da trajetória individual: modelo de 4 parâmetros logísticos e
    modelo harmônico linear_. [Software]. https://fredpedrosa.shinyapps.io/app_pt/"
  )
  
  # Exibir a interpretação com quebra de linha
  cat(interpretation, "\n")
}

interpret_HLM <- function(output_HLM) {
  # Arredondar os parâmetros do modelo harmônico
  intercept <- round(output_HLM$harmonic_parameters$intercept, 2)
  sine_coef <- round(output_HLM$harmonic_parameters$sine_coef, 2)
  cosine_coef <- round(output_HLM$harmonic_parameters$cosine_coef, 2)
  linear_trend <- round(output_HLM$harmonic_parameters$linear_trend, 2)
  frequency <- round(output_HLM$harmonic_parameters$frequency, 2)
  harmonic_R2 <- round(output_HLM$harmonic_parameters$R2_har, 2) * 100
  
  # Arredondar os parâmetros do modelo linear
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
    cat("Os dados apresentam um padrão não-harmônico. Este modelo deve ser desconsiderado.\n")
    return()  # Sai da função sem fazer mais nada
  }
  
  # Verificação de harmonia
  harmony_check <- "indicando que os dados seguem um padrão harmônico, mas também possuem uma tendência linear."
  
  # Verificação da normalidade dos resíduos
  normality_check <- if (is.na(normality_pvalue)) {
    "o resultado do teste de normalidade não está disponível."
  } else if (normality_pvalue < 0.05) {
    "o modelo deve ser interpretado com cautela, pois pode haver problemas com suas suposições."
  } else {
    "o modelo não apresenta problemas com suas suposições."
  }
  
  # Determinar a direção da tendência geral
  trend_direction <- if (linear_trend > 0) {
    "a trajetória do indivíduo indica crescimento."
  } else {
    "a trajetória do indivíduo indica declínio."
  }
  
  # Construção da interpretação final
  interpretation <- paste0(
    "Os dados apresentam um padrão harmônico e ", trend_direction, "\n\n", 
    "O Modelo Harmônico Linear (Pedrosa, 2024) explicou ", combined_R2, "% da variação nos dados, e o teste de 
    normalidade dos resíduos informou que ", normality_check, " 
    O coeficiente de seno é ", sine_coef, ", e o coeficiente de cosseno é ", cosine_coef, ", apontando um padrão cíclico nos dados. 
    A tendência linear é ", linear_trend, ", e a frequência é ", frequency, ", 
    ", harmony_check, "\n\n",
    "Referência:", "\n",
    "Pedrosa, F.G. (2024). Harmonic Linear Model [R]. https://github.com/FredPedrosa/HarmonicLinearModel", "\n\n",
    "Como citar:", "\n",
    "Pedrosa, F.G. (2024). _Estimação da trajetória individual: modelo de 4 parâmetros logísticos e
    modelo harmônico linear_. [Software]. https://fredpedrosa.shinyapps.io/app_pt/"
  )
  
  # Exibir a interpretação com espaçamento adequado
  cat(interpretation, "\n")
}



# Define UI for application
ui <- fluidPage(
  titlePanel("Estimação da trajetória individual: modelo de 4 parâmetros logísticos 
  e modelo harmônico linear"),
  
  # Adicionando o subtítulo
  tags$h4("Nota: Como o modelo 4-PL gera múltiplos ajustes durante a análise para selecionar o melhor modelo, o processo pode levar alguns minutos. 
          Para assegurar uma estimativa precisa da trajetória individual, recomenda-se o uso de instrumentos de medida que apresentem evidências robustas de validade no contexto intrapessoal."),
  
  sidebarLayout(
    sidebarPanel(
      textInput("scores", "Insira os escores (separados por vírgula):", value = "35, 36, 35, 31, 32, 35, 41, 47, 43"),
      numericInput("min_limit", "Limite Mínimo dos Escores:", value = 0),
      numericInput("max_limit", "Limite Máximo dos Escores:", value = 68),
      
      # Adicionando seleção de modelo
      checkboxGroupInput("model_choice", "Escolha o(s) modelo(s) a ser(em) rodado(s):", 
                         choices = list("Modelo Harmônico Linear (HLM)" = "HLM", 
                                        "Modelo de 4 Parâmetros Logísticos (4PL)" = "4PL"),
                         selected = c("HLM", "4PL")),
      
      actionButton("run", "Executar Análise")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Modelo Harmônico Linear (HLM)",
                 plotOutput("plot_hlm"),
                 div(style = "white-space: pre-wrap;", verbatimTextOutput("interpret_hlm"))
        ),
        tabPanel("Modelo de 4 Parâmetros Logísticos (4PL)",
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
    # Lê os escores inseridos como uma string e os converte em um vetor numérico
    scores <- as.numeric(unlist(strsplit(input$scores, ",")))
    
    # Verifica se todos os valores são numéricos
    if (any(is.na(scores))) {
      output$interpret_hlm <- renderPrint({
        "Erro: Alguns valores inseridos não são números válidos. Por favor, insira escores numéricos separados por vírgulas."
      })
      output$interpret_4pl <- renderPrint({
        "Erro: Alguns valores inseridos não são números válidos. Por favor, insira escores numéricos separados por vírgulas."
      })
      return()
    }
    
    min_limit <- input$min_limit
    max_limit <- input$max_limit
    
    # Executa o HLM se selecionado
    if ("HLM" %in% input$model_choice) {
      result_HLM <- EstIT_HLM(scores)
      
      output$plot_hlm <- renderPlot({
        result_HLM$plot + ylim(min_limit, max_limit)
      })
      
      output$interpret_hlm <- renderPrint({
        cat("========== Resultados do Modelo Harmônico Linear ==========\n")
        cat("R²: ", round(result_HLM$R2,2), "\n")
        cat("Normalidade dos resíduos (p-value): ", round(result_HLM$shapiro_p_value,2), "\n\n")
        cat("Parâmetros Otimizados do HML:\n")
        cat("Intercepto: ", round(result_HLM$harmonic_parameters$intercept,2), "\n")
        cat("Coeficiente senoidal: ", round(result_HLM$harmonic_parameters$sine_coef,2), "\n")
        cat("Coeficiente cossenoidal: ", round(result_HLM$harmonic_parameters$cosine_coef,2), "\n")
        cat("Tendência Linear: ", round(result_HLM$harmonic_parameters$linear_trend,2), "\n")
        cat("Frequência: ", round(result_HLM$harmonic_parameters$frequency,2), "\n")
        cat("R² harmônico: ", round(result_HLM$harmonic_parameters$R2_har,2), "\n\n")
        cat("Parâmetros da Modelo Linear:\n")
        cat("Intercepto Linear: ", round(result_HLM$linear_parameters$intercept,2), "\n")
        cat("Inclinação (Slope): ", round(result_HLM$linear_parameters$slope_coef,2), "\n")
        cat("R² linear: ", round(result_HLM$linear_parameters$R2_lin,2), "\n\n")
        
        interpret_HLM(result_HLM)
      })
    }
    
    # Executa o 4PL se selecionado
    if ("4PL" %in% input$model_choice) {
      result_4pl <- EstIT_4pl(scores)
      
      output$plot_4pl <- renderPlot({
        result_4pl$plot + ylim(min_limit, max_limit)
      })
      
      output$interpret_4pl <- renderPrint({
        cat("================== 4PL model results ==================\n")
        cat("R² do melhor modelo: ", round(result_4pl$best_model$R2,2), "\n")
        cat("Parâmetros do melhor modelo:\n")
        cat("Assíntota inferior: ", round(result_4pl$best_model$lower.asymp,2), "\n")
        cat("Assíntota superior: ", round(result_4pl$best_model$upper.asymp,2), "\n")
        cat("Ponto de inflexão: ", round(result_4pl$best_model$inflec,2), "\n")
        cat("Inclinação da curva: ", round(result_4pl$best_model$hill,2), "\n")
        cat("Resíduos: ", round(result_4pl$best_model$resid,2), "\n\n")
        
        interpret_4pl(result_4pl)
      })
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)