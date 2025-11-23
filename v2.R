library(shiny)
library(shinyvalidate)
library(plotly)

ui <- fluidPage(
  titlePanel("PRS Simulation Tool"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("m", "Number of SNPs:", value = 100, min = 1, step = 1),
      numericInput("n1", "Training sample size (n1):", value = 1000, min = 1, step = 1),
      numericInput("n2", "Testing sample size (n2):", value = 500, min = 1, step = 1),
      
      numericInput("w_min", "Overlap proportion min (0–1):", value = 0.1, min = 0, max = 1),
      numericInput("w_max", "Overlap proportion max (0–1):", value = 0.5, min = 0, max = 1),
      numericInput("w_step", "Overlap step:", value = 0.1, min = 0.0001),
      
      numericInput("nsim", "Number of simulations per Overlap proportion:", value = 10, min = 1, step = 1),
      
      numericInput("seed", "Random seed:", value = 123),
      
      hr(),
      # ---- 下载按钮分两行 ----
      div(
        style = "margin-bottom:15px;",
        downloadButton("download_plot_png", "Download Plot (PNG)", class = "btn-success",
                       style="margin-right:10px;"),
        downloadButton("download_plot_pdf", "Download Plot (PDF)", class = "btn-warning")
      ),
      div(
        style = "margin-bottom:20px;",
        downloadButton("download_data", "Download Data (CSV)", class = "btn-info")
      ),
      div(
        style = "margin-bottom:20px;",
        actionButton("reset", "Reset to Default", class = "btn-danger")
      ),
      
      
      hr(),
      actionButton("run", "Run Simulation", class = "btn-primary")
    ),
    
    mainPanel(
      h3("Simulation Results"),
      plotlyOutput("corPlot", height = "500px"),
      verbatimTextOutput("resultText")
    )
  )
)


server <- function(input, output, session) {
  
  # =========================================================
  #  ShinyValidate --- 实时输入检查
  # =========================================================
  iv <- InputValidator$new()
  iv$add_rule("m", sv_gt(0, message = "Number of SNPs must be > 0"))
  iv$add_rule("n1", sv_gt(0, message = "n1 must be > 0"))
  iv$add_rule("n2", sv_gt(0, message = "n2 must be > 0"))
  iv$add_rule("w_min", sv_between(0, 1, message = "Overlap proportion min must be between 0 and 1"))
  iv$add_rule("w_max", sv_between(0, 1, message = "Overlap proportion max must be between 0 and 1"))
  iv$add_rule("w_step", sv_between(0, 1, message = "Overlap step must be between 0 and 1"))
  iv$add_rule("nsim", sv_gt(0, message = "Number of simulations per Overlap proportion must be > 0"))
  iv$enable()
  
  
  # =========================================================
  #  恢复默认参数按钮
  # =========================================================
  observeEvent(input$reset, {
    updateNumericInput(session, "m", value = 100)
    updateNumericInput(session, "n1", value = 1000)
    updateNumericInput(session, "n2", value = 500)
    updateNumericInput(session, "w_min", value = 0.1)
    updateNumericInput(session, "w_max", value = 0.5)
    updateNumericInput(session, "w_step", value = 0.1)
    updateNumericInput(session, "nsim", value = 10)
    updateNumericInput(session, "seed", value = 123)
  })
  
  
  # =========================================================
  #  主模拟逻辑
  # =========================================================
  
  bar_colors <- c("#4C4C4C", "#4DA6FF")
  
  results <- eventReactive(input$run, {
    
    # ---- 阻止无效输入运行 ----
    req(iv$is_valid())
    
    set.seed(input$seed)
    
    m <- input$m
    n1 <- input$n1
    n2 <- input$n2
    h2 <- 0.5
    
    w_values <- seq(input$w_min, input$w_max, by = input$w_step)
    observed_mean <- numeric(length(w_values))
    observed_sd <- numeric(length(w_values))
    theory_vals <- numeric(length(w_values))
    
    withProgress(message="Running simulations...", value=0, {
      for (iw in seq_along(w_values)) {
        
        w <- w_values[iw]
        sim_R2 <- numeric(input$nsim)
        
        for (sim in 1:input$nsim) {
          
          n_overlap <- ceiling(n2 * w)
          n <- n1 + n2 - n_overlap
          
          p <- runif(m, 0.1, 0.9)
          generate_row <- function(p) rbinom(length(p),1,p) + rbinom(length(p),1,p)
          x_matrix <- t(replicate(n, generate_row(p)))
          x_matrix <- scale(x_matrix)
          
          x_train <- x_matrix[1:n1, ]
          x_test  <- x_matrix[(n1+1-n_overlap):n, ]
          colnames(x_train) <- paste0("Col_", seq_len(ncol(x_train)))
          
          beta <- rnorm(m, 0, 1)
          G <- x_matrix %*% beta
          ei <- rnorm(n, 0, sqrt(var(G)/h2*(1-h2)))
          y <- G + ei
          y_test <- y[(n1+1-n_overlap):n]
          
          # 单 SNP 线性模型
          data_train_df <- data.frame(Y=y[1:n1], x_train)
          individual_snp_betas <- numeric(m)
          for (i in 1:m) {
            snp_i <- colnames(x_train)[i]
            fit <- lm(as.formula(paste("Y ~", snp_i)), data=data_train_df)
            individual_snp_betas[i] <- coef(fit)[2]
          }
          y_hat <- x_test %*% individual_snp_betas
          corObs <- cor(y_hat, y_test)
          sim_R2[sim] <- corObs^2
        }
        
        observed_mean[iw] <- mean(sim_R2)
        observed_sd[iw] <- sd(sim_R2)
        theory_vals[iw] <- (h2 + w*m/n1)^2 /
          (w*h2*m/n1 + (h2 + m/n1)*(1 + w*m/n1))
        
        incProgress(1/length(w_values), detail = paste("w =", round(w,3)))
      }
    })
    
    list(
      w_values = w_values,
      observed_mean = observed_mean,
      observed_sd = observed_sd,
      theory_vals = theory_vals
    )
  })
  
  
  # =========================================================
  #  文本输出
  # =========================================================
  output$resultText <- renderText({
    req(results())
    res <- results()
    paste(
      "Overlap proportion values:", paste(round(res$w_values,3), collapse=", "),
      "\n\nObserved R² mean:", paste(round(res$observed_mean,4), collapse=", "),
      "\nObserved R² sd:", paste(round(res$observed_sd,4), collapse=", "),
      "\nTheoretical R²:", paste(round(res$theory_vals,4), collapse=", ")
    )
  })
  
  
  # =========================================================
  #  交互式 Plotly 图
  # =========================================================
  output$corPlot <- renderPlotly({
    req(results())
    res <- results()
    
    plot_ly(x = round(res$w_values,3)) %>%
      add_trace(y=res$observed_mean, type="bar", name="Observed R²",
                marker=list(color=bar_colors[1]),
                error_y=list(type="data", array=res$observed_sd,
                             color="black", thickness=2)) %>%
      add_trace(y=res$theory_vals, type="bar", name="Theoretical R²",
                marker=list(color=bar_colors[2])) %>%
      layout(barmode="group",
             xaxis=list(title="Overlap proportion (w)"),
             yaxis=list(title="R²"),
             legend=list(x=0.99, y=1))
  })
  
  
  # =========================================================
  # 下载 CSV
  # =========================================================
  output$download_data <- downloadHandler(
    filename = function() paste0("PRS_simulation_data_", Sys.Date(), ".csv"),
    content = function(file){
      res <- results()
      write.csv(data.frame(
        Overlap_proportion=res$w_values,
        Observed_mean=res$observed_mean,
        Observed_sd=res$observed_sd,
        Theory=res$theory_vals
      ), file, row.names=FALSE)
    }
  )
  
  
  # =========================================================
  # 下载 PNG 图
  # =========================================================
  output$download_plot_png <- downloadHandler(
    filename = function() paste0("PRS_simulation_plot_", Sys.Date(), ".png"),
    content = function(file){
      res <- results()
      png(file, width=1400, height=900, res=150)
      
      bar_matrix <- rbind(res$observed_mean, res$theory_vals)
      colnames(bar_matrix) <- round(res$w_values,3)
      
      ylim_max <- max(res$observed_mean + res$observed_sd, res$theory_vals)*1.2
      
      bp <- barplot(bar_matrix, beside=TRUE, col=bar_colors,
                    ylim=c(0,ylim_max), ylab="R²", xlab="Overlap proportion (w)",
                    main="Observed R² vs Theoretical R²")
      
      arrows(bp[1,], res$observed_mean - res$observed_sd,
             bp[1,], res$observed_mean + res$observed_sd,
             angle=90, code=3, length=0.05, lwd=2, col="black")
      
      legend("topright", legend=c("Observed R²","Theoretical R²"),
             fill=bar_colors, bty="n")
      
      dev.off()
    }
  )
  
  
  # =========================================================
  # 下载 PDF 图（矢量图，高质量）
  # =========================================================
  output$download_plot_pdf <- downloadHandler(
    filename = function() paste0("PRS_simulation_plot_", Sys.Date(), ".pdf"),
    content = function(file){
      res <- results()
      pdf(file, width=10, height=7)
      
      bar_matrix <- rbind(res$observed_mean, res$theory_vals)
      colnames(bar_matrix) <- round(res$w_values,3)
      
      ylim_max <- max(res$observed_mean + res$observed_sd, res$theory_vals)*1.2
      
      bp <- barplot(bar_matrix, beside=TRUE, col=bar_colors,
                    ylim=c(0,ylim_max), ylab="R²",
                    xlab="Overlap proportion (w)",
                    main="Observed R² vs Theoretical R²")
      
      arrows(bp[1,], res$observed_mean - res$observed_sd,
             bp[1,], res$observed_mean + res$observed_sd,
             angle=90, code=3, length=0.05, lwd=2, col="black")
      
      legend("topright", legend=c("Observed R²","Theoretical R²"),
             fill=bar_colors, bty="n")
      
      dev.off()
    }
  )
}

shinyApp(ui = ui, server = server)
