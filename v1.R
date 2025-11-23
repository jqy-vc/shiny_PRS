library(shiny)
library(plotly)

ui <- fluidPage(
  titlePanel("PRS Simulation Tool"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("m", "Number of SNPs (loci):", value = 100, min = 1, step = 1),
      
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
  
  # 统一颜色（Plotly + PDF + PNG）
  bar_colors <- c("#4C4C4C", "#4DA6FF")
  
  # ---- 输入检查（最推荐方式）-----
  observeEvent(input$run, {
    validate(
      need(input$m > 0, "Error: Number of SNPs (m) must be > 0."),
      need(input$n1 > 0, "Error: Training sample size (n1) must be > 0."),
      need(input$n2 > 0, "Error: Testing sample size (n2) must be > 0."),
      need(input$w_min >= 0, "Overlap proportion min must be >= 0."),
      need(input$w_max >= 0, "Overlap proportion max must be >= 0."),
      need(input$w_step > 0, "Overlap step must be > 0."),
      need(input$nsim > 0, "Number of simulations must be > 0.")
    )
  })
  
  # ---- 模拟结果 ----
  results <- eventReactive(input$run, {
    
    set.seed(input$seed)
    
    m <- input$m
    n1 <- input$n1
    n2 <- input$n2
    h2 <- 0.5
    
    w_values <- seq(input$w_min, input$w_max, by = input$w_step)
    
    observed_mean <- numeric(length(w_values))
    observed_sd <- numeric(length(w_values))
    theory_vals <- numeric(length(w_values))
    
    withProgress(message = "Running simulations...", value = 0, {
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
          
          beta <- rnorm(m,0,1)
          G <- x_matrix %*% beta
          ei <- rnorm(n,0,sqrt(var(G)/h2*(1-h2)))
          y <- G + ei
          y_test <- y[(n1+1-n_overlap):n]
          
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
        
        incProgress(1/length(w_values), detail=paste("w =", round(w,3)))
      }
    })
    
    list(
      w_values = w_values,
      observed_mean = observed_mean,
      observed_sd = observed_sd,
      theory_vals = theory_vals
    )
  })
  
  # ---- 文本输出 ----
  output$resultText <- renderText({
    req(results())
    res <- results()
    paste(
      "Overlap proportion values:", paste(round(res$w_values, 3), collapse=", "),
      "\n\nObserved R² mean:", paste(round(res$observed_mean,4), collapse=", "),
      "\nObserved R² sd:", paste(round(res$observed_sd,4), collapse=", "),
      "\nTheoretical R²:", paste(round(res$theory_vals,4), collapse=", ")
    )
  })
  
  # ---- 交互式 Plotly 图 ----
  output$corPlot <- renderPlotly({
    req(results())
    res <- results()
    
    plot_ly(x = round(res$w_values,3)) %>%
      add_trace(y = res$observed_mean, type='bar', name='Observed R²',
                error_y=list(type='data', array=res$observed_sd,
                             color='black', thickness=2, width=8),
                marker=list(color=bar_colors[1])) %>%
      add_trace(y = res$theory_vals, type='bar', name='Theoretical R²',
                marker=list(color=bar_colors[2])) %>%
      layout(barmode='group',
             xaxis=list(title='Overlap proportion (w)'),
             yaxis=list(title='R²'),
             legend=list(
               x=1, y=1,
               bgcolor="rgba(255,255,255,0.8)"
             ))
  })
  
  # ---- 下载 CSV ----
  output$download_data <- downloadHandler(
    filename = function() paste0("PRS_simulation_data_", Sys.Date(), ".csv"),
    content = function(file){
      res <- results()
      df <- data.frame(
        w = res$w_values,
        observed_mean = res$observed_mean,
        observed_sd = res$observed_sd,
        theory = res$theory_vals
      )
      write.csv(df, file, row.names=FALSE)
    }
  )
  
  # ---- 下载 PNG ----
  output$download_plot_png <- downloadHandler(
    filename = function() paste0("PRS_simulation_plot_", Sys.Date(), ".png"),
    content = function(file){
      res <- results()
      
      png(file, width=1200, height=800, res=120)
      
      bar_matrix <- rbind(res$observed_mean, res$theory_vals)
      colnames(bar_matrix) <- round(res$w_values, 3)
      
      ylim_max <- max(res$observed_mean + res$observed_sd, res$theory_vals) * 1.2
      
      bp <- barplot(
        bar_matrix,
        beside = TRUE,
        col = bar_colors,
        ylim = c(0, ylim_max),
        ylab = "R²",
        xlab = "Overlap proportion (w)",
        main = "Observed R² vs Theoretical R²"
      )
      
      # ---- 加误差棒 ----
      arrows(
        bp[1,], res$observed_mean - res$observed_sd,
        bp[1,], res$observed_mean + res$observed_sd,
        angle = 90, code = 3, length = 0.05, lwd = 2, col = "black"
      )
      
      # ---- 图例（Legend）----
      legend(
        "topright",
        legend = c("Observed R²", "Theoretical R²"),
        fill = bar_colors,
        border = NA,
        bty = "n",
        cex = 1.2
      )
      
      dev.off()
    }
  )
  
  
  # ---- 下载 PDF（高质量矢量图）----
  output$download_plot_pdf <- downloadHandler(
    filename = function() paste0("PRS_simulation_plot_", Sys.Date(), ".pdf"),
    content = function(file){
      res <- results()
      
      pdf(file, width=10, height=7)
      
      bar_matrix <- rbind(res$observed_mean, res$theory_vals)
      colnames(bar_matrix) <- round(res$w_values, 3)
      
      ylim_max <- max(res$observed_mean + res$observed_sd, res$theory_vals) * 1.2
      
      bp <- barplot(
        bar_matrix,
        beside = TRUE,
        col = bar_colors,
        ylim = c(0, ylim_max),
        ylab = "R²",
        xlab = "Overlap proportion (w)",
        main = "Observed R² vs Theoretical R²"
      )
      
      # ---- 加误差棒 ----
      arrows(
        bp[1,], res$observed_mean - res$observed_sd,
        bp[1,], res$observed_mean + res$observed_sd,
        angle = 90, code = 3, length = 0.05, lwd = 2, col = "black"
      )
      
      # ---- 图例（Legend）----
      legend(
        "topright",
        legend = c("Observed R²", "Theoretical R²"),
        fill = bar_colors,
        border = NA,
        bty = "n",
        cex = 1.2
      )
      
      dev.off()
    }
  )
  
}

shinyApp(ui = ui, server = server)
