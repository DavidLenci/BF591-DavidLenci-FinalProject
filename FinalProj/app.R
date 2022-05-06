# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(colourpicker)
library(DT)
library(ComplexHeatmap)
library(DESeq2)
library(igraph)


# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel("BF591 Final Project Spring 2022"),
    tabsetPanel(
        tabPanel("Sample", fluid = TRUE,
            sidebarLayout(
                sidebarPanel(fileInput("sample_file", "Input Sample Meta Data (csv)", accept = ".csv")
                ),
                mainPanel(
                    tabsetPanel(
                        tabPanel("Summary", fluid = TRUE, 
                                tableOutput("table1")),
                        tabPanel("Table", fluid = TRUE,
                                 DT::dataTableOutput("table2")),
                        tabPanel("Plot", fluid = TRUE, 
                                sidebarLayout(
                                    sidebarPanel(
                                        radioButtons(
                                            "var",
                                            "Select Variable",
                                            choices = "",
                                            inline = T
                                        )
                                    ),
                                    mainPanel(plotOutput("histogram"))
                                )
                        ))))),
        tabPanel("Counts", fluid=TRUE,
                 sidebarLayout(
                     sidebarPanel(fileInput("count_file", "Input Count Data (csv)", accept = ".csv"),
                                   sliderInput("var_slider", min = 0, max = 100,
                                               label ="Percent Variance", value = 75, step = 5),
                                   sliderInput("zero_slider", min = 1, max = 50,
                                               label ="Zeros", value = 50, step = 1)
                               ),
                 mainPanel(
                     tabsetPanel(
                         tabPanel("Summary", fluid=TRUE,
                                  tableOutput("sum_table2")),
                         tabPanel("Scatter", fluid=TRUE,
                                  plotOutput("med_vs_var"),
                                  plotOutput("med_vs_zero")),
                         tabPanel("Heatmap", fluid=TRUE,
                                  sidebarLayout(
                                      sidebarPanel(
                                          radioButtons(
                                              "logT",
                                              "Log Transformation?",
                                              choices = c("False", "True"),
                                              inline = T
                                          )
                                      ),
                                      mainPanel(plotOutput("heatmap"))
                                  )),
                         tabPanel("PCA", fluid=TRUE,
                                  sidebarLayout(
                                      sidebarPanel(
                                          radioButtons(
                                              "PC1",
                                              "Select PC",
                                              choices = c("1", "2", "3", "4",
                                                          "5", "6"),
                                              inline = T
                                          ),
                                          radioButtons(
                                              "PC2",
                                              "Select PC",
                                              choices = c("1", "2", "3", "4",
                                                          "5", "6"),
                                              selected = "2",
                                              inline = T
                                          )
                                      ),
                                      mainPanel(plotOutput("PCA"))
                                  )))))),
        tabPanel("Differential Expression", fluid = TRUE,
                 sidebarLayout(
                     sidebarPanel(fileInput("de_file", "Input DE Results (csv)", accept = ".csv"),
                                  hr(),
                                  actionButton("runDE", "Perform DE Using Norm Counts"),
                                  radioButtons(
                                      "de_var",
                                      "Select Variable for DE",
                                      choices = "",
                                      inline = T
                                  )
                                  
                     ),
                     mainPanel(
                         sidebarLayout(
                             sidebarPanel(
                                 radioButtons("var_x", "X Variable", 
                                              c("log2FoldChange", "stat", "lfcSE", "padj", "pvalue")),
                                 radioButtons("var_y", "Y Variable",
                                              c("padj", "pvalue", "log2FoldChange", "stat", "lfcSE")),
                                 sliderInput("slider", min = -300, max = 0,
                                             label ="significance", value = -100, step = 10),
                                 colourInput("nonSig", "Select Color 1", value = "black",
                                             showColour = "text"),
                                 colourInput("sig", "Select Color 2", value = "red",
                                             showColour = "text")
                             ),
                             mainPanel(
                                 tabsetPanel(
                                     tabPanel("Volcano Plot", plotOutput("volcano")),
                                     tabPanel("Table", tableOutput("de_table"))
                                 )))))),
        tabPanel("Correlation Network Analysis",fluid = TRUE,
                 sidebarLayout(
                     sidebarPanel(fileInput("cor_file", "Input Norm Counts (csv)", 
                                            accept = ".csv"),
                                  textAreaInput("cor_genes", "Genes of Interest", 
                                                width = "1000px"),
                                  actionButton("runCor", 
                                               "Perform Correlation Analysis"),
                                  sliderInput("cor_slider", min = 0, max = 1,
                                              label ="correlation", value = 0.5, step = 0.1)),
                     mainPanel(
                         tabsetPanel(
                             tabPanel("Heatmap", plotOutput("cor_hm")),
                             tabPanel("Correlation Graph", plotOutput("cor_graph")),
                             tabPanel("Metrics", tableOutput("cor_metrics"))
                         ))))
        
        
    )
)
                        
        


                        
# Define server logic required to draw a histogram
server <- function(input, output, session) {
    options(shiny.maxRequestSize=30*1024^2)
    #Samples
    #-------------------------------------------------------------------------->
    load_data <- eventReactive(input$sample_file, {
        
        if(!is.null(input$sample_file$datapath)){
            data <- read.csv(input$sample_file$datapath,
                            header = T)
        }else{
            return(NULL)
        }
        
        if(!is.null(data)){
            observe({
                vchoices <- get_radio_names(data)
                updateRadioButtons(session, "var", choices = vchoices)
            })
            
            observe({
                dechoices <- get_de_names(data)
                updateRadioButtons(session, "de_var", choices = dechoices)
            })
            return(data)
        }
    })
    helper <- function(x){
        
        type <- class(x)
        value = ""
        
        if(type=="numeric"||type=="integer"){
            av <- as.character(format(round(mean(na.omit(x)), 2), nsmall = 2))
            sd <- as.character(format(round(sd(na.omit(x)), 2), nsmall = 2))
            value <- str_c(c(av,"(+/-", sd, ")"), collapse = "")
        }
        if(type=="character"){
            value <- unique(x)
            value <- str_c(value, collapse=", ")
        }
        return(list(type, value))
    }
    
    summary_tb <- function(df_in){
        df <- df_in
        types <- lapply(df, function(x){helper(x)[[1]]})
        values <- lapply(df, function(x){helper(x)[[2]]})
        
        types <- as.data.frame(types)
        values <- as.data.frame(values)
        
        dataf <- as.data.frame(t(rbind(types, values)))
        dataf <- dataf %>%
            filter(nchar(V2)<50)
        colnames(dataf) <- c("Type", "Mean (sd) or Distinct Values")
        dataf <- tibble::rownames_to_column(dataf, "Column Names")
        
        return(dataf)
    }
    
    get_radio_names <- function(dataf){
        dataf <- dataf %>%
            select_if(function(x) is.integer(x) || is.numeric(x))
        names <- colnames(dataf)
        return(names)
    }
    
    sample_histogram <- function(dataf, var) {
        #clean data
        plot <- dataf %>%
            ggplot(aes(x=!!sym(var))) + 
            geom_histogram(stat="count") +
            theme_bw()
        return(plot)
    }
    #counts
    #-------------------------------------------------------------------------->
    
    load_counts <- eventReactive(input$count_file, {
        
        if(!is.null(input$count_file$datapath)){
            data <- read.csv(input$count_file$datapath,
                             header = T)
            rownames(data) <- data[,1]
            data <- data[,-1]
            return(data)
        }else{
            return(NULL)
        }
    })
    
    filter_data <- function(data, thresh, zeros){
        thresh <- thresh/100
        variances <- apply(data, 1, var)
        thresh <- quantile(variances, probs = thresh)
        
        keep <- apply(data, 1, function(x){var(x)>thresh&sum(x==0)<zeros})
        data <- data[keep,]
        
        return(data)
    }
    
    summary_counts <- function(data, thresh, zeros){
        filter <- filter_data(data, thresh, zeros)
        #filter/get final dims
        initial_genes <- nrow(data)
        samples <- as.character(ncol(data))
        filter_genes <- nrow(filter)
        perc <- as.character(format(round((filter_genes/initial_genes)*100), 2), nsmall = 2)
        no_pass <- initial_genes-filter_genes
        no_pass_perc <- as.character(format(round((no_pass/initial_genes)*100), 2), nsmall =2)
        
        summary <- list(c("Genes", "Samples", "Passing Filter", "Percent Passing Filter",
                          "Not Passing Filter", "Percent Not Passing Filter"),
                        c(initial_genes, samples, filter_genes, perc, no_pass, no_pass_perc))
        
        summary <- as.data.frame(summary)
        colnames(summary) <- NULL
        
        return(summary)
    }
    
    scatter_plot <- function(data, thresh, zeros_n){
        thresh <- thresh/100
        variances <- apply(data, 1, var)
        thresh <- quantile(variances, probs = thresh)
        
        zeros <- apply(data, 1, function(x){sum(x==0)})
        medians <- apply(data, 1, median)
        
        plotData <- tibble(variance = variances,
                           zero = zeros,
                           median = medians)
        
        plot1 <- plotData %>%
            ggplot(aes(x = log(median+1), y = log(variance), colour = variance<thresh)) +
            geom_point(stat = "identity", show.legend = FALSE) +
            scale_colour_manual(values = c("black", "grey")) +
            ylab("Variance") +
            xlab("log(Median+1)") +
            ggtitle("log(Median Norm Counts) vs log(Variance Norm Counts)")+
            theme_minimal()
        
        plot2 <- plotData %>%
            ggplot(aes(x = log(median+1), y = zero, colour = zero>zeros_n)) +
            geom_point(stat = "identity", show.legend = FALSE) +
            scale_colour_manual(values = c("black", "grey")) +
            ylab("Number of Zeros") +
            xlab("log(Median+1)") +
            ggtitle("log(Median Norm Counts) vs Number of Zero Values")+
            theme_minimal()
        
        return(list(plot1, plot2))
    }
    
    heatmap_f <- function(data, thresh, zeros_n, logT){
        filter_df <- filter_data(data, thresh, zeros_n)
        
        plotData <- as.matrix(filter_df)
        
        if(logT == "False"){
            h <- Heatmap(plotData, name = 'Counts', cluster_rows = FALSE,
                         column_gap = unit(0, "mm"), 
                         show_row_names = FALSE, show_column_names = FALSE, 
                         border = TRUE)
        }else{
            log_plot_data <- log(plotData+1)
            
            h <- Heatmap(log_plot_data, name = 'log(Counts)', cluster_rows = FALSE,
                         column_gap = unit(0, "mm"), 
                         show_row_names = FALSE, show_column_names = FALSE, 
                         border = TRUE)
            
        }
        return(h)
    }
    
    PCA <- function(data, thresh, zeros_n, X, Y){
        X_val <- as.integer(X)
        Y_val <- as.integer(Y)
        
        filter_df <- filter_data(data, thresh, zeros_n)
        pca_results <- prcomp(filter_df)
        #add metadata to pca output for access while plotting
        pca_data <- as.data.frame(pca_results$x)
        #calculate percent variances for each PC
        variances <- sapply(pca_results$sdev, function(x) x^2 )
        perc_variances <- sapply(variances, function(x) x/sum(variances))
        
        #generate x and y labels
        PC1_var <- as.character(round(perc_variances[X_val]*100))
        PC2_var <- as.character(round(perc_variances[Y_val]*100))
        
        y_label <- str_c("PC", Y, ": ", PC2_var, "% variance")
        x_label <- str_c("PC", X, ": ", PC1_var, "% variance")
        #plot PC1 and 2 using ggplot
        
        pc1 <- str_c("PC", X)
        pc2 <- str_c("PC", Y)
        plot <- pca_data %>%
            ggplot(aes(x=as.numeric(!!sym(pc1)), y=as.numeric(!!sym(pc2))))+
            geom_point(stat="identity")+
            xlab(x_label)+
            ylab(y_label) +
            theme_minimal()
        
        return(plot)
    }
    
    #DE
    #-------------------------------------------------------------------------->
    
    load_de <- eventReactive(input$de_file, {
        
        if(!is.null(input$de_file$datapath)){
            data <- read.csv(input$de_file$datapath,
                             header = T)
            rownames(data) <- data[,1]
            data <- data[,-1]
            return(data)
        }else{
            return(NULL)
        }
    })
    
    get_de_names <- function(dataf){
        dataf <- dataf %>%
            select_if(function(x) is.character(x))
        names <- colnames(dataf)
        return(names)
    }
    
    run_deseq <- function(count_dataframe, meta, condition){
        
        samples <- colnames(count_dataframe)
        design <- meta[,condition]
        genes <- rownames(count_dataframe)
        
        coldata <- data.frame(samples=samples,
                              diagnosis=design)
        
        count_dataframe <- as.data.frame(apply(count_dataframe, 2, as.integer))
        #gen deseq object
        dds <- DESeqDataSetFromMatrix(countData = count_dataframe,
                                      colData = coldata,
                                      design = ~ diagnosis)
        
        #perform de analysis
        dds <- DESeq(dds)
        res <- as.data.frame(results(dds))
        rownames(res) <- genes
        
        return(res)
    }
    
    
    volcano_plot <-
        function(dataf, x_name, y_name, slider, color1, color2) {
            #clean data
            dataf <- dataf %>%
                drop_na(!!sym(x_name)) %>%
                drop_na(!!sym(y_name))
            #create plot
            plot <- dataf %>%
                ggplot(aes(!!sym(x_name), -log10(!!sym(y_name)), 
                           colour = log10(!!sym(y_name)) < slider)) +
                geom_point(stat="identity")+
                scale_colour_manual(values=c(color1, color2)) +
                xlab(x_name)+
                ylab(y_name)+
                theme_bw()
            
            return(plot)
        }
    
    draw_table <- function(dataf, slider) {
        #clean data
        dataf <- dataf %>%
            drop_na(padj) %>%
            filter(log10(padj)<slider)%>%
            arrange(padj)
        
        dataf$padj <- formatC(dataf$padj, digits = 5)
        dataf$pvalue <- formatC(dataf$pvalue, digits = 5)
        
        return(dataf)
    }
    
    #Correlation Plot
    #-------------------------------------------------------------------------->
    
    load_cor <- eventReactive(input$cor_file,{
        
        if(!is.null(input$cor_file$datapath)){
            data <- read.csv(input$cor_file$datapath,
                             header = T)
            rownames(data) <- data[,1]
            data <- data[,-1]
            return(data)
        }else{
            return(NULL)
        }
    })
    
    filter_cor <- function(data, genes){
        genes <- strsplit(genes, "\n")
        genes <- genes[[1]]
        data <- data[genes,]
        
        return(data)
    }
    
    heatmap_cor <- function(data, genes){
        data_s <- filter_cor(data, genes)
        plotData <- as.matrix(data_s)
        
        h <- Heatmap(plotData, name = 'Counts', cluster_rows = FALSE,
                    column_gap = unit(0, "mm"), 
                    show_row_names = FALSE, show_column_names = FALSE, 
                    border = TRUE)
        
        return(h)
    }
    
    corr_network <- function(data, thresh, genes){
        data <- filter_cor(data, genes)
        mat <- cor(t(data))
        
        mat[mat<thresh] <- 0
        
        # Make an Igraph object from this matrix:
        network <- graph_from_adjacency_matrix( mat, weighted=T, 
                                                mode="undirected", diag=F)
        
        # Basic chart
        plot <- plot(network)
        
        stats <- data.frame(degree = degree(network),
                            closeness = closeness(network),
                            betweenness = betweenness(network))
        
        return(list(plot, stats))
        
    }
    
    
    
    #Output
    #-------------------------------------------------------------------------->
    #Samples
    output$table1 <- renderTable({summary_tb(load_data())},
                                bordered = T)
    
    output$table2 <- DT::renderDataTable(load_data(),
                                         options = list(scrollX = TRUE),
                                         rownames = FALSE)
    
    output$histogram <- renderPlot({sample_histogram(load_data(), input$var)}) 
    
    #Counts
    output$sum_table2 <- renderTable({summary_counts(load_counts(),
                                                     input$var_slider,
                                                     input$zero_slider)},
                                     bordered=T)
    
    output$med_vs_var <- renderPlot({scatter_plot(load_counts(),
                                                  input$var_slider,
                                                  input$zero_slider)[1]})
    
    output$med_vs_zero <- renderPlot({scatter_plot(load_counts(),
                                                  input$var_slider,
                                                  input$zero_slider)[2]})
    
    
    output$heatmap<-renderPlot({heatmap_f(load_counts(),
                                          input$var_slider,
                                          input$zero_slider,
                                          logT=input$logT)})
    
    output$PCA <- renderPlot({PCA(load_counts(),
                                  input$var_slider,
                                  input$zero_slider,
                                  input$PC1,
                                  input$PC2)})
    
    #DE
    
    observeEvent(input$runDE,{
        DESeq_result <- run_deseq(load_counts(), load_data(), input$de_var)
        
        output$volcano <- renderPlot({volcano_plot(DESeq_result, 
                                                   input$var_x, 
                                                   input$var_y,
                                                   input$slider,
                                                   input$nonSig,
                                                   input$sig)
        })
        
        output$de_table <- renderTable({draw_table(DESeq_result, input$slider)},
                                    bordered = T) 
        
        
    })
    
    output$volcano <- renderPlot({volcano_plot(load_de(), 
                                               input$var_x, 
                                               input$var_y,
                                               input$slider,
                                               input$nonSig,
                                               input$sig)
    })
    
    output$de_table <- renderTable({draw_table(load_de(), input$slider)},
                                bordered = T) 
    
    #Correlation
    
    observeEvent(input$runCor, {
        output$cor_hm <- renderPlot({heatmap_cor(load_cor(),
                                                input$cor_genes)
            })
    
        output$cor_graph <- renderPlot({corr_network(load_cor(),
                                                    input$cor_slider,
                                                    input$cor_genes)[1]})
    
        output$cor_metrics <- renderTable({corr_network(load_cor(),
                                                        input$cor_slider,
                                                        input$cor_genes)[2]}, 
                                          bordered = T)
    })
    
}

# Run the application
shinyApp(ui = ui, server = server)
