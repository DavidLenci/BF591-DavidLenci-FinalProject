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
                                               label ="Percent Variance", value = 50, step = 1),
                                   sliderInput("zero_slider", min = 1, max = 100,
                                               label ="Allowed Zeros Counts per Gene", 
                                               value = 100, step = 1),
                                  radioButtons(
                                      "logT",
                                      "Log Transformation for Heatmap:",
                                      choices = c("False", "True"),
                                      inline = T
                                  )
                               ),
                 mainPanel(
                     tabsetPanel(
                         tabPanel("Summary", fluid=TRUE,
                                  tableOutput("sum_table2")),
                         tabPanel("Scatter", fluid=TRUE,
                                  plotOutput("med_vs_var"),
                                  plotOutput("med_vs_zero")),
                         tabPanel("Heatmap", fluid=TRUE,
                                  plotOutput("heatmap")
                                  ),
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
                         tabsetPanel(
                             tabPanel("DE Data", 
                                      DT::dataTableOutput("DE_table")),
                             tabPanel("DE Analysis",
                                      sidebarLayout(
                                          sidebarPanel(
                                              radioButtons("var_x", "X Variable", 
                                                           c("log2FoldChange", "stat", "lfcSE", "padj", "pvalue")),
                                              radioButtons("var_y", "Y Variable",
                                                           c("padj", "pvalue", "log2FoldChange", "stat", "lfcSE")),
                                              sliderInput("slider", min = -150, max = 0,
                                                           label ="significance", value = -50, step = 5),
                                              colourInput("nonSig", "Select Color 1", value = "black",
                                                           showColour = "text"),
                                              colourInput("sig", "Select Color 2", value = "red",
                                                           showColour = "text")
                                          ),
                                          mainPanel(
                                              tabsetPanel(
                                                  tabPanel("Volcano Plot", plotOutput("volcano")),
                                                  tabPanel("Table", DT::dataTableOutput("de_table"))
                                 )))))))),
        tabPanel("Correlation Network Analysis",fluid = TRUE,
                 sidebarLayout(
                     sidebarPanel(fileInput("cor_file", "Input Norm Counts (csv)", 
                                            accept = ".csv"),
                                  textAreaInput("cor_genes", "Genes of Interest", 
                                                width = "1000px"),
                                  actionButton("runCor", 
                                               "Perform Correlation Analysis"),
                                  sliderInput("cor_slider", min = 0, max = 1,
                                              label ="correlation", value = 0.5, step = 0.05),
                                  verbatimTextOutput("b_genes")),
                     mainPanel(
                         tabsetPanel(
                             tabPanel("Heatmap", plotOutput("cor_hm")),
                             tabPanel("Correlation Graph", plotOutput("cor_graph")),
                             tabPanel("Metrics", tableOutput("cor_metrics"))
                         )))),
        
        tabPanel("Help Page", fluid = T,
                 mainPanel(
                     h2("Sample Tab"),
                     h3("Data Input"),
                     p("- The expected input for this section of the shiny app is
                       an associated meta data file. The rows are expected to be 
                       samples and the columns are variables describing the data. There
                       are not expectations for the contents of this meta data file.
                       Although, if the users wishes to perform DESeq2 through the 
                       shiny app then there weill need to be atleast one column 
                       that can be used to separate the samples into groups for
                       the anlysis (a design parameter to give to DESeq). Lastly,
                       the input must be in csv format."),
                     br(),
                     h3("Tabs"),
                     p("- The sample tab allows users to get a high level analysis
                       of any meta data associated with their sequencing data. The
                       expected input is a csv file with rows as the samples and
                       columns the different variables."), 
                     p("- The first tab will generate a table summarizing the contents
                       of each column in the meta data. "),
                     p("- The table tab generates a data table that allows users
                       to sort by specific columns and search the table."),
                     p("- The plot tab will allow the user to generate histograms
                       for some of the data from the meta data. These columns are
                       chosen based on the number of unique values in that column.
                       If a column has too many unique values (equal to the number 
                       of samples) or too little unique values then it is not 
                       included as an option."),
                     br(),
                     br(),
                     hr(),
                     h2("Counts Tab"),
                     h3("Data Input"),
                     p("- Similarly to the sample input, here the app expects a
                       csv file, but with the rows as genes and the columns as 
                       samples. The counts data shouldn't consist of any character
                       values and should only consist of numeric values. "),
                     br(),
                     h3("Tabs"),
                     p("- The summary tab for the counts section generates a table 
                       that tells the user high level information about the count
                       file they uploaded. Specifically, it allows them to see what
                       affect the filters they are applying through the corresponding
                       sliders is having on their data. i.e what percent of genes
                       in the counts file is passing and not passing the applied
                       filters."),
                     p("- The scatter plot tab generates two plots from the input
                       counts file. Both are scatter plots; the first shows the
                       variance vs the median counts and the second plot shows the number
                       of samples with zero counts for a gene vs the median counts. 
                       Genes that pass the applied filters are black, and those that 
                       do not are a light grey. The figures are reactive to changes
                       in filter values."),
                     p("- The next tab generates a heatmap clustering samples based
                       on the count data that is passing the applied filters. There
                       is also an option to apply a log transformation to the data
                       prior to generating the heatmap."),
                     p("- The PCA tab allows users to perform principle component
                       analysis on the counts that passed the applied filter. 
                       Users have the option to select from the first 5 principle
                       components for the x and y axis. "),
                     br(),
                     br(),
                     hr(),
                     h2("Differential Expression Tab"),
                     h3("Data Input"),
                     p("- For this input the expectations of the rows and, more
                       specifically, column variables is more stringent. The file
                       should have a log 2 Fold change column, p-value column, 
                       p-value adjusted column, log fold change standard error coloum,
                       and a statistic column. The expected names for these
                       columns are log2FoldChange, pvalue, padj, lfcSE, and stat
                       respectively. This is the exact output from DESeq2 currently,
                       and is what will be generated if the user instead chooses to
                       generate the DE results with the uploaded counts and sample
                       data files."),
                     br(),
                     h3("Tabs"),
                     p("- The first tab generates a scatter plot of the log2 Fold
                       Change values, and colors points based on whether or not
                       they pass the signficance threshold applied in the slider
                       on the side of the tab."),
                     p("- The table tab generates a table with DE results for 
                       genes that passed the significance threshold."),
                     br(),
                     br(),
                     hr(),
                     h2("Correlation Analysis Tab"),
                     h3("Data Input"),
                     p("- For this tab users input the same counts file tab, or
                       a different tab following the description in the counts
                       tab data input description. Additionally, users can 
                       specify there genes of interest in generating the correlation
                       analysis results. The analysis can not be run without specifying
                       genes of interest. For this text input each gene is on a 
                       different line in the input box."),
                     h3("Tabs"),
                     p("- Generates a heatmap clustering samples based on the expression
                       for the genes specified. "),
                     p("- Takes the genes specified in the text input, and generates
                       a correlation matrix to be used for creating a correlation
                       network graph."),
                     p("- Pulls closeness, degree, and betweeness values from the 
                       network and presents them in a table."),
                     br(),
                     br()
                 ))
        
        
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
        sample_n <- length(x)
        
        if(type=="numeric"||type=="integer"){
            av <- as.character(format(round(mean(na.omit(x)), 2), nsmall = 2))
            sd <- as.character(format(round(sd(na.omit(x)), 2), nsmall = 2))
            value <- str_c(c(av,"(+/-", sd, ")"), collapse = "")
        }
        if(type=="character"){
            value <- unique(x)
            if(length(value)==sample_n){
                value = "Identifier"
            }else{
                value <- str_c(value, collapse=", ")
            }
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
        
        colnames(dataf) <- c("Type", "Mean (sd) or Distinct Values")
        dataf <- tibble::rownames_to_column(dataf, "Column Names")
        
        return(dataf)
    }
    
    get_radio_names <- function(dataf){
        dataf <- dataf %>%
            select_if(function(x) is.integer(x) || is.numeric(x)) %>%
            select_if(function(x) length(unique(x)) < (0.9 * nrow(dataf)) & length(unique(x)) > 2)
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
                         show_row_names = FALSE, 
                         border = TRUE)
        }else{
            log_plot_data <- log(plotData+1)
            
            h <- Heatmap(log_plot_data, name = 'log(Counts)', cluster_rows = FALSE,
                         column_gap = unit(0, "mm"), 
                         show_row_names = FALSE, #show_column_names = FALSE, 
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
            select_if(function(x) is.character(x)) %>%
            select_if(function(x) length(unique(x)) < 5 & length(unique(x)) > 1)
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
        
        genes <- intersect(genes, rownames(data))
        genes <- unlist(genes)
    
        data <- data[genes,]
        
        return(data)
    }
    
    heatmap_cor <- function(data, genes){
        data_s <- filter_cor(data, genes)
        plotData <- as.matrix(data_s)
        
        h <- Heatmap(log(plotData+1), name = 'log(counts)', cluster_rows = FALSE,
                    column_gap = unit(0, "mm"), 
                    #show_row_names = FALSE, show_column_names = FALSE, 
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
        plot <- plot(network,
                     vertex.size=18,
                     vertex.color="lightblue",
                     vertex.label.cex = 1,
                     vertex.label.dist = 2,
                     vertex.label.color = "black")
        
        genes <- strsplit(genes, "\n")
        genes <- genes[[1]]
        genes <- intersect(genes, rownames(data))
        
        stats <- data.frame(genes = genes,
                            degree = degree(network),
                            closeness = closeness(network),
                            betweenness = betweenness(network))
        
        return(list(plot, stats))
        
    }
    
    find_genes <- function(data, gene){
        
        if(gene %in% rownames(data)){
        }
        else{
            return(gene)
        }
        
    }
    
    bad_genes <- function(data, genes){
        
        genes <- strsplit(genes, "\n")
        genes <- genes[[1]]
        
        b_genes <- lapply(genes, function(x){find_genes(data, x)})
        
        b_genes <- b_genes %>% 
            discard(is.null)
        
        if(length(b_genes)>0){
            b_genes <- str_c(b_genes, collapse=", ")
            return(str_c("Genes not in data: ", b_genes))
        }
        else{
            return("")
        }
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
        
        output$DE_table <- DT::renderDataTable(DESeq_result,
                                               options = list(scrollX = TRUE),
                                               rownames = TRUE)
        
        output$volcano <- renderPlot({volcano_plot(DESeq_result, 
                                                   input$var_x, 
                                                   input$var_y,
                                                   input$slider,
                                                   input$nonSig,
                                                   input$sig)
        })
        
        output$de_table <- DT::renderDataTable(draw_table(DESeq_result, input$slider),
                                               options = list(scrollX = TRUE),
                                               rownames = TRUE)
        
        
        
    })
    
    output$DE_table <- DT::renderDataTable(load_de(),
                                           options = list(scrollX = TRUE),
                                           rownames = TRUE)
    
    output$volcano <- renderPlot({volcano_plot(load_de(), 
                                               input$var_x, 
                                               input$var_y,
                                               input$slider,
                                               input$nonSig,
                                               input$sig)
    })
    
    output$de_table <- output$de_table <- DT::renderDataTable(draw_table(load_de(), input$slider),
                                                              options = list(scrollX = TRUE),
                                                              rownames = TRUE)
    
    #Correlation
    
    observeEvent(input$runCor, {
        
        output$b_genes <- renderText({bad_genes(load_cor(),
                                                input$cor_genes)})
        
        output$cor_hm <- renderPlot({heatmap_cor(load_cor(),
                                                input$cor_genes)
            })
    
        output$cor_graph <- renderPlot({corr_network(load_cor(),
                                                    input$cor_slider,
                                                    input$cor_genes)[1]},
                                       width = 750,
                                       height = 750)
    
        output$cor_metrics <- renderTable({corr_network(load_cor(),
                                                        input$cor_slider,
                                                        input$cor_genes)[2]}, 
                                          bordered = T)
    })
    
}

# Run the application
shinyApp(ui = ui, server = server)
