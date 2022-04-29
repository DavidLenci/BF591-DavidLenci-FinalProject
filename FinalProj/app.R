# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(colourpicker)
library(DT)# you might need to install this.


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
                                               label ="Percent Variance", value = 0, step = 5),
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
                                  plotOutput("heatmap")),
                         tabPanel("PCA", fluid=TRUE,
                                  sidebarLayout(
                                      sidebarPanel(
                                          radioButtons(
                                              "X",
                                              "Select PC",
                                              choices = c("PC1", "PC2", "PC3", "PC4"),
                                              inline = T
                                          ),
                                          radioButtons(
                                              "Y",
                                              "Select PC",
                                              choices = c("PC1", "PC2", "PC3", "PC4"),
                                              inline = T
                                          )
                                      ),
                                      mainPanel(plotOutput("PCA"))
                                  ))))))
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
    
    summary_tb <- function(df){
        types <- as.data.frame(lapply(df, function(x){helper(x)[[1]]}))
        values <- as.data.frame(lapply(df, function(x){helper(x)[[2]]}))
        
        dataf <- as.data.frame(t(rbind(types, values)))
        dataf <- dataf %>%
            filter(nchar(V2)<30)
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
            ylab("Number of Zeros") +
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
}

# Run the application
shinyApp(ui = ui, server = server)
