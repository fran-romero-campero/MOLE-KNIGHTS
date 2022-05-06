
## Increase max file size allowed to upload to 100MB
options(shiny.maxRequestSize=100*1024^2)

# Paquetes
library(shiny)
library(ggplot2, logical.return = T)
#install.packages("MetBrewer")
library(MetBrewer, logical.return = T)
#install.packages("cowplot")
library(cowplot, logical.return = T)
library(plotly)
library(ggrepel)
library(DT)
library(limma)

# Funciones:
barplot.gene <- function(gene.id,gene.name,gene.expression)
{
  expression.ll <- mean(unlist(gene.expression[,c("LL_1", "LL_2")]))
  expression.hl <- mean(unlist(gene.expression[,c("HL_1", "HL_2")]))
  means <- c(expression.ll,expression.hl)

  expression.sd.ll <- sd(unlist(gene.expression[,c("LL_1", "LL_2")]))
  expression.sd.hl <- sd(unlist(gene.expression[,c("HL_1", "HL_2")]))
  sds <- c(expression.sd.ll,expression.sd.hl)


  par(lwd=3)
  xpos <- barplot(means,col=c("blue","firebrick2"),
                  names.arg = c("LL","HL"),las=2,cex.names = 1.5,
                  ylim=c(0,max(means+sds)*1.2),cex.axis = 1.5,lwd=3,
                  main=paste(c(gene.name, "-", gene.id),collapse=" "),
                  cex.main=2)
  arrows(x0 = xpos,y0 = means + sds, x1 = xpos, y1 = means - sds,length = 0.1,
         code = 3,angle=90,lwd=2)
  return(list(means,means[2]/means[1]))
}



# Aplicación para boxplots
ui <- fluidPage(

    # Título
    titlePanel("Boxplot for RNAseq"),

        # Show a plot of the generated distribution
        mainPanel(align="center",
           plotlyOutput("distPlot"),
           plotlyOutput("Volcano"),
           tags$br(), tags$br(),
           tags$br(), tags$br(),
           tags$br(), tags$br(),
           tags$br(), tags$br(),
           tags$br(), tags$br(),
           tags$br(), tags$br(),
           tags$br(), tags$br(),
           DT::dataTableOutput("table"),
           plotOutput("Barplot")
        )
        
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
    # Cargamos los datos y asjutamos los datos
    gene.expression_total <- read.table(file = "gene_expression.tsv",header=T,as.is=T)
    gene.ids <- gene.expression_total$geneID
    gene.expression <- as.matrix(gene.expression_total[,2:ncol(gene.expression_total)])
    rownames(gene.expression) <- gene.ids
    
    log.gene.expression <-read.table(file = "log_gene-expression.tsv",header=T,as.is=T)
    de.results<-read.table(file = "de_results.tsv",header=T,as.is=T)
    resultado_final<-read.table(file = "resultado_final.tsv",header=T,as.is=T)
    head(resultado_final)
  

    # Variable que almacena los colores a usar en el volcano plot
    volcol <- c(met.brewer("Egypt",3)[1], met.brewer("Egypt",3)[2],"grey33")
    
    # Se asigna a cada uno de los colores los posibles niveles de expresión de los genes
    # de acuerdo a la columna "gene_type" agregada al marco de datos de los resultados
    # de los contrastes
    names(volcol) <- c("activado","reprimido","ns")
    
    
  
    output$distPlot <- renderPlotly({

        # Boxplots
        A <- ggplot(stack(as.data.frame(gene.expression)), aes(x = ind, y = values, fill = ind)) +
          #stat_boxplot(geom = 'errorbar') +
          geom_boxplot(outlier.shape = NA) +
          scale_y_continuous(limits = quantile(gene.expression[,1], c(0.1, 0.8))) +
          scale_fill_manual(values=c(met.brewer("Egypt",2)[1],met.brewer("Egypt",2)[1],
                                     met.brewer("Egypt",2)[2],met.brewer("Egypt",2)[2])) +
          theme(legend.position="none", axis.text.x = element_text(angle = 90),
                plot.title = element_text(hjust = 0.5, size = 13),
                panel.background = element_rect(fill = "white"),
                axis.title = element_text(size = 11),
                axis.text = element_text(size=9),
                panel.grid.major = element_line(colour = "white"),
                panel.grid.minor = element_line(colour = "white"),
                axis.line.x.bottom = element_line(color = 'black'),
                axis.line.y.left   = element_line(color = 'black'),
                panel.border = element_blank()
          ) +
          ggtitle("Datos crudos") +
          xlab("Muestras") +
          ylab("FPKM")

        B <- ggplot(stack(as.data.frame(log.gene.expression)), aes(x = ind, y = values, fill = ind)) +
          #stat_boxplot(geom = 'errorbar') +
          geom_boxplot(outlier.shape = NA) +
          scale_y_continuous(limits = quantile(log.gene.expression[,1], c(0.1, 0.8))) +
          scale_fill_manual(values=c(met.brewer("Egypt",2)[1],met.brewer("Egypt",2)[1],
                                     met.brewer("Egypt",2)[2],met.brewer("Egypt",2)[2])) +
          theme(legend.position="none", axis.text.x = element_text(angle = 90),
                plot.title = element_text(hjust = 0.5, size = 13),
                panel.background = element_rect(fill = "white"),
                axis.title = element_text(size = 11),
                axis.text = element_text(size=9),
                panel.grid.major = element_line(colour = "white"),
                panel.grid.minor = element_line(colour = "white"),
                axis.line.x.bottom = element_line(color = 'black'),
                axis.line.y.left   = element_line(color = 'black'),
                panel.border = element_blank()
          ) +
          ggtitle("Datos normalizados") +
          xlab("Muestras") +
          ylab("FPKM")

        # Plot interactivo de cada uno
        boxplot_A<-ggplotly(A)
        boxplot_B<-ggplotly(B)

        # Arreglando outliers
        boxplot_A$x$data <- lapply(boxplot_A$x$data, FUN = function(x){
          # hacerlos transparentes
          x$marker = list(opacity = 0)
          # hacerlos pequeños
          #x$marker = list(size = 0.001)
          # Con el cursor no sale la etiqueta del los outliers, pero si la de max
          x$hoverinfo= list("y"< quantile(gene.expression[,1],probs = 0.75))
          return(x)
       })


        # Arreglando outliers
        boxplot_B$x$data <- lapply(boxplot_B$x$data, FUN = function(x){
          # hacerlos transparentes
          x$marker = list(opacity = 0)
          # hacerlos pequeños
          #x$marker = list(size = 0.001)
          # Con el cursor no sale la etiqueta del los outliers, pero si la de max
          x$hoverinfo= list("y"< quantile(gene.expression[,1],probs = 0.75))
          return(x)
        })

        # Plot interactivo con los dos
        subplot(boxplot_A, boxplot_B) %>%
          layout(title = 'Datos Pre y Post-Normalizados')

    })
    output$Volcano <- renderPlotly({
      # Volcano plot estatico con ggplot2
      volcano_A<-ggplot(as.data.frame(de.results), aes(x=logFC, y=-log10(adj.P.Val),
                                                         color=gene_type,
                                                         text= paste0("</br> Gen: " ,gene_name,
                                                                      "</br> Anotacion: ",annotation))) +
        geom_point(size = 1) +
        scale_colour_manual(values = volcol) + 
        theme(legend.position = "none",
              panel.background = element_rect(fill = "white"),
              panel.grid.major = element_line(colour = "white"),
              panel.grid.minor = element_line(colour = "white"),
              axis.line.x.bottom = element_line(color = 'black'),
              axis.line.y.left   = element_line(color = 'black'),
              panel.border = element_blank(),
              plot.title = element_text(hjust = 0.5,size = 19)) + 
        ggtitle("Volcano plot")
      
      # # Volcano interactivo con plotly
      Volcano<-ggplotly(volcano_A,x = ~logFC, y = ~adj.P.Value,
               tooltip = "text",
               width = 640, height = 700,source = "Volcano")
    })
    
      
    
      # Tabla segun click en el volcano
    data <- reactive({
      de.results
    })
    output$table<- DT::renderDataTable({
      # Si no hay click, no hay na
      event.data <- event_data("plotly_click", source = "Volcano")
      print(event.data)
      if(is.null(event.data)) { return(NULL)}
      # Con click sale la  fila con informacion de ese gen
      result <- data()[data()$logFC == event.data$x ,]
      DT::datatable(result)
    })
    
        # Barplot segun click
      barplot_data<-reactive({
        resultado_final
      })

    output$Barplot<-renderPlot({
      # Si no hay click, no hay na
      event.data <- event_data("plotly_click", source = "Volcano")
      print(event.data)
      if(is.null(event.data)) { return(NULL)}

      # Seleccionamos la fila de interés segun el logFC
      result2 <- barplot_data()[barplot_data()$logFC == event.data$x ,]
      gene.id<-result2[,"geneID"]
      gene.name<-result2[,"annotation"]
      gene.expression<-result2[,c("LL_1","LL_2","HL_1","HL_2")]

      barplot.gene(gene.id,gene.name,gene.expression)

    })

 
    
    
      # Funcion original
      # output$Barplot<-renderPlot({
        # Funcion
        # barplot.gene <- function(gene.id,gene.name,gene.expression)
        # {
        #   expression.ll <- mean(unlist(gene.expression[gene.id,
        #                                                c("LL_1", "LL_2")]))
        #   expression.hl <- mean(unlist(gene.expression[gene.id,
        #                                                c("HL_1", "HL_2")]))
        #   means <- c(expression.ll,expression.hl)
        #   
        #   expression.sd.ll <- sd(unlist(gene.expression[gene.id,
        #                                                 c("LL_1", "LL_2")]))
        #   expression.sd.hl <- sd(unlist(gene.expression[gene.id,
        #                                                 c("HL_1", "HL_2")]))
        #   sds <- c(expression.sd.ll,expression.sd.hl)
        #   
        #   
        #   par(lwd=3)
        #   xpos <- barplot(means,col=c("blue","firebrick2"),
        #                   names.arg = c("LL","HL"),las=2,cex.names = 1.5,
        #                   ylim=c(0,max(means+sds)*1.2),cex.axis = 1.5,lwd=3,
        #                   main=paste(c(gene.name, "-", gene.id),collapse=" "),
        #                   cex.main=2)
        #   arrows(x0 = xpos,y0 = means + sds, x1 = xpos, y1 = means - sds,length = 0.1,
        #          code = 3,angle=90,lwd=2)
        #   return(list(means,means[2]/means[1]))
        # }
        # Plot con click
        # 
        # paste0("x=", input$plot_click$x, "\ny=", input$plot_click$y)
      }
      
    

# Run the application 
shinyApp(ui = ui, server = server)

