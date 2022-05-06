#########################################################################
## Authors: Emma Serrano Pérez                                          #
## Contact & Maintainer: Emma Serrano Pérez <emmserper@alum.us.es>      #
#########################################################################

# Necessary packages:
library(shiny) # app
library(shinythemes) # app
library(ggplot2, logical.return = T) # plots
library(MetBrewer, logical.return = T) #palette
library(cowplot, logical.return = T) # plots
library(plotly) # interactive plot
#library(ggrepel) #
library(DT) # interactive table


# Uploading data:
gene.expression_total <- read.table(file = "data/gene_expression.tsv",header=T,as.is=T)
gene.ids <- gene.expression_total$geneID
gene.expression <- as.matrix(gene.expression_total[,2:ncol(gene.expression_total)])
rownames(gene.expression) <- gene.ids

log.gene.expression <-read.table(file = "data/log_gene-expression.tsv",header=T,as.is=T)
de.results <-read.table(file = "data/de_results.tsv",header=T,as.is=T)
# resultado_final<-read.table(file = "data/resultado_final.tsv",header=T,as.is=T)
# head(resultado_final)


# Functions:
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





# Define UI
ui <- shinyUI(fluidPage(theme = shinytheme("flatly"),

  fluidRow(
    column(
      width = 2,
      img(src='logo_1.png', align = "center", width=200),
      tags$br(),
      radioButtons(inputId = "navigation_bar", width="100%",selected="home",
                   label="",
                   choices=c(
                     "Home" = "home",
                     "Relevant Information" = "intro",
                     "Global Transcriptome Statistics" = "globaltrans",
                     "Specific Gene Analysis" = "gene",
                     "Global Metabolomic Statistics" = "globalmet",
                     "Specific Metabolite Analysis" = "metabolite",
                     "Tutorials" = "tutorials",
                     "GitHub repository" = "github",
                     "Citation and Contact" = "citation"
                   ))),
    column(
      width = 8,
      tags$div(align = "center", 
               tags$h1(tags$b("Bona Nitens"), tags$br()),
               tags$h2("microALGAE Multiomic Data Exploration for MicroAlgae Transcriptomic and Metabolomic AnalysiS")),
      tags$br(),tags$br(),
      conditionalPanel(condition = "input.navigation_bar == 'home'",
                       tags$div(align = "justify", "Welcome to", tags$b("Bona Nitens"),"a web based tool for the exploration of ", 
                                tags$b("Transcriptomic"), "and ", tags$b("Metabolomic"), "data. Specially, this web tool focuses on the study of the early response to light stress
                                in the microalgae specie",tags$i("Klebsormidium nitens"), "as a trigger for the development of different protection mechanisms at the molecular and 
                                celular level. We present an integrative analysis between these transcriptomic and metabolic data to characterize
                                the molecular mechanisms involved in the response to high light in Klebsormidium and relate them to the mechanisms used by land plants
                                or Embryophyta."),
                       tags$br(),
                       tags$div(align = "justify","For more information, select from the left side navigation bar the type of exploration in which you are interested and follow the relevant intructions.
                       In addition, a", tags$b("video tutorial"), "has been provided that describes in detail the information collected in each section. The code used to generate this websites is available on", 
                       tags$b("Github."), "If you find this work interesting and useful in your research, please
                       do not hesitate to cite us."),
               tags$br(),tags$br(),
               tags$div(align ="center",img(src="optico.png", align = "center", width=600))
               # 
      ),
      conditionalPanel(condition = "input.navigation_bar == 'intro'",
                       tags$h4(tags$b(align="right","Introduction to our case study:")),
                       tags$br(),
                       tags$div(align = "justify", "The evolutionary history of the ",tags$b("green plants kingdom")," or", tags$i("Viridiplantae"),"splits into two different lineages", 
                       tags$i("Chlorophyta")," and",tags$i("Streptophyta"),".",tags$i("Chlorophyta"),"are primarily constituted by marine and freshwater green microalgae. In turn,", tags$i("Streptophyta")," are
                       divided into two different clades",tags$i("Charophyta"),"and",tags$i("Embryophyta"),". Whereas",tags$i("Embryophyta"),"comprises mainly land plants,",tags$i("Charophyta"),"are still 
                       considered algae with a preference for freshwater and with some facultative terrestrial species.  Present-day",tags$i("Charophyta")," are generally 
                       accepted as the extant algal species most closely related to the aquatic ancestors of land plants or",tags$i("Embriophyta"),". Accordingly, the molecular 
                       systems that potentially allowed this group of photosynthetic organisms to evolve towards terrestrial land plants are under intense analysis."),
                       tags$br(),
                       tags$div(align = "justify", "During this transition, the evolution of response molecular systems to",tags$b(" terrestrial environmental stresses")," was 
                                critical. Some terrestrial physiological adaptations, such as desiccation resistance and tolerance to UV radiation are present in 
                                ",tags$i("Charophyta")," from which current land plant mechanisms supposedly evolved. Other system found in",tags$i("Embryophyta"),"such as auxin transport, 
                                photoprotective capacity and adaptation to transient light changes have been identified in ",tags$i("Charophyta")," as",tags$i("Zygnema  circumcarinatum.")," Whereas
                                these studies focus mainly on genomic data, the",tags$b("lack of multi-omic data"),"such as transcriptomic and metabolomic 
                                data for  ",tags$i("Charophyta")," under specific conditions relevant to the terrestralization process is preventing the full characterization of
                                the molecular systems that promoted the transition to the first land plants. "),
                       tags$br(),tags$br(),
                       tags$div(align ="center",img(src="filo.png", align = "center", width=600), tags$br(),
                                tags$b("Figure 1."),"Phylogenetic tree of the evolution of green algae."),
                       tags$br(),tags$br(),
                       tags$div(align = "justify", "In this study, we have chosen the freshwater facultative terrestrial ",tags$i("Charophyta Klebsormidium nitens")," as",tags$b("model organism"),"
                                to study the",tags$b("transcriptomic and metabolomic response to high light intensity"),"recreating at least one of the most critical
                                environmental changes faced by plants during terrestralization. ",tags$i("K. nitens")," cultures consist of multicellular and non-branching
                                filaments without specialized cells with a single chloroplast. Many",tags$i("Klebsormidium"),"species are cosmopolitan distributed in terrestrial
                                environments as soil crusts and rocks as well as freshwater habitats like streams and rivers. "),
                       tags$br(),
                       tags$div(align = "justify", "Their presence in these environments expose cells to extreme conditions including ",tags$b("high light irradiance."),"
                                Physiological studies under such conditions have been carried out reporting photosynthetic resistance against intense light 
                                meditated by the presence of photoprotective mechanisms dissipating energy as heat (non- photochemical quenching, NPQ) and/or
                                by the activation of alternative electron routes to reduce reactive oxygen species (ROS) production. Several comparative genomic
                                analyses have been carried out providing evidence about ",tags$i("Klebsormidium")," possessing fundamental molecular mechanisms required for the
                                adaptation and survival in terrestrial environments including wax-related genes,phytohormone signaling and transcription factors 
                                involved in resistance to high light and UV radiation. Nonetheless, there are very few transcriptomic studies integrating gene 
                                expression with physiological data aiming at the characterization of ",tags$i("Klebsormidium"),"responses to abiotic stresses."),
                       tags$br(),tags$br(),
                       tags$h4(tags$b(align="right","Experimental design:")),
                       tags$br(),
                       tags$div(align = "justify",tags$i("Klebsormidium nitens")," (strain NIES-2285) was obtained from the",tags$b(" National Institute for Environmental Studies (Japan)."),"
                                Cells were grown photoautotrophically in Bold’s Basal Medium using photobioreactors containing 0.8 L of cell suspension and bubbled with
                                air supplemented with 1% (v/v) CO 2 as carbon source. Photobioreactors were continuously illuminated with white light lamps at 50 μE m
                                -2 s -1 and maintained at 20ºC. Defoamer (Antifoam 204) was added to avoid the contamination of the aeration systems. Cultures at 
                                exponential phase with 45 μg/ml chlorophyll content were used in our experiments. Control cultures were kept under a control light 
                                irradiance of 50 μE m -2 s -1 whereas high light cultures were illuminated for three hours with an irradiance of 1500 μE m -2 s -1. Six independent
                                biological replicates were considered for low and high light irradiance metabolomic data generation. Cells were collected, 
                                washed with PBS  and stored at -80ºC. "),
                       tags$br(),
                       tags$div(align = "justify", "METER VIDEO DE FRAN DE TIKTOK"),
                       tags$br(),
                       tags$div(align = "justify", "Using these cells, ",tags$b("RNA extraction"),"was performed to obtein purified RNA and the computational pipeline",tags$b("MARACAS"),"
                                was used to determine differentially expressed genes according to a log2FC of ± 1 and a adjusted p-value or FDR 
                                (False Discobery Rate) threshold of 0.05. he software tool",tags$b("AlgaeFUN"),"was used to perform functional 
                                enrichment analysis based on Gene Ontology (GO) terms and Kyoto Encyclopedia of Genes and Genomes (KEGG) pathways over the sets of differentially
                                expressed genes."),
                       tags$br(),
                       tags$div(align = "justify", "For",tags$b("metabolite content determination (primary metabolite, phytohormone and carotenoids)"),", cell pellets were lyophilized and the 
                                determination was carried out by ultra high performance liquid chromatography system coupled with mass spectrometry (UPLC/MS), or HPLC (High-performance
                                liquid chromatography) coupled to an UV-visible scanning spectrophotometer in case of carotenoids.")
                       # 
      ),
      conditionalPanel(condition = "input.navigation_bar == 'globaltrans'",
                       tags$div(align = "justify", "Todo análisis exploratorio debe comenzar con un estudio básico. En primer lugar, 
                                debe asegurarse que mimimimimi Boxplot."),
                       tags$br(),tags$br(),
                       tags$div(align = "justify", "mimimimimi PCA.")
                       # 
      ),
      conditionalPanel(condition = "input.navigation_bar == 'gene'",
                        tags$div(align="justify", tags$b("AlgaeFUN"), "allows researchers to perform", tags$b("functional annotation"), 
                       "over gene sets.", tags$b("Gene Ontology (GO) enrichment"), "analysis as well as", tags$b("KEGG (Kyoto Encyclopedia
                       of Genes and Genomes) pathway enrichment"), "analysis are supported. The gene set of interest can be obtained, for example,
                       as the result of a differential expression analysis carried out using", tags$b("MARACAS."), " See our", tags$b("video tutorial"),
                       "for details or follow the next steps to perform your analysis:",
                                tags$ol(
                                  tags$li("In the left panel choose your ", tags$b("microalgae")," of interest, the type of enrichment analysis 
                                          to perform and the", tags$b("p-value threshold.")),
                                  tags$li("Insert your ", tags$b("gene set"), " in the text box or load it from a file using the",
                                          tags$b("Browse …"), " button. An example can be loaded by clicking on the ",  tags$b("Example"), " button. 
                                          Click on", tags$b("Clear"), " button to remove the loaded gene set."),
                                  tags$li("Users can choose between the default", tags$b("background"), " gene provided by AlgaeFUN of a custom one 
                                          that can be specified."),
                                  tags$li("Click on the ", tags$b("Have Fun"), " button to perform the specified functional enrichment analysis. The
                                          results will be shown in the different tabs below.")
                                 )
      )),
      
      conditionalPanel(condition = "input.navigation_bar == 'globalmet'",
                       tags$div(align="justify", tags$b("AlgaeFUN"), "allows researchers to perform", tags$b("annotation analysis 
                                of genomic loci or regions."), "These are typically generated from", tags$b("ChIP-seq"), "studies 
                                of the genome-wide distribution of", tags$b("epigenetic marks or transcription factor binding sites."),
                                "Our tool", tags$b("MARACAS"), "can be used to perform this type of analysis. The set of marked genes 
                                can be obtained as well as the distribution of the genomic loci overlapping specific genes parts. Also, individual marked genes 
                                and the average signal level around the TSS (Transcription Start Site) and TES (Transcription
                                End Site) over the complete set of marked genes can be visualized.", " See our", tags$b("video tutorial"),
                                "for details or follow the next steps to perform your analysis:",
                                tags$ol(
                                  tags$li("In the left panel choose your ", tags$b("microalgae")," of interest, the gene", 
                                          tags$b("promoter length"), "and the",  tags$b("gene parts"), "that will be considered
                                          when determining the marked genes."),
                                  tags$li("Insert in the text box your ", tags$b("set of genomic regions"), " as a table consisting 
                                          of three tab-separated columns representing the chromosome, the start and end position of 
                                          the regions. An example can be loaded by clicking on the ",  tags$b("Example"), " button. 
                                          Click on", tags$b("Clear"), " button to remove the loaded gene set. Alternatively, using the",
                                          tags$b("Browse..."), "button, the genomic regions can be uploaded from a file in BED format as 
                                          described previously containing as least three columns. This file can be obained using our tool", 
                                          tags$b("MARACAS.")),
                                  tags$li("Optionally, users can upload the genome wide signal level of a epigenetic mark or transcription 
                                          factor binding in a BigWig file. This file can be obained using our tool", tags$b("MARACAS.")),
                                  tags$li("Click on the ", tags$b("Have Fun"), " button to perform the specified analysis. The
                                          results will be shown in the different tabs below.")
                                )
                       )),
      conditionalPanel(condition = "input.navigation_bar == 'metabolite'",
                       tags$div(align="justify", tags$b("AlgaeFUN"), "allows researchers to perform", tags$b("annotation analysis 
                                of genomic loci or regions."), "These are typically generated from", tags$b("ChIP-seq"), "studies 
                                of the genome-wide distribution of", tags$b("epigenetic marks or transcription factor binding sites."),
                                "Our tool", tags$b("MARACAS"), "can be used to perform this type of analysis. The set of marked genes 
                                can be obtained as well as the distribution of the genomic loci overlapping specific genes parts. Also, individual marked genes 
                                and the average signal level around the TSS (Transcription Start Site) and TES (Transcription
                                End Site) over the complete set of marked genes can be visualized.", " See our", tags$b("video tutorial"),
                                "for details or follow the next steps to perform your analysis:",
                                tags$ol(
                                  tags$li("In the left panel choose your ", tags$b("microalgae")," of interest, the gene", 
                                          tags$b("promoter length"), "and the",  tags$b("gene parts"), "that will be considered
                                          when determining the marked genes."),
                                  tags$li("Insert in the text box your ", tags$b("set of genomic regions"), " as a table consisting 
                                          of three tab-separated columns representing the chromosome, the start and end position of 
                                          the regions. An example can be loaded by clicking on the ",  tags$b("Example"), " button. 
                                          Click on", tags$b("Clear"), " button to remove the loaded gene set. Alternatively, using the",
                                          tags$b("Browse..."), "button, the genomic regions can be uploaded from a file in BED format as 
                                          described previously containing as least three columns. This file can be obained using our tool", 
                                          tags$b("MARACAS.")),
                                  tags$li("Optionally, users can upload the genome wide signal level of a epigenetic mark or transcription 
                                          factor binding in a BigWig file. This file can be obained using our tool", tags$b("MARACAS.")),
                                  tags$li("Click on the ", tags$b("Have Fun"), " button to perform the specified analysis. The
                                          results will be shown in the different tabs below.")
                                )
                       )),

      conditionalPanel(condition = "input.navigation_bar == 'github'",
                       tags$div(align = "justify", tags$b("AlgaeFUN,"), "is entirely developed using 
        the R package", tags$b( tags$a(href="https://shiny.rstudio.com/", "shiny.")), "The 
        source code is released under", tags$b("GNU General Public License v3.0"), "and is hosted at",
                                tags$b("GitHub."), "If you experience any problem using AlgaeFUN please create an", 
                                tags$b(tags$a(href="https://github.com/fran-romero-campero/AlgaeFUN/issues","issue")), 
                                "in GitHub and we will address it."),
                       tags$div(align="center",tags$h1(tags$b(tags$a(href="https://github.com/fran-romero-campero/AlgaeFUN",
                                                                     "AlgaeFUN at GitHub")))),
                       tags$br(),
                       tags$br(),
                       
                       tags$div(align = "justify", tags$b("MARACAS,"), "is developed using bash scripting and several 
        bioconductor R packages. The source code is released under", tags$b("GNU General Public License v3.0"), "and is hosted at",
                                tags$b("GitHub."), "If you experience any problem using AlgaeFUN please create an", 
                                tags$b(tags$a(href="https://github.com/fran-romero-campero/MARACAS/issues","issue")), 
                                "in GitHub and we will address it."),
                       tags$div(align="center",tags$h1(tags$b(tags$a(href="https://github.com/fran-romero-campero/MARACAS",
                                                                     "MARACAS at GitHub")))),
      ),
      
      conditionalPanel(condition = "input.navigation_bar == 'citation'",
                       tags$div(align = "justify", "We are strongly committed to", tags$b("open access software"), 
                                "and", tags$b("open science."),"Following our philosophy we have deposited our GitHub code 
                       into", tags$a(href="https://zenodo.org/record/4754516#.YJxLPSaxUws", target="_blank",tags$b("Zenodo")), ", a
                       general-purpose open-access repository developed under the", 
                                tags$a(href="https://www.openaire.eu/", target="_blank", tags$b("European OpenAIRE program.")), "Meanwhile we publish 
                       our work in a journal if you find", tags$b("AlgaeFUN with MARACAS"), "useful in your research we would be most grateful if you cite 
                       our GitHub repository with a,", tags$b("DOI"),  "as follows:",
                                tags$br(),
                                tags$br(),
                                tags$div(tags$b("Romero-Losada, A.B., Arvanitidou, C., de los Reyes, P., 
                                García-González, M., Romero-Campero, F.J. (2021) AlgaeFUN with MARACAS, microAlgae FUNctional 
                                enrichment tool for MicroAlgae RnA-seq and Chip-seq AnalysiS v1.0, Zenodo, doi:10.5381/zenodo.4754516 doi:10.5381/zenodo.4752818"))),
                       
                       tags$br(),
                       tags$br(),
                       #tags$div(align="center", img(src='smiley.png', align = "center", width=200,hight=200)),
                       tags$br()
                       
      ),
      
      conditionalPanel(condition = "input.navigation_bar == 'tutorials'",
                       tags$div(align="center",uiOutput("video_tutorial")),
                       tags$div(align = "justify", 
                                tags$br(),
                                tags$br(),
                                tags$div(tags$h4(tags$b("Above you can find a video tutorial on how to use the different tools implemented 
                                in AlgaeFUN with MARACAS."))))
                       
      ),
      
    ),
    column(
      width = 2,
      img(src='logo_ibvf.jpg', align = "center", width=100),
      img(src='logo_us.png', align = "center", width=100),
      tags$br(),tags$br(),tags$br(),
      img(src='logo_csic.jpg', align = "center", width=100),
      tags$br(),tags$br(),
      tags$div(align="center",width=60,
               HTML("<script type=\"text/javascript\" src=\"//rf.revolvermaps.com/0/0/8.js?i=5jamj0c2y0z&amp;m=7&amp;c=ff0000&amp;cr1=ffffff&amp;f=arial&amp;l=33\" async=\"async\"></script>"))
    )
  ),

  
  tags$br(),tags$br(),
  
  #Interface where the user can choose his/her preferencies, separated by columns
  fluidRow(
      column(width = 4,
      
     conditionalPanel(condition = "input.navigation_bar == 'basic'",
                      tags$div(align="center",uiOutput("boxplot")),
                      fileInput(inputId = "data_matrix",label = "Data Matrix:", width= "100%")
        ),
      )
  ),
))
      

server <- shinyServer(function(input, output, session) {

        
})

# Run the application 
shinyApp(ui = ui, server = server)