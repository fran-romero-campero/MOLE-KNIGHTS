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
library(tidyr) # Para el PCA


# Uploading data:
gene.expression_total <- read.table(file = "data/gene_expression.tsv",header=T,as.is=T)
gene.ids <- gene.expression_total$geneID
gene.expression <- as.matrix(gene.expression_total[,2:ncol(gene.expression_total)])
rownames(gene.expression) <- gene.ids

log.gene.expression <-read.table(file = "data/log_gene-expression.tsv",header=T,as.is=T)
de.results <-read.table(file = "data/de_results.tsv",header=T,as.is=T)
resultado_final<-read.table(file = "data/resultado_final.tsv",header=T,as.is=T)

met.results <-read.table(file = "data/met_results.tsv",header = T, as.is = T)
metabolomic.data.1<-read.table(file = "data/metabolomic_data_1",header = T, as.is = T)
metabolomic.data.2<-read.table(file = "data/metabolomic_data_2",header = T, as.is = T)
metabolites <-read.table(file = "data/metabolites.txt",as.is = T)
metabolites <- metabolites$x


# Functions:
barplot.gene <- function(gene.id,gene.name,gene.expression)
{
  expression.ll <- mean(unlist(gene.expression[,c("LL_1", "LL_2")]))
  expression.hl <- mean(unlist(gene.expression[,c("HL_1", "HL_2")]))
  means <- c(expression.ll,expression.hl)
  
  expression.sd.ll <- sd(unlist(gene.expression[,c("LL_1", "LL_2")]))
  expression.sd.hl <- sd(unlist(gene.expression[,c("HL_1", "HL_2")]))
  sds <- c(expression.sd.ll,expression.sd.hl)
  
  
  par(lwd=1.5)
  xpos <- barplot(means,col=c(met.brewer("Egypt",2)[2],met.brewer("Egypt",2)[1]),
                  names.arg = c("LL","HL"),las=2,cex.names = 1.5,
                  ylim=c(0,max(means+sds)*1.2),cex.axis = 1.5,lwd=3,
                  main=paste(c(gene.name, gene.id,sep = "\n",collapse=" ")),
                  cex.main=1.25)
  arrows(x0 = xpos,y0 = means + sds, x1 = xpos, y1 = means - sds,length = 0.1,
         code = 3,angle=90,lwd=2)
  return(list(means,means[2]/means[1]))
}


barplot_met <-function(met.name)
{
  current.metabolite <- met.name
  hl.1 <- metabolomic.data.1[1:3,current.metabolite]
  ll.1 <- metabolomic.data.1[4:6,current.metabolite]
  mean.ll.1 <- mean(ll.1)
  norm.hl.1 <- hl.1/mean.ll.1
  norm.ll.1 <- ll.1/mean.ll.1

  hl.2 <- metabolomic.data.2[1:3,current.metabolite]
  ll.2 <- metabolomic.data.2[4:6,current.metabolite]
  mean.ll.2 <- mean(ll.2)
  norm.hl.2 <- hl.2/mean.ll.2
  norm.ll.2 <- ll.2/mean.ll.2
  

  norm.ll <- c(norm.ll.1, norm.ll.2)
  norm.hl <- c(norm.hl.1, norm.hl.2)
  means <- c(mean(norm.ll), mean(norm.hl))
  sds <- c(sd(norm.ll), sd(norm.hl))
  
  xpos <- barplot(means,col=c(met.brewer("Egypt",2)[2],met.brewer("Egypt",2)[1]),
                  names.arg = c("LL","HL"),las=2,cex.names = 1.5,
                  ylim=c(0,max(means+sds)*1.3),cex.axis = 1.5,lwd=3,
                  main=current.metabolite,
                  cex.main=2)
  arrows(x0 = xpos,y0 = means + sds, x1 = xpos, y1 = means - sds,length = 0.1,
         code = 3,angle=90,lwd=2)
}


# Define UI
ui <- shinyUI(fluidPage(theme = shinytheme("flatly"),

  fluidRow(
    column(
      width = 2,
      img(src='confocal.png', align = "center", width=150),
      tags$br(),
      radioButtons(inputId = "navigation_bar", width="100%",selected="home",
                   label="",
                   choices=c(
                     "Sweet Home Alabama!" = "home",
                     "Previous Research Information" = "intro",
                     "Global Transcriptome Statistics" = "globaltrans",
                     "Specific Gene Analysis" = "gene",
                     "Global Metabolomic Statistics" = "globalmet",
                     "Specific Metabolite Analysis" = "metabolite",
                     "Omics... Assemble!" = "integration",
                     "Tutorials" = "tutorials",
                     "GitHub repository" = "github",
                     "Citation and Contact" = "citation"
                   ))),
    column(
      width = 8,
      tags$div(align = "center", 
               tags$h1(tags$b("Zzz... BONA Nitens!"), tags$br()),
               tags$h2("Best multiOmic exploratioN of trAnscriptomic
                       and metabolomic data in klebsormidium Nitens")),
      tags$br(),tags$br(),
      # Home:
      conditionalPanel(condition = "input.navigation_bar == 'home'",
                       tags$div(align = "justify", "Welcome to", tags$b("BONA Nitens"),"a web based tool for the exploration of ", 
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
      # Relevant Information:
      conditionalPanel(condition = "input.navigation_bar == 'intro'",
                       tags$h4(tags$b(align="right","Introduction to our case study:")),
                       tags$br(),
                       tags$div(align = "justify", "The evolutionary history of the ",tags$b("green plants kingdom")," or", tags$i("Viridiplantae"),"splits into two different lineages", 
                       tags$i("Chlorophyta")," and",tags$i("Streptophyta"),".",tags$i("Chlorophyta"),"are primarily constituted by marine and freshwater green microalgae. In turn,", tags$i("Streptophyta")," are
                       divided into two different clades",tags$i("Charophyta"),"and",tags$i("Embryophyta"),". Whereas",tags$i("Embryophyta"),"comprises mainly land plants,",tags$i("Charophyta"),"are still 
                       considered algae with a preference for freshwater and with some facultative terrestrial species.  Present-day",tags$i("Charophyta")," are generally 
                       accepted as the extant algal species most closely related to the aquatic ancestors of land plants or",tags$i("Embriophyta"),tags$b("(Fig. 1)"),". Accordingly, the molecular 
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
                       tags$div(align ="center",img(src="phylogeny.png", align = "center", width=800), tags$br(),
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
                       tags$div(align = "center", tags$video(src="experiment.mp4", type = "video/mp4",height ="400px", width="400px",controls="controls")),
                       tags$br(),
                       tags$div(align = "justify", "Using these cells, ",tags$b("RNA extraction"),"was performed to obtein purified RNA and the computational pipeline",
                                tags$a(href="https://github.com/fran-romero-campero/AlgaeFUN","MARACAS"),
                                "was used to determine differentially expressed genes according to a log2FC of ± 1 and a adjusted p-value or FDR 
                                (False Discobery Rate) threshold of 0.05. he software tool",tags$a(href="https://github.com/fran-romero-campero/AlgaeFUN","AlgaeFun"),"was used to perform functional 
                                enrichment analysis based on Gene Ontology (GO) terms and Kyoto Encyclopedia of Genes and Genomes (KEGG) pathways over the sets of differentially
                                expressed genes."),
                       tags$br(),
                       tags$div(align = "justify", "For",tags$b("metabolite content determination (primary metabolite, phytohormone and carotenoids)"),", cell pellets were lyophilized and the 
                                determination was carried out by ultra high performance liquid chromatography system coupled with mass spectrometry (UPLC/MS), or HPLC (High-performance
                                liquid chromatography) coupled to an UV-visible scanning spectrophotometer in case of carotenoids.")
                       # 
      ),
      # Global Transcriptomic Statistics:
      conditionalPanel(condition = "input.navigation_bar == 'globaltrans'",
                       tags$div(align = "justify", "One of the main ways of responding to changes in the environment is by",tags$b("modifying gene expression"),"by varying the number of transcripts
                                present in cells that conform a specific transcriptome. For this reason, the",tags$b("RNA-seq technique"),"was developed, a massive cDNA sequencing technique obtained through RNA 
                                extraction using high-performance sequencing. In this case, the study focused on RNA-seq applied on eukaryotic coding RNA."),
                       tags$br(),
                       
                       tags$h4(tags$b(align="right","Previous Mathematical/Computational Analysis:")),
                       tags$br(),
                       tags$div(align = "justify", "The previous mathematical/computational analysis was carried out according to the protocol established in the",tags$a(href="https://github.com/fran-romero-campero/AlgaeFUN","MARACAS")," specialized
                                in microalgae using the",tags$b("HISAT2 mapper"),", obtaining as a result a",tags$b("normalized continuous data matrix by FPKM"),"(Fragments Per Kilobase of
                                transcripts of mapped reads) that normalizes by the length of the read and by the number of millions of reads. In addition, it applies a second upper quantile
                                normalization and log2 transformation that, as can be seen",tags$b("(Fig. 2)"),", improves the distribution of the data."),
                       tags$br(),tags$br(),
                       tags$div(align="center",plotlyOutput("distPlot"), tags$br(),
                                tags$b("Figure 2."),"Boxplot before and after normalization. Please, hover over the plot for more information."),
                       tags$br(),tags$br(),
                       tags$div(align="justify","A quick way to explain the variability between samples is using a",tags$b("PCA (Principal Component Analysis)")," where the 3 main components are represented. In this case,
                                it is observed that the first 2 components would be sufficient since they collect most of the variability of the data."),
                       tags$br(),
                       tags$div(align="center",plotlyOutput("PCA_trans")),
                       
                       tags$h4(tags$b(align="right","Differentially Expressed Genes:")),
                       tags$br(),
                       tags$div(align="justify", "Moreover, this software allows to determinate",tags$b("differentially expressed genes"),". In this case, expression was detected
                       in in",tags$b("68.4 %"),"of the 17290 genes in
                       the current Klebsormidium genome annotation. According to a log2FC of ± 1 and a q-value or FDR (False Discobery Rate) threshold of 0.05,"
                       ,tags$b("significant differential gene expression"),"
                       was determined after three hours of high light treatment in",tags$b("7.84 %"),"of the entire Klebsormidium genome with respect to low light conditions. 
                       Specifically, we identified",tags$b("667 activates")," 
                       and",tags$b(" 678 repressed genes (Fig. 3).")),
                       tags$br(),tags$br(),
                       tags$div(splitLayout(cellWidths = c("50%", "50%"), plotlyOutput("Volcano"), plotOutput("Barplot")), tags$br(),align="center",
                                tags$b("Figure 3."),"Volcano plot of DEGs. Please, hover over the plot for more information and select any point of your interest."),
                       tags$br(),
                       DT::dataTableOutput("table"),
                       tags$br(),tags$br(),
                       tags$h4(tags$b(align="right","Functional Enrichment Analysis:")),
                       tags$br(),
                       tags$div(align="justify","The software tool",tags$a(href="https://github.com/fran-romero-campero/AlgaeFUN","AlgaeFUN")," was used to perform",tags$b("functional enrichment analysis"),"
                                based on",tags$b("Gene Ontology (GO)"),"terms, to identify the cellular components and biological processes significantly affected by high light,
                                and",tags$b("Kyoto Encyclopedia of Genes and Genomes (KEGG)"),"pathways over the sets of differentially expressed genes",tags$b("Fig. 4.")),
                       tags$br(),
                       tags$div(align="justify","The proteins encoded by differentially expressed genes,",tags$b("both activated and repressed genes"),", were significantly localized
                                in the",tags$b("chloroplast thylakoid membranes"),"indicating the initiation of a major chloroplast reprogramming. "),
                       tags$br(),
                       tags$div(align ="center",img(src="figure_4.png", align = "center", width=600), tags$br(),
                                tags$b("Figure 4.")," GO terms enrichment."),
                       tags$br(),
                       tags$div(align="justify","Specifically, proteins encoded by",tags$b("repressed genes"),"were significantly associated with",tags$b("photosystems"),"and
                                cellular structures present during",tags$b("cell division"),"such as condensed nuclear chromosomes and microtubules ",tags$b("Fig. 4."),". Accordingly,
                                photosynthesis, hexose biosynthesis, cell cycle and DNA metabolism were significantly enriched processes in the repressed
                                genes. This points to an arrest in the photosynthetic machinery and cell cycle progression as response to high light. "),
                       tags$br(),
                       tags$div(align="justify", "Proteins encoded by",tags$b("activated genes"),"are, in turn, significantly localized in cellular structures involved
                                in ",tags$b("de novo protein biosybthesis")," such as preribosomes and translation initiation factor 3’ complex ",tags$b("Fig. 4.")," In particular, categories
                                encompassing ribosome biogenesis, cytoplasmic translation initiation and protein folding were significantly enriched in the
                                activated genes. Moreover, response to oxidative stress, response to high light intensity, tetraterpenoid and carotenoid
                                metabolism were identified as significantly activated processes. This suggests an",tags$b("activation of repair and protective mechanisms")," 
                                to damages caused by high light.")
                       # 
      ),
      # Specific Gen Study:
      conditionalPanel(condition = "input.navigation_bar == 'gene'",
                        tags$div(align="justify", "In the",tags$b("Global Transcriptomic Statistics")," tab, a global exploration of the transcriptomic data is performed, but it may be interesting
                                 to study a single gene due to possible implications in different aspects. To do this, enter the",tags$b("gene identifier")," in the 
                                 lower box and you will be able to see a barplot with the difference in gene expression between conditions."),
                       tags$br(),tags$br(),
                       textInput("gen", label="Please, choose a gen:", placeholder = ""),
                       actionButton("do","Update Now!",icon ("sync")),
                       tags$br(),tags$br(),
                       plotOutput("barplot_selection")
                       
                       
                      
      ),
      # Global Metabolomic Statistics:
      conditionalPanel(condition = "input.navigation_bar == 'globalmet'",
                       tags$div(align="justify","",tags$b("Metabolomics"),"is understood as the discipline that analyzes the metabolites of a living organism and tries to find the interaction between
                                metabolic pathways as well as quantify the largest amount of metabolites present. In the omics context, it offers an additional view on the
                                characteristics of a metabolism, giving an idea of points of regulation and association of functions for unknown genes. In this case, a",tags$b("metabolic profile"),"
                                was performed where interrelated compounds were studied, generating a specific",tags$b("metabolome")," for that sample using",tags$b("mass spectrometry.")),
                       tags$br(),
                       tags$h4(tags$b(align="right","Previous Mathematical/Computational Analysis:")),
                       tags$br(),
                       tags$div(align="justify","Analogously to the Global Transcriptomic Statistics section, the metabolic data were normalized. Specifically, from the data obtained by the
                                mass spectrometer, the amount of metabolite present was determined according to the area under the curve of the peak obtained and assigned to each metabolite.
                                They were relativized based on the total weight of the sample and the presence of a standard compound: paracetamol."),
                       tags$br(),tags$br(),
                       tags$h4(tags$b(align="right","Differentially Expressed Metabolites:")),
                       tags$br(),
                       tags$div(align="justify","Six independent biological replicates were
                                considered for both, high and low light conditions. We detected 69 different primary and secondary metabolites including most amino acids and some phytohormones.
                                Significant differentially abundant metabolites were identified by performing the non-parametric Wilcoxon test using a p-value threshold
                                of 0.05. We found 12 significantly more abundant and 8 less abundant metabolites under high light when compared to low light (Figure 5). For instance, under
                                high light, we detected significant changes in specific carotenoids, accumulation of the amino acid tryptophan and the phytohormone indole-3- acetic acid (IAA)"),
                       tags$br(),tags$br(),
                       tags$div(splitLayout(cellWidths = c("50%", "50%"), plotlyOutput("Volcano_met"), plotOutput("Barplot_met")), tags$br(),align="center",
                                tags$b("Figure 5."),"Volcano plot of Metabolites. Please, hover over the plot for more information and select any point of your interest."),
                       tags$br(),
                       DT::dataTableOutput("table_met"),
                       tags$br()
                       
                       ),
      
      # Specific Metabolite Study:
      conditionalPanel(condition = "input.navigation_bar == 'metabolite'",
                       tags$div(align="justify", "In the",tags$b("Global Metabolomic Statistics")," tab, a global exploration of the metabolomic data is performed, but it may be interesting
                                 to study a single metabolite due to possible implications in different aspects. To do this, enter the",tags$b("metabolite name")," in the 
                                 lower box and you will be able to see a barplot with the difference in abundance between conditions."),
                       tags$br(),tags$br(),
                       selectInput("metabolite",label="Please, choose a metabolite:",choices = metabolites,selected = "Arginine"),
                       actionButton("do_met","Update Now!",icon ("sync")),
                       tags$br(),tags$br(),
                       plotOutput("barplot_selection_metabolites")
                       ),
      
      # Total integration
      conditionalPanel(condition = "input.navigation_bar == 'integration'",
                       tags$div(align="justify","Until now we have explored information from different omics separately, but the
                                real insight lies in the integration of both. Bioinformatics offers the possibility of obtaining
                                a more complete view of the behavior of biological systems, as will be seen in the following results.
                                If you are interested in more details, you can get it in our complete article."),
                       tags$br(),tags$br(),
                       tags$h4(tags$b(align="right","An activation of the carotenoid biosynthesis β-branch and xantophyll cycle is observed:")),
                       tags$br(),
                       tags$div(align="justify", "Here, we present an integrated 
                       transcriptomic and metabolomic analysis of this specific photoprotective response to high light in 
                       Klebsormidium",tags$b("(Figure 6)"),". As can be seen, carotenoid biosynthesis is favored for the generation of beta-carotenoid
                       by the activation of several enzymes of the pathway, while the ε-branch leading to lutein is repressed. Although, β–carotene
                       content was similar under low and high light conditions, the cycle of xanthophylls towards the formation of zeaxanthin from
                       violaxanthin is favored. Violaxanthin content decreased 4.73 fold whereas antheraxanthin and zeaxanthin contents were 
                       increased 3.44 and 41.5 fold respectively under high light when compared to low light. Accordingly, 
                       the gene encoding the enzyme involved in the xanthophyll cycle, violaxanthin de-epoxidase (VDE, 
                       kfl00604_0070) converting violaxanthin into antheraxanthin and zeaxanthin was activated 1.86 fold.
                       Furthermore, the gene encoding zeaxanthin epoxidase (ZEP, kfl00092_0060) that catalyzes the
                       synthesis of violaxanthin from zeaxanthin and antheraxanthin was 3.84 fold repressed under high 
                       light."),
                       tags$br(),
                       tags$div(align ="center",img(src="figure6.png", align = "center", width=600), tags$br(),
                                tags$b("Figure 6.")," Carotenoids."),
                       tags$br(),
                       tags$div(align="justify","In the xanthophyll cycle, the interconversion of violaxanthin into antheraxanthin and 
                       zeaxanthin, constitutes one of the major photoprotective mechanism in Embryophyta and 
                       Chlorophyta. High light induces the mobilization of violaxanthin to zeaxanthin whereas low light or 
                       darkness produce the reverse reaction (Goss and Jakob, 2010; Latowski et al., 2011). De-epoxidation 
                       of violaxanthin to zeaxanthin enhances dissipation of excess excitation energy (non-photochemical 
                       quenching, NPQ) in the photosystem II (PSII) antenna, thereby preventing inactivation and damage 
                       to the photosynthetic apparatus. NPQ is considered a fundamental mechanism for Streptophyta 
                       adaptation to terrestrial habitats (Pierangelini et al., 2017). Here, we specifically show that the 
                       xanthophyll cycle is part of the early transcriptomic and metabolomic response to high light intensity 
                       in the Charophyta Klebsormidium."),
                       tags$br(),tags$br(),
                       tags$h4(tags$b(align="right","Chloroplast retrograde signaling triggered by oxidative stress and protein misfolding is 
                                      identified as a response to high light:")),
                       tags$br(),
                       tags$div(align="justify","Under high light conditions exceeding photosynthetic capacity, production of harmful reactive 
                       oxygen species (ROS) is unavoidable associated with electron transport in the photosystems. Excess 
                       electron leakage to molecular oxygen and incomplete water oxidation produce singlet oxygen, 
                       superoxide, hydrogen peroxide and hydroxyl radical (Pospíšil, 2016). This 
                       triggers a signaling cascade communicating the chloroplast state to the nucleus termed retrograde 
                       signaling that ultimately induces the expression of nuclear genes. The evolution of this system has played a central
                       role in plant terrestralization (Zhao et al., 2019; Calderon and Strand, 2021)."),
                       tags$br(),
                       tags$div(align="justify","
                       Indeed, response to oxidative stress 
                       was one of the most significant GO term in our functional enrichment analysis over the activated 
                       genes in a response to high light treatment in Klebsormidium", tags$b("Figure 7.")," Under 
                       these conditions proteins suffer oxidative damage specifically but not limited to the active thiol 
                       groups of cysteine residues, which are oxidized to disulfide bonds (Cejudo et al., 2021). This 
                       produces major modifications in protein structure that can lead to misfolding and loss of function. 
                       The accumulation in the chloroplast of aberrant misfolded proteins also contributes to initiate 
                       retrograde signaling (Dogra et al., 2019a). Moreover, we found the activation of multiple chloroplast
                       targeted chaperones, co-chaperones and chaperonins that would contribute to restore misfolded proteins."),
                       tags$br(),
                       tags$div(align="center",img(src="figure7.png", align = "center", width=600), tags$br(),
                                tags$b("Figure 7.")," Chloroplasts."),
                       tags$br(),
                       tags$div(align="justify","Concomitant to the activation of protein repair mechanisms we found significant activation of 
                       ribosome biogenesis and cytoplasmic translation initiation. These strongly activated processes are required for de novo 
                       protein synthesis and, together with the previously described protein repair mechanisms, constitute 
                       part of the response to high light in Klebsormidium, contributing to maintain proteome homeostasis under this stress.
                       Besides, the retrograde signaling pathways induced by ROS and aberrant misfolded proteins discussed above, there exists
                       another pathway regulated by the accumulation of 3′- phosphoadenosine-5′-phosphate (PAP). The inositol polyphosphate 
                       1-phosphatase SAL1 removes PAP preventing its accumulation. The gene encoding this enzyme kfl00096_0240 was 2 fold 
                       repressed indicating a possible accumulation of PAP and an activation of the SAL1-PAP retrograde signaling pathway,
                       as a response to high light intensity in Klebsormidium.")
      
                       
                       
                       
        
      ),
      
      # Github
      conditionalPanel(condition = "input.navigation_bar == 'github'",
                       tags$div(align="justify", tags$b("Bona Nitens,"), "is entirely developed using 
        the R package", tags$b( tags$a(href="https://shiny.rstudio.com/", "shiny.")), "The 
        source code is released under", tags$b("GNU General Public License v3.0"), "and is hosted at",
                                tags$b("GitHub."), "If you experience any problem using Bona Nitens please create an", 
                                tags$b(tags$a(href="https://github.com/fran-romero-campero/AlgaeFUN/issues","issue")), 
                                "in GitHub and we will address it."),
                       tags$div(align="center",tags$h1(tags$b(tags$a(href="https://github.com/fran-romero-campero/AlgaeFUN",
                                                                     "Bona Nitens at GitHub")))),
                       tags$br(),
                       tags$br(),
                       tags$div(align="justify","Here we present some of the other main programs that we have used in
                                this exploratory tool:"),
                       tags$br(),
                       tags$div(align = "justify", tags$b("AlgaeFUN,"), "is entirely developed using 
        the R package", tags$b( tags$a(href="https://shiny.rstudio.com/", "shiny.")), "The 
        source code is released under", tags$b("GNU General Public License v3.0"), "and is hosted at",
                                tags$b("GitHub.")),
                       tags$div(align="center",tags$h4(tags$b(tags$a(href="https://github.com/fran-romero-campero/AlgaeFUN",
                                                                     "AlgaeFUN at GitHub")))),
                       tags$br(),
                       
                       tags$div(align = "justify", tags$b("MARACAS,"), "is developed using bash scripting and several 
        bioconductor R packages. The source code is released under", tags$b("GNU General Public License v3.0"), "and is hosted at",
                                tags$b("GitHub.")), 
                       tags$div(align="center",tags$h4(tags$b(tags$a(href="https://github.com/fran-romero-campero/MARACAS",
                                                                     "MARACAS at GitHub"))))
      ),
      # Citations
      conditionalPanel(condition = "input.navigation_bar == 'citation'",
                       tags$div(align = "justify","Recently we published
                       our work in a journal if you find", tags$b("this information"), "useful in your research we would be most grateful if you cite 
                       our GitHub repository with a,", tags$b("DOI"), "as follows:",
                                tags$br(),
                                tags$br(),
                                tags$div(tags$b("Serrano-Pérez E, Romero-Losada AB, Morales-Pineda M, García-Gómez ME, Couso I, García-González M
                                and Romero-Campero FJ (2022) Transcriptomic and Metabolomic Response to High Light in the Charophyte Alga Klebsormidium nitens.
                                                Front. Plant Sci. 13:855243. doi: 10.3389/fpls.2022.855243"))),
                       
                       tags$br(),
                       tags$br(),
                       tags$br()
                       
      ),
      # Tutorial
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
      )
  ),

  
  tags$br(),tags$br(),
  
  
))
      

server <- shinyServer(function(input, output, session) {
  
# Boxplots before and after
  output$distPlot <- renderPlotly({
    
    # Boxplots
    A <- ggplot(stack(as.data.frame(gene.expression)), aes(x = ind, y = values, fill = ind)) +
      #stat_boxplot(geom = 'errorbar') +
      geom_boxplot(outlier.shape = NA) +
      scale_y_continuous(limits = quantile(gene.expression[,1], c(0.1, 0.8))) +
      scale_fill_manual(values=c(met.brewer("Egypt",2)[2],met.brewer("Egypt",2)[2],
                                 met.brewer("Egypt",2)[1],met.brewer("Egypt",2)[1])) +
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
      # ggtitle("Datos crudos") +
      xlab("Samples before normalization") +
      ylab("FPKM")
    
    B <- ggplot(stack(as.data.frame(log.gene.expression)), aes(x = ind, y = values, fill = ind)) +
      #stat_boxplot(geom = 'errorbar') +
      geom_boxplot(outlier.shape = NA) +
      scale_y_continuous(limits = quantile(log.gene.expression[,1], c(0.1, 0.8))) +
      scale_fill_manual(values=c(met.brewer("Egypt",2)[2],met.brewer("Egypt",2)[2],
                                 met.brewer("Egypt",2)[1],met.brewer("Egypt",2)[1])) +
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
      # ggtitle("Datos normalizados") +
      xlab("Samples after normalization") +
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
    subplot(boxplot_A, boxplot_B,titleX = T,titleY = T) # %>%
      # layout(title = 'Before and After Normalization') %>%  
  })
  
# PCA trans
  output$PCA_trans<-renderPlotly({
    pca_all_samples <- prcomp(x=t(log.gene.expression),center = T,scale. = F)
    PCs <- as.data.frame(pca_all_samples$x)
    
    x_axis <- 1
    y_axis <- 2
    z_axis <- 3
    
    targets <- colnames(log.gene.expression) # samples
    targets <- data.frame(targets) # samples como df
    targets <- separate(targets, col = targets, sep = '_', into = c('condition','replicate')) # Separamos la condicion de la replica
    targets$sample <- colnames(log.gene.expression) # columna con la muestra completa
    
    color_variable <- targets$condition
    id <- targets$sample
    
    xlab <- paste0('PC_',x_axis,'(', round((pca_all_samples$sdev^2)[x_axis]/sum(pca_all_samples$sdev^2)*100, 2), '%)')
    ylab <- paste0('PC_',y_axis,'(', round((pca_all_samples$sdev^2)[y_axis]/sum(pca_all_samples$sdev^2)*100, 2), '%)')
    zlab <- paste0('PC_',z_axis,'(', round((pca_all_samples$sdev^2)[z_axis]/sum(pca_all_samples$sdev^2)*100, 2), '%)')
    
    trace <- list()
    i=1
    
    for (x in unique(color_variable)){
      
      trace[[i]] <- list(mode="markers", 
                         name = x,
                         type = "scatter3d",
                         x = PCs[,1][color_variable==x],
                         y = PCs[,2][color_variable==x],
                         z = PCs[,3][color_variable==x],
                         text = id[color_variable==x]
      )
      
      i <- i + 1
    }
    
    layout <- list(
      scene = list(
        xaxis = list(
          title = xlab, 
          showline = FALSE
        ), 
        yaxis = list(
          title = ylab, 
          showline = FALSE
        ), 
        zaxis = list(
          title = zlab, 
          showline = FALSE
        )
      ), 
      title = "PCA Transcriptomic Data"
    )
    
    # Creamos el plot interactivo vacío
    p <- plot_ly()
    for (x in trace){
      p <- add_trace(p, mode=x$mode, name=x$name, type=x$type, x=x$x, y=x$y, z=x$z, text=x$text)
    }
    p <- layout(p, scene=layout$scene, title=layout$title)
    
    
  })
  
  
  
  
# Volcano plot
  volcol <- c(met.brewer("Egypt",3)[1], met.brewer("Egypt",3)[2],"grey33")
  
  # Se asigna a cada uno de los colores los posibles niveles de expresión de los genes
  # de acuerdo a la columna "gene_type" agregada al marco de datos de los resultados
  # de los contrastes
  names(volcol) <- c("activado","reprimido","ns")
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
            plot.title = element_text(hjust = 0.5,size = 19))# + 
      #ggtitle("Volcano plot")
    
    # # Volcano interactivo con plotly
    Volcano<-ggplotly(volcano_A,x = ~logFC, y = ~adj.P.Value,
                      tooltip = "text",
                      width = 350, height = 400,source = "Volcano")
  })

  # Tabla que salga con el click
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
    DT::datatable(result, options = list(searching = FALSE, scrollX = T,lengthChange = FALSE))%>%
      
      formatStyle( 0, target= 'row',color = 'black', lineHeight='80%')
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
  
  # Seleccion de barplots
  gen_name<-eventReactive( input$do ,{
    input$gen
  })
  output$barplot_selection<-renderPlot({
    
    # Seleccionamos la fila de interés segun el logFC
    # kfl00093_0070
    # kfl00458_0030
    # kfl00361_0130
    if(gen_name()!="")
    {
      result2<-resultado_final[resultado_final$gene_name== gen_name(),]
      
      gene.id<-result2[,"geneID"]
      gene.name<-result2[,"annotation"]
      gene.expression<-result2[,c("LL_1","LL_2","HL_1","HL_2")]
      
      barplot.gene(gene.id,gene.name,gene.expression)
    }
    
  })
  
# Volcano Metabolite
  volcol <- c(met.brewer("Egypt",3)[1], met.brewer("Egypt",3)[2],"grey33")
  names(volcol) <- c("activado","reprimido","ns")
  
  output$Volcano_met <- renderPlotly({
    # Volcano plot estatico con ggplot2
    volcano_met <- ggplot(as.data.frame(met.results), aes(x=metabolites.fold.change, y=-log10(metabolites.p.value),
                                                          color=met_type,
                                                          text= paste0("</br> Metabolilte: " ,met_name))) +
      geom_point(size = 1.2) +
      scale_colour_manual(values = volcol) + 
      theme(legend.position = "none",
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            axis.line.x.bottom = element_line(color = 'black'),
            axis.line.y.left   = element_line(color = 'black'),
            panel.border = element_blank(),
            plot.title = element_text(hjust = 0.5,size = 15)) + 
      ggtitle("Volcano plot from Metabolites")+ 
      labs(x = "logFC",y="adj.P.Value") +
      xlim(-1, 3) +
      geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "grey33",size=0.3)
    
    
    # # Volcano interactivo con plotly
    Volcano_met <- ggplotly(volcano_met,x = ~logFC, y = ~adj.P.Value,
                                    tooltip = "text",
                                    width = 400, height = 300, source = "Volcano_met")
  })
  
  # Tabla que salga con el click
  data_met <- reactive({
    met.results
  })
  output$table_met<- DT::renderDataTable({
    # Si no hay click, no hay na
    event.data_met <- event_data("plotly_click", source = "Volcano_met")
    print(event.data_met)
    if(is.null(event.data_met)) { return(NULL)}
    # Con click sale la  fila con informacion de ese gen
    result_met <- data_met()[data_met()$metabolites.fold.change == event.data_met$x ,]
    DT::datatable(result_met, options = list(searching = FALSE, scrollX = T,lengthChange = FALSE))%>%

      formatStyle( 0, target= 'row',color = 'black', lineHeight='80%')
  })
  
  # # Barplot metabolitos segun click
  barplot_data_met<-reactive({
    met.results
  })

  output$Barplot_met<-renderPlot({
    # Si no hay click, no hay na
    event.data_met <- event_data("plotly_click", source = "Volcano_met")
    print(event.data_met)
    if(is.null(event.data_met)) { return(NULL)}

    # Seleccionamos la fila de interés segun el logFC de los metabolitos
    result_met <- barplot_data_met()[barplot_data_met()$metabolites.fold.change == event.data_met$x ,]
    met_name<-result_met[,"met_name"]
    barplot_met(met_name)

  })
  
  # Seleccion de barplots metabolites
  met_name_bar<-eventReactive( input$do_met ,{
    input$metabolite
  })
  
  output$barplot_selection_metabolites<-renderPlot({
      barplot_met(met_name_bar())
  })
  
   
  
# parentesis y corchete del server          
})

# Run the application 
shinyApp(ui = ui, server = server)