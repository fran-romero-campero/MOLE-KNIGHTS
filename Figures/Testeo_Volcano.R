

# Cargamos los datos y asjutamos los datos
gene.expression <- read.table(file = "gene_expression.tsv",header=T,as.is=T)
head(gene.expression)
gene.ids <- gene.expression$geneID
gene.expression <- as.matrix(gene.expression[,2:ncol(gene.expression)])
rownames(gene.expression) <- gene.ids
head(gene.expression)


# Para las distintas representaciones de los datos, se utilizan funciones del
# paquete ggplot2, y una de las paletas de colores interpretables por personas
# con daltonismo ofrecidas en el paquete MetBrewer
library(ggplot2, logical.return = T)
#install.packages("MetBrewer")
library(MetBrewer, logical.return = T)
#install.packages("cowplot")
library(cowplot, logical.return = T)
library(plotly)

## Aplicamos la normalización:
upper.quantiles <- vector(mode="numeric",length=ncol(gene.expression))

for(i in 1:ncol(gene.expression))
{
  upper.quantiles[i] <- quantile(gene.expression[,i],probs=0.75)
}

mean.upper.quantiles <- mean(upper.quantiles)

for(i in 1:ncol(gene.expression))
{
  gene.expression[,i] <- (gene.expression[,i] / upper.quantiles[i]) * mean.upper.quantiles
}

## Log2 transformation:
log.gene.expression <- log2(gene.expression+1)




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

# Con la función plot_grid, del paquete cowplot, se pueden representar juntos
# diferentes plots generados con ggplot
plot_grid(A, B, labels=c("A", "B"), ncol = 2, nrow = 1)

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

boxplot_A

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

boxplot_B

# Plot interactivo con los dos
fig <- subplot(boxplot_A, boxplot_B) %>% 
  layout(title = 'Datos Pre y Post-Normalizados') 
fig





# Representación datos sin normalizar sin ggplot ni plotly:
par(mfrow=c(1,2))
boxplot(gene.expression,col=rainbow(ncol(gene.expression)),ylab="FPKM"
        ,cex.lab=1.5,outline=F,las=2,main="NO-Normalized Gene Expression")

# Boxplot de la normalizacion sin ggplot ni plotly
boxplot(log.gene.expression,col=rainbow(ncol(gene.expression)),ylab="log2(FPKM + 1)",
        cex.lab=1.5,outline=F,las=2,main="Normalized Gene Expression")





# Análisis con limma:
library(limma)

## Specification of the experimental design
limma.experimental.design <- model.matrix(~ -1+factor(c(1,1,2,2)))
colnames(limma.experimental.design) <- c("LL","HL")

## Linear model fit
linear.fit <- lmFit(log.gene.expression, limma.experimental.design)

## Contrast specification and computation
contrast.matrix <- makeContrasts(contrasts = "HL-LL",
                                 levels=c("LL","HL"))

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)

## Extract results
de.results <- topTable(contrast.results, number=nrow(gene.expression),coef=1,sort.by="logFC")

# fold.change <- de.results$logFC
# q.values <- de.results$adj.P.Val
# genes.ids <- rownames(de.results)
# 
# names(fold.change) <- genes.ids
# names(q.values) <- genes.ids
# 
# activated.genes <- genes.ids[fold.change > log2(2) & q.values < 0.05]
# repressed.genes <- genes.ids[fold.change < - log2(2) & q.values < 0.05]
# 
# length(activated.genes)
# length(repressed.genes)
# 
# log10.qval <- -log10(q.values)

# Volcano de MARACAS:
##################################################################
# plot(fold.change,log10.qval,pch=19,cex=0.7,col="grey", xlab="Fold Change", ylab="-log10(q-value)",cex.lab=1.5)
# points(fold.change[activated.genes],log10.qval[activated.genes],cex=0.7,col="red",pch=19)
# points(fold.change[repressed.genes],log10.qval[repressed.genes],cex=0.7,col="blue",pch=19)
##################################################################

# Volcano ggplot
library(ggrepel) # Para añadir la etiqueta sin que solape con el punto correspondiente

# Se añade una nueva variable de acuerdo al nivel de expresión del gen para 
# la posterior representación en volcano plots.
de.results[,"gene_type"] <- "ns" # Se rellena toda la columna con: No significativo ni sustancial
head(de.results)
de.results[,"gene_name"] <- rownames(de.results)

# Se modifica ns por activado o reprimido:
de.results[which(de.results$logFC > 2 & de.results$adj.P.Val < 0.05),"gene_type"] <- "activado"
de.results[which(de.results$logFC< -2 & de.results$adj.P.Val < 0.05),"gene_type"] <- "reprimido"

genes.ids.de.results <- rownames(de.results)

# Añadimos la anotacion:
anotacion<-read.table("131203_gene_list.txt",header = T,sep = "\t")
head(anotacion)
de.results[,"annotation"]<-"hypothetical protein"
head(de.results,20)
nrow(de.results)

# Metemos la anotación en la tabla de de.results
for (i in 1:nrow(anotacion))
{
  coincidencia<-subset(anotacion,anotacion[,1] == de.results[i,"gene_name"])
  if(is.na(coincidencia[1,3])== T)
  {
    de.results[i,"annotation"]<-"hypothetical protein"
  } else {
    de.results[i,"annotation"]<-coincidencia[1,3]
  }
}

### CODIGO SUCIO ###########################################################
coincidencia<-subset(anotacion,anotacion[,1] == de.results[1,"gene_name"])
if(is.na(coincidencia[1,3])== T)
{
  de.results[1,"annotation"]<-"hypo"
} else {
  de.results[1,"annotation"]<-coincidencia[1,3]
}
############################################################




# Variable que almacena los colores a usar en el volcano plot
volcol <- c(met.brewer("Egypt",3)[1], met.brewer("Egypt",3)[2],"grey33")

# Se asigna a cada uno de los colores los posibles niveles de expresión de los genes
# de acuerdo a la columna "gene_type" agregada al marco de datos de los resultados
# de los contrastes
names(volcol) <- c("activado","reprimido","ns")

# volcano_A <- ggplot(as.data.frame(de.results), aes(x=logFC, y=-log10(adj.P.Val),
#                                                    color=gene_type,text=gene_name)) +
#   geom_point(size = 1) +
#   scale_colour_manual(values = volcol) + 
#   theme(legend.position = "none",
#         panel.background = element_rect(fill = "white"),
#         panel.grid.major = element_line(colour = "white"),
#         panel.grid.minor = element_line(colour = "white"),
#         axis.line.x.bottom = element_line(color = 'black'),
#         axis.line.y.left   = element_line(color = 'black'),
#         panel.border = element_blank(),
#         plot.title = element_text(hjust = 0.5,size = 19)) + 
#   ggtitle("Volcano plot")
# 
# volcano_A


volcano_A <- ggplot(as.data.frame(de.results), aes(x=logFC, y=-log10(adj.P.Val),
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


# Plot interactivo:
# interactive_voclano <- ggplotly(volcano_A,x = ~logFC, y = ~adj.P.Value,
#                                 text = ~paste(gene_name),
#                                 width = 640, height = 700)

interactive_voclano <- ggplotly(volcano_A,x = ~logFC, y = ~adj.P.Value,
                                tooltip = "text",
                                width = 640, height = 700)


interactive_voclano

# htmlwidgets::saveWidget(as.widget(interactive_voclano), "interactive_volcano.html")
