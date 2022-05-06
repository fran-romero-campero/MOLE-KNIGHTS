

# Cargamos los datos y asjutamos los datos
gene.expression_total <- read.table(file = "gene_expression.tsv",header=T,as.is=T)
head(gene.expression_total)
gene.ids <- gene.expression_total$geneID
gene.expression <- as.matrix(gene.expression_total[,2:ncol(gene.expression_total)])
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
library(tidyr) # Para el PCA

# Scatterplot

scatterplot_function<-function(x,y)
{
  ggplot(as.data.frame(log2(gene.expression)+1), 
         aes(x=log2(gene.expression[,x])+1, y=log2(gene.expression[,y])+1)) + 
    geom_point(size=1, color=met.brewer("Cassatt2",1) [1],size=0.4,
               aes(text= paste0("</br> X: ",round(log2(gene.expression[,x])+1,digits = 5),
                                "</br> Y: ", round(log2(gene.expression[,y])+1,digits = 5)))) +
    xlim(0,NA) +
    ylim(0,NA) +
    geom_smooth(method=lm,color="darkgreen",size=0.3) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "white"),
          panel.grid.minor = element_line(colour = "white"),
          axis.line.x.bottom = element_line(color = 'black'),
          axis.line.y.left   = element_line(color = 'black'),
          panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 13)) +
    ggtitle("Scatterplot") +
    # xlab("Replica 1") +
    # ylab("Réplica 2") +
    annotate(geom = "text", x = 1.5, y = max(c(log2(gene.expression[,x]+1),log2(gene.expression[,y]+1)))-2, 
             label = paste(c(round(100*cor(gene.expression[,x],
                                           gene.expression[,y]),
                                   digits = 2),
                             "%"), collapse=""),size=2.5)
}
  

A_1<-ggplotly(scatterplot_function(1,1),tooltip = "text")
A_2<-ggplotly(scatterplot_function(1,2),tooltip = "text")
A_3<-ggplotly(scatterplot_function(1,3),tooltip = "text")
A_4<-ggplotly(scatterplot_function(1,4),tooltip = "text")

B_1<-ggplotly(scatterplot_function(2,1),tooltip = "text")
B_2<-ggplotly(scatterplot_function(2,2),tooltip = "text")
B_3<-ggplotly(scatterplot_function(2,3),tooltip = "text")
B_4<-ggplotly(scatterplot_function(2,4),tooltip = "text")

C_1<-ggplotly(scatterplot_function(3,1),tooltip = "text")
C_2<-ggplotly(scatterplot_function(3,2),tooltip = "text")
C_3<-ggplotly(scatterplot_function(3,3),tooltip = "text")
C_4<-ggplotly(scatterplot_function(3,4),tooltip = "text")

D_1<-ggplotly(scatterplot_function(4,1),tooltip = "text")
D_2<-ggplotly(scatterplot_function(4,2),tooltip = "text")
D_3<-ggplotly(scatterplot_function(4,3),tooltip = "text")
D_4<-ggplotly(scatterplot_function(4,4),tooltip = "text")

subplot(A_1,A_2,A_3,A_4,B_1,B_2,B_3,B_4,C_1,C_2,C_3,C_4,D_1,D_2,D_3,D_4,
        nrows = 4,shareX = F,shareY = F)



# PCA 3D

# Sacamos las componentes principales para nuestra matriz de expresion y las almacenamos en df
pca_all_samples <- prcomp(x=t(log.gene.expression),center = T,scale. = F)
PCs <- as.data.frame(pca_all_samples$x)
head(PCs)

# Definimos las componentes que nos interesan
x_axis <- 1
y_axis <- 2
z_axis <- 3

# Creamos un df con informacion importante para el plot
targets <- colnames(log.gene.expression) # samples
targets <- data.frame(targets) # samples como df
targets <- separate(targets, col = targets, sep = '_', into = c('condition','replicate')) # Separamos la condicion de la replica
targets$sample <- colnames(log.gene.expression) # columna con la muestra completa
head(targets)

# Definimos el color de cada punto del PCA segun la condition
color_variable <- targets$condition
# Creamos una variable para los nombres de cada muestra
id <- targets$sample


# Definimos los nombres de los ejes
xlab <- paste0('PC',x_axis,'\n(', round((pca_all_samples$sdev^2)[x_axis]/sum(pca_all_samples$sdev^2)*100, 2), '%)')
ylab <- paste0('PC',y_axis,'\n(', round((pca_all_samples$sdev^2)[y_axis]/sum(pca_all_samples$sdev^2)*100, 2), '%)')
zlab <- paste0('PC',z_axis,'\n(', round((pca_all_samples$sdev^2)[z_axis]/sum(pca_all_samples$sdev^2)*100, 2), '%)')

# Creamos una traza vacia que usaremos en plotly
trace <- list()
i=1

# Hacemos un bucle for para crear una traza con informacion de cada color
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

trace[[2]]

# Definimos el estilo para el plot
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
  title = "PCA (3D)"
)

# Creamos el plot interactivo vacío
p <- plot_ly()

# Rellenamos el plot con la información que hemos guardado en la traza
for (x in trace){
  p <- add_trace(p, mode=x$mode, name=x$name, type=x$type, x=x$x, y=x$y, z=x$z, text=x$text)
}
# Rellenamos con el estilo del plot
p <- layout(p, scene=layout$scene, title=layout$title)
p

# Guardamos el plot en html porque no se ve en el previsualizador
setwd('~/Escritorio/') # salida del html
htmlwidgets::saveWidget(as_widget(p), "PCA3D_testeo.html", selfcontained = T)


############## PRUEBA 2 PCA 3D ###################################

prin_comp <- prcomp(t(log.gene.expression), rank. = 3)

components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components <- cbind(components, targets$condition)
head(components)

fig <- plot_ly(components, x = ~PC1, y = ~PC2, z = ~PC3, color = ~targets$condition, colors = c('#636EFA','#EF553B','#00CC96') ) %>%
  add_markers(size = 12)


fig <- fig %>%
  layout(
    title = "prueba no luis",
    scene = list(bgcolor = "#e5ecf6")
  )

fig

setwd('~/Escritorio/') # salida del html
htmlwidgets::saveWidget(as_widget(fig), "PCA3D_testeo.html", selfcontained = T)

##################################################################





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
head(log.gene.expression)
write.table(log.gene.expression, file = "log_gene-expression.tsv",sep = "\t",quote = F )
head(read.table(file = "log_gene-expression.tsv",header=T,as.is=T))



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
head(de.results)


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

head(de.results)
write.table(de.results, file = "de_results.tsv",sep = "\t" )
head(read.table(file = "de_results.tsv",header=T,as.is=T))

### CODIGO SUCIO ###########################################################
coincidencia<-subset(anotacion,anotacion[,1] == de.results[1,"gene_name"])
if(is.na(coincidencia[1,3])== T)
{
  de.results[1,"annotation"]<-"hypo"
} else {
  de.results[1,"annotation"]<-coincidencia[1,3]
}
#############################################################################

# Unimos la gene.expresion con la de.results
resultado_final<-matrix(ncol = 14)
colnames(resultado_final)<-c(colnames(gene.expression_total),colnames(de.results))

for (i in 1:nrow(gene.expression_total))
{
  coincidencia<-subset(gene.expression_total,gene.expression_total[,1] == de.results[i,"gene_name"])
  fila<-cbind(coincidencia, de.results[i,])
  resultado_final<-rbind(resultado_final, fila)
}
resultado_final<-resultado_final[-1,]
rownames(resultado_final)<-NULL
head(resultado_final)

write.table(resultado_final, file = "resultado_final.tsv",sep = "\t")
head(read.table(file = "resultado_final.tsv",header=T,as.is=T))



# Redondeamos FC y adj.P.Val: no lo necesitas
# for (i in 1:nrow(resultado_final))
# {
#   new_value<-round(resultado_final[i,"logFC"] ,digits = 6)
#   resultado_final[i,"logFC"]<-new_value
# }
# 
# for (i in 1:nrow(resultado_final))
# {
#   new_value<-round(resultado_final[i,"adj.P.Val"] ,digits = 9)
#   resultado_final[i,"adj.P.Val"]<-new_value
# }

# head(resultado_final)


# Buscar segun logFC y adj.P.Val
# fila_interesante<-subset(resultado_final,resultado_final$logFC == 5.179841)
# fila_interesante
# 
# fila_interesante2<-subset(resultado_final,resultado_final$adj.P.Val == 0.09376373)
# fila_interesante2




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
                                width = 600, height = 600)


interactive_voclano

library("xml2")
htmlwidgets::saveWidget(as.widget(interactive_voclano), "interactive_volcano.html")
read_html("interactive_volcano.html", skip = 0, remove.empty = TRUE, trim = TRUE)

# Representación de Barplot

# Funcion Fran
# barplot.gene <- function(gene.id,gene.name,gene.expression)
# {
  # expression.ll <- mean(unlist(gene.expression[gene.id,
  #                                              c("LL_1", "LL_2")]))
  # expression.hl <- mean(unlist(gene.expression[gene.id,
  #                                              c("HL_1", "HL_2")]))
  # means <- c(expression.ll,expression.hl)
  # 
  # expression.sd.ll <- sd(unlist(gene.expression[gene.id,
  #                                               c("LL_1", "LL_2")]))
  # expression.sd.hl <- sd(unlist(gene.expression[gene.id,
  #                                               c("HL_1", "HL_2")]))
  # sds <- c(expression.sd.ll,expression.sd.hl)

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

# testeo barplot
plotly(barplot.gene(gene.id="kfl00096_0240",gene.name="SAL1",gene.expression=gene.expression))

# Mi Barplot
expression.ll <- mean(unlist(gene.expression[1,c("LL_1", "LL_2")]))
expression.hl <- mean(unlist(gene.expression[1,c("HL_1", "HL_2")]))
means<-c(expression.ll,expression.hl)
expression.sd.ll <- sd(unlist(gene.expression[1,c("LL_1", "LL_2")]))
expression.sd.hl <- sd(unlist(gene.expression[1,c("HL_1", "HL_2")]))
sds<-c(expression.sd.ll,expression.sd.hl)
data<-data.frame(means)
condition<-c("LL","HL")
data<-cbind(data,condition)

positions <- c("LL", "HL")
par(lwd=3)

barcol <- c(met.brewer("Egypt",2)[1], met.brewer("Egypt",2)[2])
names(barcol) <- c("HL","LL")



p<-ggplot(data=data, aes(x=condition, y=means)) +
  scale_x_discrete(limits = positions)+
  geom_errorbar(aes(ymin = means+sds, ymax = means-sds), width = 0.2) +
  geom_bar(stat="identity",width = 0.8,fill=c("#0f7ba2","#dd5129"))+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,size = 19))+
  ggtitle(paste(c("kfl", "-", "halo"),collapse=" "))+
  labs(x="Conditions",y= "Expression")

p  

