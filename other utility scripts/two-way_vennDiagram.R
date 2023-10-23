library(VennDiagram)

set1=geneSets$DE_Genes
set2=geneSets$DCGs_DiffCoExp
grid.newpage()
draw.pairwise.venn(area2=length(set1), 
                   area1=length(set2), 
                   cross.area=length(intersect(set2, set1)), 
                   category = c("Differentially\nExpressed Genes\n(n=1,055)", 
                                "Differentially\nCo-Expressed Genes\n(n=777)"), 
                   col=c("black","black"),
                   lwd=c(2,2),
                   fill = c("lightblue", "pink"), 
                   alpha = rep(0.5, 2), 
                   cex=c(2,2,2),
                   fontfamily = rep("sans", 3),
                   cat.pos = c(315,45), 
                   cat.cex = c(1.5,1.5),
                   cat.dist = rep(0.2, 2), 
                   cat.fontfamily = rep("sans", 2),
                   scaled = FALSE)
