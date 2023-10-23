library(VennDiagram)


## Three Way Venn Diagram
set1 <- genes[which(modgenes_df$Unstim=="M2")]
set2 <- geneSets$LVMASS.1Mb
set3 <- geneSets$DCGs_DiffCoExp
grid.newpage()
venn <- draw.triple.venn(
  area1=length(set1),
  area2=length(set2),
  area3=length(set3),
  n12=length(intersect(set1, set2)),
  n23=length(intersect(set2, set3)),
  n13=length(intersect(set1, set3)),
  n123=length(intersect(set1, intersect(set2, set3))),
  category = c("Baseline Module 2 Genes\n(n=2,397)", 
               "LV Mass GWAS Genes\n(n=1,169)",
               "Differentially CoExpressed Genes\n(n=777)"),
  col=c("black","black","black"),
  lwd=c(2,2,2),
  fill = c("yellow", "lightblue","pink"), 
  alpha = rep(0.5, 3), 
  cex=rep(2, 7),
  fontfamily = rep("sans", 7),
  cat.pos = c(345,15,180), 
  cat.cex = c(1.5,1.5,1.5),
  cat.dist = rep(0.095, 3), 
  cat.fontfamily = rep("sans", 3),
  euler.d = F,
  scaled = FALSE,
)