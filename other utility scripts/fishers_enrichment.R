## fishers test
m=length(query.genes) ## set1 genes
k2=length(compare.genes) ##set2 genes
x=length(intersect(query.genes, compare.genes))  ##overlap bt set1 and set2
n=length(genes)-m #universe minus set1 genes
##fisher exact test
fshr.mtx <- c(x, k2-x, m-x, n-(k2-x))
ft <-fisher.test(matrix(data=fshr.mtx, nrow=2), alternative="g")  ## default is 'two.sided', 'g'=enrichment, 'l'=depletion



## permutation test
n.total=length(genes)  ## universe genes
n.set1=length(query.set) ## set 1
n.set2=length(compare.set) ## set 2
n.overlap=length(intersect(query.set, compare.set)) ## overlap bt set1 and set2
nits=100000; ## number of iterations to run
permute.out=c(); 
for (i in 1:nits) {
  permute.out=c(permute.out, length(intersect(sample(n.total, n.set1), sample(n.total, n.set2))))
}
hist(permute.out, breaks=50) ## histogram of permutation overlaps
length(which(permute.out>n.overlap))  ## how many permutation overlaps are greater than the set1/set2 overlap?
length(which(permute.out>n.overlap))/nits ## permutation pval





       
