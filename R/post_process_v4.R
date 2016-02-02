path_resultat=commandArgs(TRUE)[1]
setwd(path_resultat)


library(phyloseq)
library(vegan)



#Données brutes

otu_table=read.table("OTU.txt",h=T,sep="\t")
tax_tableNN=read.table("taxo_NN.csv", h=T, sep=";",, fill=TRUE, quote = "\"")
rownames(tax_tableNN)=rownames(otu_table)
otu_nn=phyloseq(otu_table(otu_table, taxa_are_rows=T), tax_table(as.matrix(tax_tableNN[,1:7])))	

write.table(round(estimate_richness(otu_nn),2), "../Richness_diversity.tmp",sep="\t", quote=F, col.names=T, row.names=T)

jpeg("Richness_Diversity_OTUs.jpg")
plot_richness(otu_nn,  measures=c("Chao1", "ACE", "Shannon"))
	
jpeg("Rarefaction_curves_OTUs.jpg")
rarecurve(t(otu_table), step=100,  xlab = "Sequences", ylab = "OTUs")


if(ncol(otu_table)>2)
{
jpeg("Heatmap_OTUs.jpg")
plot_heatmap(otu_nn, "MDS", "bray", low = "#66CCFF", high = "#000033", na.value = "white") 
}



#Données normalisées

otu_table_norm_NN=read.table("OTU_normalized_NN.txt",h=T,sep="\t")
#otu_table_norm_LCA=read.table("OTU_normalized_LCA.txt",h=T,sep="\t")
otu_norm_nn=phyloseq(otu_table(otu_table_norm_NN, taxa_are_rows=T))


#otu_norm_nn=phyloseq(otu_table(otu_table_norm_NN, taxa_are_rows=T), tax_table(as.matrix(tax_tableNN[,1:7]))) # Ne fonctionne pas car la table normaliése n'a pas le même  nb de lignes
#otu_norm_lca=phyloseq(otu_table(otu_table_norm_LCA, taxa_are_rows=T)) # données identiques aux précédentes


if(sum(otu_table_norm_NN))
{
write.table(round(estimate_richness(otu_norm_nn),2), "../Richness_diversity_norm.tmp",sep="\t", quote=F, col.names=T, row.names=T)

jpeg("Rarefaction_curves_OTUs_normalized.jpg")
rarecurve(t(otu_table_norm_NN), step=100,  xlab = "Sequences", ylab = "OTUs")

jpeg("Richness_Diversity_OTUs_normalized.jpg")
plot_richness(otu_norm_nn,  measures=c("Chao1", "ACE", "Shannon"))
}

if(ncol(otu_table_norm_NN)>2)
{
jpeg("Heatmap_OTUs_normalized.jpg")
plot_heatmap(otu_norm_nn,"MDS", "bray", low = "#66CCFF", high = "#000033", na.value = "white")
}



dev.off()

save.image()


# derniers changements : le 29/08/2014
