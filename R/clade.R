library(picante)


clade_name = commandArgs(TRUE)[1]
path_phylogeny=commandArgs(TRUE)[2]	### phylo d'un profil
path_otuFile=commandArgs(TRUE)[3]	### otus d'un profil
profile_fig=commandArgs(TRUE)[4]

phylogeny = read.tree(path_phylogeny)	
otuFile = read.table(path_otuFile)

#mntd_CladeTree=commandArgs(TRUE)[5]
#depth_deepest=commandArgs(TRUE)[6]

file = commandArgs(TRUE)[5]
############

otuVector = as.vector(otuFile)  ## récupérer les noms des OTUs affiliés à un arbre phylo
t_otuVector = t(otuVector)	## transposée
otuNumber = length(t_otuVector)	## nombre d'otus

otuMatrix=matrix(1,1,length(t_otuVector))	## générer une matrice [1, nombre d'OTUs], avec des uns comme valeurs (abondance = 1)
colnames(otuMatrix)= t_otuVector	### renommer les colonnes de la matrice avec les noms des OTUs


#############

CladeTree=drop.tip(phylogeny, which(!phylogeny$tip.label %in% t_otuVector)) ## extraire le sous arbre des otus
patristic_CladeTree = cophenetic(CladeTree)

mntd_CladeTree = round(mntd(otuMatrix, patristic_CladeTree, abundance.weighted=T), digits=4)
depth_deepest = round(sum(CladeTree$edge.length)/length(CladeTree$tip.label), digits=4)

sink(file)
print (mntd_CladeTree)
print (depth_deepest)
sink()

#### fig

svg(profile_fig)
zoom(phylogeny, t_otuVector,  subtree=T, col="red", no.margin=F, root.edge=T, show.node.label = F, cex=0.4)
dev.off()

###########


