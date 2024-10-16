### Investigating MSH7

msh7<-potential[grepl("MUTS homolog 7",arabidopsis_name) & grepl("DNA mismatch repair protein Msh6", human_name)]
msh7$Organism <-msh7$organism
interpro_msh7<-interpro[Protein %in% msh7$Protein]

msh7_plants<-organisms_in_node(tree, msh7, c("Viridiplantae"))
msh7_plants[Protein %in% interpro_msh7[grepl("Tudor",Signature_Description)]$Protein]

msh7_Sar<-organisms_in_node(tree, msh7, c("Sar"))
msh7_Sar[Protein %in% interpro_msh7[grepl("Tudor",Signature_Description)]$Protein]


tree<-read.tree("data/phyloT_generated_tree_Aug12_23.txt") # NCBI Taxonomy
tree$tip.label<-gsub("_"," ",tree$tip.label)
tree<-drop.tip(tree, tree$tip.label[!tree$tip.label %in% msh6_good$Organism])

x <- as_tibble(tree)
x$MSH7<-sapply(x$label, function(n){
  if(n %in% tree$node.label){
    clade<-extract.clade(tree, n)
    sub<-msh7[Organism %in% clade$tip.label]
    MSH7<-uniqueN(sub$Organism)
  } else {
    sub<-msh7[Organism %in% n]
    MSH7<-uniqueN(sub$Organism)
  }
  return(MSH7)
})

x$N<-sapply(x$label, function(n){
  if(n %in% tree$node.label){
    clade<-extract.clade(tree, n)
   N<-length(clade$tip.label)
  } else {

    N<-1
  }
  return(N)
})

tree_small<-drop.tip(tree, tree$tip.label[!tree$tip.label %in% small_tree$Organism])
y <- as_tibble(tree_small)

d<-x %>% select(label, MSH7, N)
z <- left_join(y, unique(d), by = 'label')


pdf("Figures/MSH6_trees_daylight_new.pdf", width=3.5, height=3.5)
# Plot the pruned tree (your plot commands)
p1 <- ggtree(as.treedata(z), size = 0.1, layout = "circular", color="gray90") +
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))

p2 <- ggtree(as.treedata(z), size = 0.4, layout = "circular", aes(col=MSH7/N, alpha=log(MSH7))) +
  scale_color_continuous(high='red', low='gray90', guide="none") +
  scale_alpha_continuous(range = c(0,1), guide="none") +
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))
ggdraw() +
  draw_plot(p1) +
  draw_plot(p2)

dev.off()
z_dt<-data.table(z)

z_dt$MSH7_pct<-z_dt$MSH7/z_dt$N

z_dt[label=="Metazoa"]
z_dt[label=="Fungi"]
z_dt[label=="Sar"]
z_dt[label=="Viridiplantae"]

