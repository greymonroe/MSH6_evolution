#MSH6_tree<-read.tree("~/Downloads/msh6_domain_sequences_aligned_rapidnj.tree")
#MSH6_tree<-read.tree("files/msh6_domain_sequences_plants_aligned.contree")
MSH6_tree<-read.tree("~/Downloads/msh6_domain_sequences_aligned.contree")

predicted_msh6_blast_annotated$lab<-paste0(predicted_msh6_blast_annotated$protein,"-", gsub(" ", "_", predicted_msh6_blast_annotated$Organism))

annot_plants<-organisms_in_node(tree, predicted_msh6_blast_annotated, "Viridiplantae")
annot_other<-organisms_in_node(tree, predicted_msh6_blast_annotated, c("Fungi","Metazoa","Viridiplantae"), exclude = T)
annot_sar<-organisms_in_node(tree, predicted_msh6_blast_annotated, c("Sar"))

MSH6_tree<-drop.tip(MSH6_tree, MSH6_tree$tip.label[!MSH6_tree$tip.label %in% predicted_msh6_blast_annotated$lab])
MSH6_tree<-drop.tip(MSH6_tree, MSH6_tree$tip.label[!MSH6_tree$tip.label %in% c(annot_plants$lab, annot_other$lab)])

# grep("Arab", MSH6_tree$tip.label, value=T) %in%
# grep("Arab", annot_plants$lab, value=T)

x <- as_tibble(MSH6_tree)
x$plant<-x$label %in% annot_plants$lab
x$other<-x$label %in% annot_other$lab
x$SAR<-x$label %in% annot_sar$lab

x$tudor<-x$label %in% predicted_msh6_blast_annotated[Tudor_domain==T]$lab

pdf("Figures/msh6_treeplant.pdf", width=7, height=7)

p0<-ggtree(as.treedata(x), size = 0.1, layout = "circular", color="gray80")+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))+
  geom_tree(size = 0.1, layout = "circular", color="gray80")+
  geom_tippoint(aes(col=tudor), size=0.15)+
  #geom_tippoint(data=as.treedata(x[x$tudor==T,]), shape=21, col="black", size=0.15, fill=NA)+
  scale_color_manual(values = c("gray90","green3"), name="Tudor")+
  geom_tiplab(size=0.5, aes(col=tudor))+
  theme(legend.position = "none")
p0

dev.off()

