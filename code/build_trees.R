# Create a  tree

tree<-read.tree("files/phyloT_generated_tree_Aug12_23.txt")
tree$tip.label<-gsub("_"," ",tree$tip.label)
tree<-drop.tip(tree, tree$tip.label[!tree$tip.label %in% msh6_good$Organism])
x <- as_tibble(tree)
x$pwwp<-sapply(x$label, function(n){
  if(n %in% tree$node.label){
  clade<-extract.clade(tree, n)
  sub_msh6_good<-msh6_good[Organism %in% clade$tip.label]
  pwwp<-sum(sub_msh6_good$tipcolor=="PWWP")
  } else {
    sub_msh6_good<-msh6_good[Organism %in% n]
    pwwp<-sum(sub_msh6_good$tipcolor=="PWWP")
  }
  return(pwwp)
})

x$tudor<-sapply(x$label, function(n){
  if(n %in% tree$node.label){
    clade<-extract.clade(tree, n)
    sub_msh6_good<-msh6_good[Organism %in% clade$tip.label]
    tudor<-sum(sub_msh6_good$tipcolor=="Tudor")
  } else {
    sub_msh6_good<-msh6_good[Organism %in% n]
    tudor<-sum(sub_msh6_good$tipcolor=="Tudor")
  }
  return(tudor)
})

x$none<-sapply(x$label, function(n){
  if(n %in% tree$node.label){
    clade<-extract.clade(tree, n)
    sub_msh6_good<-msh6_good[Organism %in% clade$tip.label]
    None<-sum(sub_msh6_good$tipcolor=="None")
  } else {
    sub_msh6_good<-msh6_good[Organism %in% n]
    None<-sum(sub_msh6_good$tipcolor=="None")
  }
  return(None)
})

x$N<-sapply(x$label, function(n){
  if(n %in% tree$node.label){
    clade<-extract.clade(tree, n)
    sub_msh6_good<-msh6_good[Organism %in% clade$tip.label]
    N<-nrow(sub_msh6_good)
  } else {
    sub_msh6_good<-msh6_good[Organism %in% n]
    N<-nrow(sub_msh6_good)
  }
  return(N)
})

tree_small<-drop.tip(tree, tree$tip.label[!tree$tip.label %in% small_tree$Organism])
y <- as_tibble(tree_small)

d<-x %>% select(label, pwwp, tudor, none, N)
z <- left_join(y, unique(d), by = 'label')


pdf("figures/MSH6_trees_daylight_new.pdf", width=6.3, height=6.3)

p1<-ggtree(as.treedata(z), size = 0.1, layout = "daylight", color="gray80")+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))

p2<-ggtree(as.treedata(z), size = 0.4, layout = "daylight", aes(col=pwwp/N, alpha=pwwp/N))+
  scale_color_continuous(high='orange', low='gray70', guide="none")+
  scale_alpha_continuous(range = c(0,1), guide="none")+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))

p3<-ggtree(as.treedata(z), size = 0.4, layout = "daylight", aes(col=tudor/N, alpha=tudor/N))+
  scale_color_continuous(high='green4', low='gray70', guide="none")+
  scale_alpha_continuous(range = c(0,1), guide="none")+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))

p4<-ggtree(as.treedata(z), size = 0.1, layout = "daylight", color="gray80", alpha=0)+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))+
  geom_tiplab(size=0.2, alpha=1)+
  geom_nodelab(size=0.2, col="black")

p1
p2
p3
p4
dev.off()


pdf("~/Documents/MSH6_trees_circular_new.pdf", width=6.3, height=6.3)

p0<-ggtree(as.treedata(z), size = 0.1, layout = "circular", color="gray80")+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))+
  geom_cladelab(node=3314, label="Fungi", align=T,
                geom='text', fill='lightblue', offset=.5, vjust=-3, col="gray")
  geom_hilight(node=2606,alpha=1,fill="lightblue")+
  geom_hilight(node=3575,alpha=1,fill="green")+
  geom_tree(size = 0.1, layout = "circular", color="gray80")
p0

p1<-ggtree(as.treedata(z), size = 0.1, layout = "circular", color="gray80")+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))+
  geom_tiplab(size=0.2)

p2<-ggtree(as.treedata(z), size = 0.4, layout = "circular", aes(col=pwwp/N, alpha=pwwp/N))+
  scale_color_continuous(high='orange2', low='gray70', guide="none")+
  scale_alpha_continuous(range = c(0,1), guide="none")+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))+
  geom_tiplab(size=0.2)


p3<-ggtree(as.treedata(z), size = 0.4, layout = "circular", aes(col=tudor/N, alpha=tudor/N))+
  scale_color_continuous(high='green4', low='gray70', guide="none")+
  scale_alpha_continuous(range = c(0,1), guide="none")+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))+
  geom_tiplab(size=0.2)
p4<-ggtree(as.treedata(z), size = 0.1, layout = "circular", color="gray80", alpha=0)+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))+
  geom_tiplab(size=0.2, alpha=0)+
  geom_nodelab(size=0.2, col="black")

p0
p1
p2
p3
p4
dev.off()

pdf("~/Documents/MSH6_trees_circular_big_new.pdf", width=15, height=15)

p1<-ggtree(as.treedata(z), size = 0.1, layout = "circular", color="gray80")+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))+
  geom_tiplab(size=0.2)

p2<-ggtree(as.treedata(z), size = 0.4, layout = "circular", aes(col=pwwp/N, alpha=pwwp/N))+
  scale_color_continuous(high='orange2', low='gray70', guide="none")+
  scale_alpha_continuous(range = c(0,1), guide="none")+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))+
  geom_tiplab(size=0.2)


p3<-ggtree(as.treedata(z), size = 0.4, layout = "circular", aes(col=tudor/N, alpha=tudor/N))+
  scale_color_continuous(high='green4', low='gray70', guide="none")+
  scale_alpha_continuous(range = c(0,1), guide="none")+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))+
  geom_tiplab(size=0.2)
p4<-ggtree(as.treedata(z), size = 0.1, layout = "circular", color="gray80", alpha=0)+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))+
  geom_tiplab(size=0.2, alpha=0)+
  geom_nodelab(size=0.2, col="black")

p1
p2
p3
p4
dev.off()








pdf("~/Documents/MSH6_trees_circular_big_new2.pdf", width=18, height=18)

p1<-ggtree(as.treedata(z), size = 0.2, layout = "circular", color="gray80")+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))+
  geom_hilight(node=dtz$node[dtz$label=="Viridiplantae"],alpha=0.1,fill="forestgreen")+
  geom_hilight(node=dtz$node[dtz$label=="Fungi"],alpha=0.1,fill="purple")+
  geom_hilight(node=dtz$node[dtz$label=="Metazoa"],alpha=0.1,fill="orange")+
  geom_tiplab(size=0.3)+
  geom_nodelab(size=.3, col="black")

p1

p1<-ggtree(as.treedata(z), size = 0.2, layout = "daylight", color="gray80")+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))+
  geom_hilight(node=dtz$node[dtz$label=="Viridiplantae"],alpha=0.1,fill="forestgreen", type="encircle")+
  geom_hilight(node=dtz$node[dtz$label=="Fungi"],alpha=0.1,fill="purple", type="encircle")+
  geom_hilight(node=dtz$node[dtz$label=="Metazoa"],alpha=0.1,fill="orange", type="encircle")+
  geom_tiplab(size=0.3)+
  geom_nodelab(size=.3, col="black")

p1


dev.off()

