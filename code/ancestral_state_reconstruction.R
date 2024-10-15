library(data.table)
library(ggplot2)
library(ggtree)
library(ape)
library(parallel)
library(dplyr)
library(tidytree)
library(bio3d)
library(phytools)

msh6_good<-fread("tables/S1_msh6_domains_annotated.csv")

tree_time<-read.tree("files/time_tree_may4.nwk") #from http://www.timetree.org/ # created May 4, 2024
tree_time$tip.label<-gsub("_"," ",tree_time$tip.label)
tree_time <- multi2di(tree_time)
tree_time$edge.length <- tree_time$edge.length + 1e-5

tree_time_sub <- drop.tip(tree_time, setdiff(tree_time$tip.label, msh6_good$Organism))

# pdf("Figures/time_tree.pdf", width=6.3, height=6.3)
#
# p1<-ggtree(tree_time_sub, size = 0.1, layout = "daylight", color="gray80")+
#   theme(plot.background = element_rect(fill = "transparent"),
#         panel.background = element_rect(fill = "transparent"))
# p1
#
# dev.off()


tree <- drop.tip(tree_time, setdiff(tree_time$tip.label, msh6_good$Organism))

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

tree_small<-drop.tip(tree, tree$tip.label[!tree$tip.label %in% msh6_good$Organism])
y <- as_tibble(tree_small)

d<-x %>% select(label, pwwp, tudor, none, N)
z <- left_join(y, unique(d), by = 'label')


pdf("Figures/time_treecircular.pdf", width=6.3, height=6.3)

p1<-ggtree(as.treedata(z), size = 0.1, layout = "circular", color="gray80")+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))+
  geom_tiplab(size=0.2, alpha=1)

p2<-ggtree(as.treedata(z), size = 0.4, layout = "circular", aes(col=pwwp/N, alpha=pwwp/N))+
  scale_color_continuous(high='orange2', low='gray70', guide="none")+
  scale_alpha_continuous(range = c(0,1), guide="none")+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))+
  geom_tiplab(size=0.2, alpha=1)

p3<-ggtree(as.treedata(z), size = 0.4, layout = "circular", aes(col=tudor/N, alpha=tudor/N))+
  scale_color_continuous(high='green4', low='gray70', guide="none")+
  scale_alpha_continuous(range = c(0,1), guide="none")+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))+
  geom_tiplab(size=0.2, alpha=1)

p4<-ggtree(as.treedata(z), size = 0.1, layout = "circular", color="gray80", alpha=0)+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))+
  geom_tiplab(size=0.2, alpha=1)

p1
p2
p3
p4
dev.off()


# pdf("Figures/time_tree.pdf", width=6.3, height=6.3)
#
# p1<-ggtree(as.treedata(z), size = 0.1, layout = "daylight", color="gray80")+
#   theme(plot.background = element_rect(fill = "transparent"),
#         panel.background = element_rect(fill = "transparent"))
#
# p2<-ggtree(as.treedata(z), size = 0.4, layout = "daylight", aes(col=pwwp/N, alpha=pwwp/N))+
#   scale_color_continuous(high='orange2', low='gray70', guide="none")+
#   scale_alpha_continuous(range = c(0,1), guide="none")+
#   theme(plot.background = element_rect(fill = "transparent"),
#         panel.background = element_rect(fill = "transparent"))
#
# p3<-ggtree(as.treedata(z), size = 0.4, layout = "daylight", aes(col=tudor/N, alpha=tudor/N))+
#   scale_color_continuous(high='green4', low='gray70', guide="none")+
#   scale_alpha_continuous(range = c(0,1), guide="none")+
#   theme(plot.background = element_rect(fill = "transparent"),
#         panel.background = element_rect(fill = "transparent"))
#
# p4<-ggtree(as.treedata(z), size = 0.1, layout = "daylight", color="gray80", alpha=0)+
#   theme(plot.background = element_rect(fill = "transparent"),
#         panel.background = element_rect(fill = "transparent"))+
#   geom_tiplab(size=0.2, alpha=1)+
#   geom_nodelab(size=0.2, col="black")
#
# p1
# p2
# p3
# p4
# dev.off()


pdf("Figures/ancestralstate.pdf", width=6.3, height=6.3)

x<-msh6_good[Organism %in% tree_time$tip.label]$tip=="Tudor"
names(x)<-msh6_good[Organism %in% tree_time$tip.label]$Organism
XX <- contMap(tree, x)
plot(XX, type = "fan", lwd=0.25,fsize=.1,outline=FALSE, plot=F)


x<-msh6_good[Organism %in% tree_time$tip.label]$tip=="PWWP"
names(x)<-msh6_good[Organism %in% tree_time$tip.label]$Organism
XX <- contMap(tree, x)
plot(XX, type = "fan", lwd=0.25,fsize=.1,outline=FALSE, plot=F)

dev.off()



plot_ancestral_states<-function(x, colors=c("gray", "orange")){
  names(x)<-msh6_good[Organism %in% tree_time$tip.label]$Organism
  testtree<-drop.tip(tree, tree$tip.label[!tree$tip.label %in% names(x)])

  fitace<-ace(x,testtree,type="discrete",model="ARD")
  AIC(fitace)

  aic_value<-AIC(fitace)
  plotTree(testtree, setEnv = TRUE, offset = 5, type="fan", fsize=0.1, lwd=0.25)
  nodelabels(node = as.numeric(rownames(fitace$lik.anc)), pie = fitace$lik.anc,
             piecol = colors, cex = 0.25, lwd=0.25)
  title(main="ARD", sub=aic_value)
  message(paste("ARD",aic_value))

  fitace<-ace(x,testtree,type="discrete",model="ER")
  aic_value<-AIC(fitace)
  plotTree(testtree, setEnv = TRUE, offset = 5, type="fan", fsize=0.1, lwd=0.25)
  nodelabels(node = as.numeric(rownames(fitace$lik.anc)), pie = fitace$lik.anc,
             piecol = colors, cex = 0.25, lwd=0.25)
  title(main="ER", sub=aic_value)
  message(paste("ER",aic_value))
}

pdf("Figures/ancestralstate_piecharts.pdf", width=7, height=7)

tree<-tree_time
x<-msh6_good[Organism %in% tree_time$tip.label]$tip
plot_ancestral_states(x, c("gray", "orange2", "green4"))

#PWWP
x<-msh6_good[Organism %in% tree_time$tip.label]$tip=="PWWP"
plot_ancestral_states(x, c("gray", "orange2"))

#Tudor
x<-msh6_good[Organism %in% tree_time$tip.label]$tip=="Tudor"
plot_ancestral_states(x, c("gray", "green4"))
par(mar=c(5.1,4.1,4.1,2.1))
dev.off()


# plot ancestral state of traits ------------------------------------------

library(ggtree)
library(treeio)
library(tidytree)
library(ggplot2)
library(TDbook)

ancestral_state_trait_tree<-function(trait){

  tree_time<-read.tree("files/species_for_tree2.nwk") #from http://www.timetree.org/
  tree_time$tip.label<-gsub("_"," ",tree_time$tip.label)
  tree_time <- multi2di(tree_time)
  tree_time$edge.length <- tree_time$edge.length + 1e-5

  tree_time_sub <- drop.tip(tree_time, setdiff(tree_time$tip.label, trait_data$Organism))

  svl <- as.matrix(data.frame(trait_data[,trait, with=F]))
  names(svl)<-trait_data$Organism
  svl<-svl[names(svl) %in% tree_time_sub$tip.label]
  svl<-log10(svl)
  svl<-svl[is.finite(svl)]

  tree_time_sub <- drop.tip(tree_time_sub, setdiff(tree_time_sub$tip.label, names(svl)))
  fit <- phytools::fastAnc(tree_time_sub, svl, vars=TRUE, CI=TRUE)

  td <- data.frame(node = nodeid(tree_time_sub, names(svl)),
                   trait = svl)
  td$reader<-trait_data$tip[match(row.names(td), trait_data$Organism)]
  nd <- data.frame(node = names(fit$ace), trait = fit$ace)
  nd$reader<-NA
  d <- rbind(td, nd)
  d$node <- as.numeric(d$node)
  tree <- full_join(tree_time_sub, d, by = 'node')

}

plot_ancestral_trait<-function(trait){
  tree<-ancestral_state_trait_tree(trait=trait)
  p0 <- ggtree(tree, aes(color=trait), layout = 'circular',
               ladderize = FALSE, continuous = 'colour', size=0.1) +
    scale_color_gradientn(colours=c("white", 'black')) +
    geom_tiplab(hjust = -.1, size=0) +
    theme(plot.background = element_blank())

  p1 <- ggtree(tree, aes(color=trait), layout = 'circular',
               ladderize = FALSE, continuous = 'colour', size=0.1) +
    scale_color_gradientn(colours=c("white", 'black')) +
    geom_tiplab(hjust = -.1, size=0) +
    theme(legend.position = "none", plot.background = element_blank(),
          panel.background = element_blank())

  p2 <- ggtree(tree, layout = 'circular',
               ladderize = FALSE, size=0.1, alpha=0) +
    scale_color_manual(values =c(PWWP="orange2", Tudor='green3', None='white')) +
    geom_tippoint(size=0, aes(color=reader)) +
    theme(legend.position = "none", plot.background = element_blank(),
          panel.background = element_blank())


  p3 <- ggtree(tree, aes(color=trait), layout = 'circular',
               ladderize = FALSE, continuous = 'colour', size=0.1) +
    scale_color_gradientn(colours=c("yellow", 'dodgerblue')) +
    geom_tiplab(hjust = -.1, size=0) +
    theme_void(base_size = 6)+
    theme(legend.position = "left", plot.background = element_blank(), legend.key.size = unit(0.25, "cm"))

  return(list(p0,
  p1,
  p2, p3))
}


# pdf("Figures/intron_length_tree_plot.pdf", width=1.75, height=1.75)
#
# plots<-plot_ancestral_trait("mean_intron_length")
# lapply(plots, plot)
#
# dev.off()

pdf("Figures/mean_exon_number_tree_plot.pdf", width=1.75, height=1.75)
  plots<-plot_ancestral_trait("mean_exon_number")
  lapply(plots, plot)
dev.off()

pdf("Figures/mean_CDS_length_tree_plot.pdf", width=3, height=3)
plots<-plot_ancestral_trait("mean_CDS_length")
lapply(plots, plot)
dev.off()

pdf("Figures/mean_intron_length_tree_plot.pdf", width=3, height=3)
plots<-plot_ancestral_trait("mean_intron_length")
lapply(plots, plot)
dev.off()

pdf("Figures/mean_exon_number_log10_binaryPGLMM.pdf", width=1.5, height=1.5)

  load(paste0("RData_results/", "mean_exon_number","_glm.Rda"))

  coefs <- -(out$log10_binaryPGLMM$B[,1])

  trait_values<-log10(trait_data$mean_exon_number)
  trait_values<-trait_values[is.finite(trait_values)]
  x_plot <- seq(min(trait_values), max(trait_values), by = 0.1)
  y_plot <- plogis(coefs[1] + coefs[2] * x_plot)

  plot_data <- data.frame(x_plot=10^x_plot, reader=y_plot)

  trait_data_tmp<-trait_data
  trait_data_tmp$x_plot<-(trait_data_tmp$mean_exon_number)
  trait_data_tmp$reader<-as.numeric(trait_data_tmp$tip!="None")
  ggplot(plot_data) +
    geom_line(aes(x_plot, reader), col = "purple4") +
    xlab("x") + ylab("p(Reader | log10(# exons/gene))") +
    scale_y_continuous(limits = c(0, 1)) + theme_bw()+
    scale_x_log10(name="")+
    geom_jitter(data=trait_data_tmp, aes(x=x_plot, y=reader, col=tip),width=0,
                height=0.01, shape=20, alpha=0.5, size=0.75)+
    geom_line(aes(x_plot, reader), col = "purple4") +
    scale_color_manual(values =c(PWWP="orange2", Tudor='green3', None='gray'), guide="none")+
    theme_bw(base_size = 6)+
    theme(panel.grid = element_blank())

dev.off()

