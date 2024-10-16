
plot_protein_models<-function(domain_search, predicted_msh6_blast_annotated, interpro, color){
  interpro$Organism<-predicted_msh6_blast_annotated$Organism[match(interpro$Protein, predicted_msh6_blast_annotated$protein)]
  interpro_sub<-interpro[Protein %in% predicted_msh6_blast_annotated$protein & Analysis!="SUPERFAMILY" & grepl(domain_search, Signature_Description, ignore.case = T)]
  interpro_sub<-interpro[Analysis %in% c("CDD", "Pfam","SMART", "ProSiteProfiles")
                         & Protein %in% interpro_sub$Protein
                         & (grepl("MutS", Signature_Description)
                            | grepl(domain_search, Signature_Description, ignore.case = T))]
  table(interpro_sub$Analysis, interpro_sub$Signature_Description)

  interpro_sub$Signature_Description[grepl(domain_search, interpro_sub$Signature_Description, ignore.case = T)]<-"Domain"
  interpro_sub$Name<-paste(interpro_sub$Organism, interpro_sub$Protein, sep=" | ")
  table(interpro_sub$Signature_Description)
  best_model<-interpro_sub[,.(I=sum(Signature_Description=="MutS domain I"),
                              II=sum(Signature_Description=="MutS domain II"),
                              III=sum(Signature_Description=="MutS domain III"),
                              IV=sum(Signature_Description=="MutS family domain IV"),
                              V=sum(Signature_Description=="MutS domain V"),
                              domain_search=sum(Signature_Description==domain_search),
                              Length=max(Stop_Location)), by=.(Organism, Protein)]


  # Count the number of zeros in columns I, II, III, IV, and V
  best_model[, zero_count := rowSums(.SD == 0), .SDcols = c("I", "II", "III", "IV", "V")]
  best_model[, greater_than_one_count := rowSums(.SD > 1), .SDcols = c("I", "II", "III", "IV", "V")]

  # Order by zero_count and Length
  setorder(best_model, Organism, zero_count, Length)

  # Subset to get the top row for each Organism
  result <- best_model[, .SD[1], by = Organism][zero_count<1 &greater_than_one_count<1 ]
  results_all<-best_model[zero_count<1 &greater_than_one_count<1 ]


  interpro_sub_result<-interpro_sub[Protein %in% results_all$Protein]

  domain_dist1<-interpro_sub_result[Signature_Description=="MutS domain I", .(dist=min(Start_Location)), by=.(Organism=Organism, Protein=Protein)]

  domain_distance<-interpro_sub_result[Signature_Description=="Domain", .(dist=max(Stop_Location), start=min(Start_Location)), by=.(Organism=Organism,Protein=Protein)]
  domain_distance$diff<-domain_dist1$dist-domain_distance$dist

  if(nrow(result)>200){
    result<-result[sample(1:nrow(result), 200)]
  }

  interpro_sub_result<-interpro_sub[Protein %in% result$Protein]

  domain_dist<-interpro_sub_result[Signature_Description=="MutS domain I", .(dist=min(Start_Location)), by=.(Organism=Organism)]
  interpro_sub_result$dist<-domain_dist$dist[match(interpro_sub_result$Organism, domain_dist$Organism)]

  interpro_sub_result$Signature_Description<-factor(interpro_sub_result$Signature_Description,
                                                    levels=c("Domain", "MutS domain I", "MutS domain II", "MutS domain III", "MutS family domain IV", "MutS domain V"))
  interpro_sub_result<-interpro_sub_result[order(Signature_Description)]
  domain_plot<-ggplot(interpro_sub_result, aes(x=Start_Location-dist, xend=Stop_Location-dist, y=Name, yend=Name, col=Signature_Description))+
    geom_segment(lineend = "round", size=.5)+
    theme_bw(base_size = 6)+
    theme(legend.position = "none", axis.title = element_blank(), panel.grid = element_blank())+
    scale_color_manual(values=c(color, colorRampPalette(c("lightblue","darkblue"))(5)))
  domain_plot

  domain_plot2<-ggplot(interpro_sub_result, aes(x=Start_Location-dist, xend=Stop_Location-dist, y=Name, yend=Name, col=Signature_Description))+
    geom_segment(lineend = "round", size=.1)+
    theme_bw(base_size = 6)+
    theme_void()+
    theme(legend.position = "none", axis.title = element_blank())+
    scale_color_manual(values=c(color, colorRampPalette(c("lightblue","darkblue"))(5)))
  domain_plot2

  return(list(domain_plot, domain_plot2, result=result, interpro_sub=interpro_sub, interpro_sub_result=interpro_sub_result, best_model=best_model, results_all=results_all, domain_distance=domain_distance))
}

pdf("Figures/MSH6_PWWP_models.pdf", width=2, height=1)
plot_protein_models(domain_search = "pwwp", predicted_msh6_blast_annotated, interpro, color="orange")
dev.off()

pdf("Figures/MSH6_ESCO_models.pdf", width=3, height=1)
plot_protein_models("ESCO", predicted_msh6_blast_annotated, interpro, color="red")
dev.off()

pdf("Figures/MSH6_Tudor_models.pdf", width=2, height=1)
plot_protein_models("tudor", predicted_msh6_blast_annotated, interpro, color = "green4")
dev.off()

pdf("Figures/MSH6_domains_models.pdf", width=2, height=1)
plot_protein_models("pwwp", predicted_msh6_blast_annotated, interpro, color="orange2")
plot_protein_models("tudor", predicted_msh6_blast_annotated, interpro, color = "green4")
dev.off()

# Archived
# pdf("Figures/MSH6_PWWP_models_big.pdf", width=5, height=10)
# plot_protein_models("pwwp", predicted_msh6_blast_annotated, interpro, color="orange2")
# dev.off()
#
# pdf("Figures/MSH6_ESCO_models_big.pdf",width=5, height=10)
# plot_protein_models("ESCO", predicted_msh6_blast_annotated, interpro, color="red")
# dev.off()
#
# pdf("Figures/MSH6_Tudor_models_big.pdf", width=5, height=10)
# plot_protein_models("tudor", predicted_msh6_blast_annotated, interpro, color = "green4")
# dev.off()


#compare sizes:
PWWP<-plot_protein_models(domain_search = "pwwp", predicted_msh6_blast_annotated, interpro, color="purple")
PWWP$domain_distance$domain<-"PWWP"
ggplot(PWWP$domain_distance, aes(x=diff))+
  geom_histogram()
tudor<-plot_protein_models("tudor", predicted_msh6_blast_annotated, interpro, color = "green4")
tudor$domain_distance$domain<-"Tudor"
alldist<-rbind(tudor$domain_distance, PWWP$domain_distance)
alldist$order<-rank(alldist$diff)

alldist95<-alldist[order>nrow(alldist)*.025 & order<nrow(alldist)*.975]
ggplot(alldist95, aes(x=diff, fill=domain))+
  geom_density(alpha=0.5)

t.test(alldist$diff~alldist$domain)
t.test(alldist$start~alldist$domain)

alldist95[,.(mean(diff), sd(diff)), by=domain]
tudor$domain_distance$meandiff<-abs(tudor$domain_distance$diff-mean(tudor$domain_distance$diff))
PWWP$domain_distance$meandiff<-abs(PWWP$domain_distance$diff-mean(PWWP$domain_distance$diff))


# ESCO Tree ---------------------------------------------------------------
esco_details<-plot_protein_models(domain_search = "ESCO, zinc-finger", predicted_msh6_blast_annotated, interpro, color="red")$interpro_sub

esco_species<-plot_protein_models("ESCO", predicted_msh6_blast_annotated, interpro, color="red")$result

fungi_esco<-organisms_in_node(tree, msh6_good, "Fungi")
fungi_esco$tipcolor[fungi_esco$Organism %in% esco_species$Organism]<-"ESCO"

tree<-read.tree("data/phyloT_generated_tree_Aug12_23.txt")
tree$tip.label<-gsub("_"," ",tree$tip.label)


esco_tree<-drop.tip(tree, tree$tip.label[!tree$tip.label %in% fungi_esco$Organism])

x <- as_tibble(esco_tree)
x$ESCO<-sapply(x$label, function(n){
  if(n %in% esco_tree$node.label){
    clade<-extract.clade(esco_tree, n)
    sub_msh6_good<-fungi_esco[Organism %in% clade$tip.label]
    ESCO<-sum(sub_msh6_good$tipcolor=="ESCO")
  } else {
    sub_msh6_good<-fungi_esco[Organism %in% n]
    ESCO<-sum(sub_msh6_good$tipcolor=="ESCO")
  }
  return(ESCO)
})


x$N<-sapply(x$label, function(n){
  if(n %in% esco_tree$node.label){
    clade<-extract.clade(esco_tree, n)
    sub_msh6_good<-fungi_esco[Organism %in% clade$tip.label]
    N<-nrow(sub_msh6_good)
  } else {
    sub_msh6_good<-fungi_esco[Organism %in% n]
    N<-nrow(sub_msh6_good)
  }
  return(N)
})


tree_small<-drop.tip(esco_tree, esco_tree$tip.label[!esco_tree$tip.label %in% fungi_esco$Organism])
y <- as_tibble(esco_tree)

d<-x %>% select(label, ESCO, N)
z <- left_join(y, unique(d), by = 'label')


pdf("Figures/MSH6_trees_circular_ESCO.pdf", width=3, height=3)

p1<-ggtree(as.treedata(z), size = 0.1, layout = "daylight", color="gray80")+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))

p2<-ggtree(as.treedata(z), size = 0.4, layout = "daylight", aes(col=ESCO/N, alpha=ESCO/N))+
  scale_color_continuous(high='red', low='gray70', guide="none")+
  scale_alpha_continuous(range = c(0,1), guide="none")+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))

p3<-ggtree(as.treedata(z), size = 0.1, layout = "daylight", color="gray80", alpha=0)+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))+
  geom_tiplab(size=0.2, alpha=1)+
  geom_nodelab(size=0.2, col="black")

p1
p2
p3
dev.off()
