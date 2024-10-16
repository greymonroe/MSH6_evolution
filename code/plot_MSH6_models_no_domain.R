



plot_protein_models_no_domain<-function(domain_search, predicted_msh6_blast_annotated_sub, interpro){
  interpro$Organism<-predicted_msh6_blast_annotated_sub$Organism[match(interpro$Protein, predicted_msh6_blast_annotated_sub$protein)]
  interpro_sub<-interpro[Protein %in% predicted_msh6_blast_annotated_sub$protein & Analysis!="SUPERFAMILY" & !grepl(domain_search, Signature_Description, ignore.case = T)]
  interpro_sub<-interpro[Analysis %in% c("Pfam")
                         & Protein %in% interpro_sub$Protein
                         & grepl("MutS", Signature_Description)]
  data.table(table(interpro_sub$Analysis, interpro_sub$Signature_Description))[order(N, decreasing = T)]
  interpro_sub$Name<-paste(interpro_sub$Organism, interpro_sub$Protein, sep=" | ")
  table(interpro_sub$Signature_Description)
  best_model<-interpro_sub[,.(I=sum(Signature_Description=="MutS domain I"),
                              II=sum(Signature_Description=="MutS domain II"),
                              III=sum(Signature_Description=="MutS domain III"),
                              IV=sum(Signature_Description=="MutS family domain IV"),
                              V=sum(Signature_Description=="MutS domain V"),
                              Length=max(Stop_Location)), by=.(Organism, Protein)]


  # Count the number of zeros in columns I, II, III, IV, and V
  best_model[, zero_count := rowSums(.SD == 0), .SDcols = c("I", "II", "III", "IV", "V")]
  best_model[, greater_than_one_count := rowSums(.SD > 1), .SDcols = c("I", "II", "III", "IV", "V")]

  # Order by zero_count and Length
  setorder(best_model, Organism, zero_count, Length)

  # Subset to get the top row for each Organism
  result <- best_model[, .SD[1], by = Organism][zero_count==0 &greater_than_one_count==0 ]

  print(result)
  if(nrow(result)>200){
    result<-result[sample(1:nrow(result), 200)]
  }



  interpro_sub_result<-interpro_sub[Protein %in% result$Protein]

  domain_dist<-interpro_sub_result[Signature_Description=="MutS domain I", .(dist=min(Start_Location)), by=.(Organism=Organism)]
  domain_dist_max<-interpro_sub_result[Signature_Description=="MutS domain V", .(dist=max(Stop_Location)), by=.(Organism=Organism)]
  interpro_sub_result$dist<-domain_dist$dist[match(interpro_sub_result$Organism, domain_dist$Organism)]
  interpro_sub_result$dist_max<-domain_dist_max$dist[match(interpro_sub_result$Organism, domain_dist_max$Organism)]

  interpro_sub_result$Signature_Description<-factor(interpro_sub_result$Signature_Description,
                                                    levels=c("Domain", "MutS domain I", "MutS domain II", "MutS domain III", "MutS family domain IV", "MutS domain V"))
  interpro_sub_result<-interpro_sub_result[order(Signature_Description)][(dist_max-dist)<1000]
  domain_plot<-ggplot(interpro_sub_result, aes(x=Start_Location-dist, xend=Stop_Location-dist, y=Name, yend=Name, col=Signature_Description))+
    geom_segment(lineend = "round", size=.2)+
    theme_bw(base_size = 6)+
    theme_void()+
    theme(legend.position = "none")+
    scale_color_manual(values=c(colorRampPalette(c("lightblue","darkblue"))(5)))
  domain_plot

  return(domain_plot)
}

pdf("Figures/MSH6_nodomain_models.pdf", width=1, height=2)
plot_protein_models_no_domain(domain_search = "tudor|pwwp", predicted_msh6_blast_annotated[Organism %in% small_tree[tipcolor=="None"]$Organism], interpro)
dev.off()


