
# functions and libraries -------------------------------------------------

library(polymorphology2)
library(ape)
library(phylolm)
library(geiger)
library(caper)
library(cowplot)

reader_specific_phylolm_pgls<-function(trait_data, cols, tree_time){

  traits_time<-data.frame(trait_data)

  tree_time_sub <- drop.tip(tree_time, tree_time$tip.label[!tree_time$tip.label %in% traits_time$Organism])

  phylolm_results<-rbindlist(lapply(cols, function(p){
  message(p)
    trait_data$predictor<-trait_data[[p]]
    traits_time<-data.frame(trait_data[Organism %in% tree_time$tip.label & is.finite(log10(predictor))])
    row.names(traits_time)<-traits_time$Organism
    #traits_time$tip<-as.numeric(traits_time$tip=="None")
    tree_time_sub <- drop.tip(tree_time, tree_time$tip.label[!tree_time$tip.label %in% row.names(traits_time)])
    log10_phylolm<-phylolm(log10(predictor)~factor(tip), traits_time, tree_time_sub, model="lambda")
    log10_phylolm_sum<-summary(log10_phylolm)
    rank_phylolm<-phylolm(rank(predictor)~factor(tip), traits_time, tree_time_sub, model="lambda")
    rank_phylolm_sum<-summary(rank_phylolm)

    log10_lm<-lm(log10(predictor)~factor(tip), traits_time)
    log10_lm_sum<-summary(log10_lm)
    rank_lm<-lm(rank(predictor)~factor(tip), traits_time)
    rank_lm_sum<-summary(rank_lm)


    test<-rep(c("log10_phylolm","rank_phylolm", "log10_lm","rank_lm"), each=2)
    pval=c(
           log10_phylolm_sum$coefficients[-1,4],
           rank_phylolm_sum$coefficients[-1,4],
           log10_lm_sum$coefficients[-1,4],
           rank_lm_sum$coefficients[-1,4])
    B=c(
        log10_phylolm_sum$coefficients[-1,1], rank_phylolm_sum$coefficients[-1,1],
        log10_lm_sum$coefficients[-1,1],rank_lm_sum$coefficients[-1,1])
    readers=gsub("factor\\(tip\\)","",row.names(rank_lm_sum$coefficients)[-1])
    out<-data.table(predictor=p, test, B, pval, group=rep(c("phylolm","phylolm","lm","lm"), each=2), reader=readers)
    return(out)
  }))
  return(phylolm_results)
}

run_pgls<-function(trait_data, cols, tree_time){

  traits_time<-data.frame(trait_data)
  traits_time$tip<-as.numeric(traits_time$tip=="None")

  tree_time_sub <- drop.tip(tree_time, tree_time$tip.label[!tree_time$tip.label %in% traits_time$Organism])


  tree_comp <- comparative.data(tree_time_sub, data.frame(trait_data), Organism, vcv=TRUE, vcv.dim=2)

  pgls_results<-rbindlist(lapply(cols, function(p){
    message(p)
    tree_comp$data$predictor<-tree_comp$data[[p]]
    log10_pgls <- pgls(log10(predictor) ~ tip, tree_comp, lambda='ML')
    log10_pgls_sum<-summary(log10_pgls)
    rank_pgls <- pgls(rank(predictor) ~ tip, tree_comp, lambda='ML')
    rank_pgls_sum<-summary(rank_pgls)

    trait_data$predictor<-trait_data[[p]]
    traits_time<-data.frame(trait_data[Organism %in% tree_time$tip.label & is.finite(log10(predictor))])
    row.names(traits_time)<-traits_time$Organism
    traits_time$tip<-as.numeric(traits_time$tip=="None")
    tree_time_sub <- drop.tip(tree_time, tree_time$tip.label[!tree_time$tip.label %in% row.names(traits_time)])
    log10_phylolm<-phylolm(log10(predictor)~factor(tip), traits_time, tree_time_sub, model="lambda")
    log10_phylolm_sum<-summary(log10_phylolm)
    rank_phylolm<-phylolm(rank(predictor)~factor(tip), traits_time, tree_time_sub, model="lambda")
    rank_phylolm_sum<-summary(rank_phylolm)

    log10_lm<-lm(log10(predictor)~tip, traits_time)
    log10_lm_sum<-summary(log10_lm)
    rank_lm<-lm(rank(predictor)~factor(tip), traits_time)
    rank_lm_sum<-summary(rank_lm)


    test<-c("log10_pgls","rank_pgls","log10_phylolm","rank_phylolm", "log10_lm","rank_lm")
    pval=c(log10_pgls_sum$coefficients[2,4],
           rank_pgls_sum$coefficients[2,4],
           log10_phylolm_sum$coefficients[2,4],
           rank_phylolm_sum$coefficients[2,4],
           log10_lm_sum$coefficients[2,4],
           rank_lm_sum$coefficients[2,4])
    B=c(log10_pgls_sum$coefficients[2,1], rank_pgls_sum$coefficients[2,1],
        log10_phylolm_sum$coefficients[2,1], rank_phylolm_sum$coefficients[2,1],
        log10_lm_sum$coefficients[2,1],rank_lm_sum$coefficients[2,1])
    out<-data.table(predictor=p, test, B, pval, group=c("pgls", "pgls","phylolm","phylolm","lm","lm"))
    return(out)
  }))
  return(pgls_results)
}

phylocompare2<-function(trait_data, predictor="Genome_size", tree_time=tree_time){

  message(predictor)
  #data wrangling
  trait_data$predictor<-trait_data[[predictor]]
  traits_time<-data.frame(trait_data[Organism %in% tree_time$tip.label & is.finite(log10(predictor))])
  row.names(traits_time)<-traits_time$Organism
  traits_time$tip<-as.numeric(traits_time$tip=="None")
  tree_time_sub <- drop.tip(tree_time, tree_time$tip.label[!tree_time$tip.label %in% row.names(traits_time)])

  ### statisticcal tests

  #GLM
  log10_glm<-glm(tip~log10(predictor), traits_time, family="binomial")
  rank_glm<-glm(tip~rank(predictor), traits_time, family="binomial")
  #ape::binaryPGLMM()
  log10_binaryPGLMM <- (binaryPGLMM(tip~(log10(predictor)), data=traits_time, phy=tree_time_sub))
  rank_binaryPGLMM <- (binaryPGLMM(tip~(rank(predictor)), data=traits_time, phy=tree_time_sub))

  #phylolm::phyloglm()
  log10_phyloglm<-phyloglm(traits_time$reader~(log10(traits_time$predictor)), data=traits_time, phy=tree_time_sub, btol = 1000, method="logistic_MPLE", log.alpha.bound=0)
  rank_phyloglm<-phyloglm(traits_time$reader~(rank(traits_time$predictor)), data=traits_time, phy=tree_time_sub, btol = 1000, method="logistic_MPLE", log.alpha.bound=0)

  out<-list(log10_binaryPGLMM=log10_binaryPGLMM,
            rank_binaryPGLMM=rank_binaryPGLMM,
            log10_phyloglm=log10_phyloglm,
            rank_phyloglm=rank_phyloglm,
            log10_glm=log10_glm,
            rank_glm=rank_glm)

  save(out, file = paste0("RData_results/", predictor,"_glm.Rda"))
  return(out)
}

phylocompare2_coefs<-function(out, p){
  log10_glm<-summary(out$log10_glm)$coefficients
  rank_glm<-summary(out$rank_glm)$coefficients

  log10_binaryPGLMM<-out$log10_binaryPGLMM
  rank_binaryPGLMM<-out$rank_binaryPGLMM

  test<-c("log10_glm","rankglm","log10_binaryPGLMM","rank_binaryPGLMM")
  B=c(log10_glm[2,1], rank_glm[2,1],   log10_binaryPGLMM$B[2,1],   rank_binaryPGLMM$B[2,1])

  pval=c(log10_glm[2,4], rank_glm[2,4], log10_binaryPGLMM$B.pvalue[2,1], rank_binaryPGLMM$B.pvalue[2,1])
  row<-data.table(predictor=p,
                  test, B, pval, group=ifelse(grepl("binary", test), "binaryPGLMM","GLM"))
  return(row)
}


# analyses ----------------------------------------------------------------

col_labels<-c(pct_intron_genic="% Intron (gene)",
              pctnongenic="% nonCDS (genome)",
              nongenic= "Total nonCDS bp (genome)",
              Genome_size ="Total bp (genome)",
              mean_exon_number ="# exons/gene",
              mean_intron_length = "Avg intron size",
              Protein_number = "# Protein coding genes",
              gene_lengths = "Gene size (CDS + nonCDS)",
              CDS_genome_total="Total CDS bp (genome)",
              N50="Contig N50",
              mean_total_CDS_length = "Avg. protein length",
              BUSCO="BUSCO",
              gene_percent="% genic",
              mean_CDS_length="Avg CDS exon size")

# tree from timetree
tree_time<-read.tree("files/time_tree_may4.nwk") #from http://www.timetree.org/
tree_time$tip.label<-gsub("_"," ",tree_time$tip.label)
tree_time <- multi2di(tree_time)
tree_time$edge.length <- tree_time$edge.length + 1e-5

trait_data<-fread("tables/S4_species_traits.csv")
cols<-colnames(trait_data)[3:27][-3]

# This takes hours to run,
# glms<-lapply(cols, function(p){
#   phylocompare2(trait_data, p, tree_time)
# })
# names(glms)<-cols

#out<-phylocompare2(trait_data, "mean_intron_length", tree_time)

glm_plot_data<-rbindlist(lapply(cols, function(p){
  message(p)
  if(!file.exists(paste0("RData_results/", p,"_glm.Rda"))) return(NULL)
  load(paste0("RData_results/", p,"_glm.Rda"))
  old<-phylocompare2_coefs(out, p)

  trait_data$predictor<-trait_data[[p]]
  traits_time<-data.frame(trait_data[Organism %in% tree_time$tip.label & is.finite(log10(predictor))])
  row.names(traits_time)<-traits_time$Organism
  traits_time$tip<-as.numeric(traits_time$tip=="None")
  tree_time_sub <- drop.tip(tree_time, tree_time$tip.label[!tree_time$tip.label %in% row.names(traits_time)])

  log10_phyloglm<-phyloglm(traits_time$tip~(log10(traits_time$predictor)), data=traits_time, phy=tree_time_sub, btol = 1000, method="logistic_MPLE", log.alpha.bound=0)
  rank_phyloglm<-phyloglm(traits_time$tip~(rank(traits_time$predictor)), data=traits_time, phy=tree_time_sub, btol = 1000, method="logistic_MPLE", log.alpha.bound=0)

  summary(log10_phyloglm)$coefficients[2,1]
  test<-c("log10_phyloglm","rank_phyloglm")
  B=c(summary(log10_phyloglm)$coefficients[2,1], summary(rank_phyloglm)$coefficients[2,1])

  pval=c(summary(log10_phyloglm)$coefficients[2,4], summary(rank_phyloglm)$coefficients[2,4])
  row<-data.table(predictor=p,
                  test, B, pval, group="phyloglm")
  # print(row)
  # print(old)
  result<-rbind(old, row)
  result$N<-nrow(traits_time)
  return(result)
}))

pgls_plotdata<-run_pgls(trait_data = trait_data, cols=cols, tree_time = tree_time)
glm_plot_data<-rbind(glm_plot_data, pgls_plotdata, fill=T)

glm_plot_data$group2<-ifelse(grepl("rank", glm_plot_data$test),"Rank","log10")
glm_plot_data$group3<-glm_plot_data$predictor %in% cols[1:18]
#glm_plot_data$gini<-importance$MeanDecreaseGini[match(glm_plot_data$predictor, importance$col)]
glm_plot_data_ordering<-glm_plot_data[group2=="Rank" & group=="GLM"][order(pval, decreasing = T)]
glm_plot_data$predictor<-factor(glm_plot_data$predictor, levels=glm_plot_data_ordering$predictor)
glm_plot_data$pval[-log10(glm_plot_data$pval)==Inf]<-1e-300

pdf("Figures/genome_traits_GLMs.pdf", width=2.7, height=3)
ggplot(glm_plot_data[group3==T & predictor %in% names(col_labels)], aes(x=-log10(pval), y=predictor, fill=group2, group=group2))+
  #geom_point(position=position_dodge(width=0.25))+
  geom_bar(stat="identity", position = "dodge", width=0.25)+
  facet_grid(~factor(group, levels=c("GLM","lm","phyloglm","binaryPGLMM", "pgls","phylolm")), scales="free")+
  geom_vline(xintercept = -log10(0.05), linetype="dashed", size=0.25)+
  theme_bw(base_size = 6)+
  scale_y_discrete(labels=col_labels, name="")+
  scale_fill_manual(values=c("black","gray"), name=NULL, labels=c("log10", "rank"))+
  theme(legend.position = "none", legend.key.size = unit(.25,"cm"), panel.background = element_blank(), strip.background = element_blank())

ggplot(glm_plot_data[!group %in% c("phyloglm","pgls") & group3==T & predictor %in% names(col_labels)], aes(x=-log10(pval), y=predictor, fill=group, group=group2, shape=group2))+
  geom_point(position = position_dodge(width=0.5), size=1, stroke=0.25, alpha=0.75)+
  geom_vline(xintercept = -log10(0.05), linetype="dashed", size=0.25)+
  theme_bw(base_size = 6)+
  facet_grid(~!group %in% c("GLM","lm"), scales="free")+
  scale_shape_manual(values=c(21, 22))+
  scale_y_discrete(labels=col_labels, name="")+
  #scale_fill_manual(values=c("black","gray"), name=NULL, labels=c("log10", "rank"))+
  scale_fill_viridis_d()+
  theme(legend.position = "none", legend.key.size = unit(.25,"cm"), panel.background = element_blank(), strip.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
dev.off()

pdf("Figures/organism_traits_GLMs.pdf", width=2.3, height=1.6)
ggplot(glm_plot_data[group3==F], aes(x=-log10(pval), y=factor(predictor, levels=rev(c("Ne", "body_mass", "longevity","snp_u","indel_u"))), fill=group2, group=group2))+
  #geom_point(position=position_dodge(width=0.25))+
  geom_bar(stat="identity", position = "dodge", width=0.25)+
  facet_grid(~factor(group, levels=c("lm", "GLM","phyloglm","binaryPGLMM", "pgls","phylolm")), scales="free")+
  geom_vline(xintercept = -log10(0.05), linetype="dashed",size=0.25)+
  theme_bw(base_size = 6)+
  scale_fill_manual(values=c("black","gray"), name=NULL)+
  scale_y_discrete(labels=c(body_mass= "Body mass", longevity="Max life span", Ne="Pop. size (Ne)",snp_u="SNP mut. rate", indel_u="InDel mut. rate"), name="")+
  theme(legend.position = "none", strip.background = element_blank())

ggplot(glm_plot_data[!group %in% c("phyloglm","pgls") & group3==F], aes(x=-log10(pval), y=factor(predictor, levels=rev(c("Ne", "body_mass", "longevity","snp_u","indel_u"))), fill=group, group=group2, shape=group2))+
  geom_point(position = position_dodge(width=0.5), size=1, stroke=0.25, alpha=0.75)+
  geom_vline(xintercept = -log10(0.05), linetype="dashed", size=0.25)+
  theme_bw(base_size = 6)+
  facet_grid(~!group %in% c("GLM","lm"), scales="free")+
  scale_shape_manual(values=c(21, 22))+
  scale_y_discrete(labels=c(body_mass= "Body mass", longevity="Max life span", Ne="Pop. size (Ne)",snp_u="SNP mut. rate", indel_u="InDel mut. rate"), name="")+  #scale_fill_manual(values=c("black","gray"), name=NULL, labels=c("log10", "rank"))+
  scale_fill_viridis_d()+
  theme(legend.position = "none", legend.key.size = unit(.25,"cm"), plot.background = element_blank(), panel.background = element_blank(), strip.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())

dev.off()




# comparisons boxplots ----------------------------------------------------

pdf("Figures/genome_properties_readers_metaplots_percentile_box.pdf", width=1, height=2.75)
# all_data_table<-rbindlist(lapply(all_results, function(x) x[["all_data_table"]]))
#
# all_data_table<-all_data_table[trait=="all_traits" & predictor %in% all_results_table$col]
# all_data_table$predictor<-factor(all_data_table$predictor, levels=rev(levels(all_results_table[col %in% all_data_table$predictor]$col)))
#
#
# all_data_table_mean<-all_data_table[trait=="all_traits" & predictor %in% all_results_table$col,.(mean=mean(value), se=sd(value)/sqrt(.N)*2), by=.(predictor, tip)]
#
# ggplot(all_data_table_mean, aes(x=mean, y=tip, fill=tip, col=tip))+
#   geom_point(size=0.25, shape=21)+
#   geom_errorbarh(aes(xmin=mean-se, xmax=mean+se), height=0)+
#   scale_color_manual(values=c(None="gray",PWWP="orange", Tudor="green2"), guide="none")+
#   scale_fill_manual(values=c(None="gray",PWWP="orange", Tudor="green2"), guide="none")+
#   theme_bw(base_size=6)+
#   theme(legend.position = "none", panel.grid = element_blank())+
#   facet_grid(predictor~., space="free", scales="free_x")


trait_names<-cols[cols %in% names(col_labels)]
sub<-trait_data[,c("tip", trait_names), with=F]
sub<-melt(sub, id.vars = "tip")[!is.na(value)]
sub[,rel_rank:=rank(value)/.N, by=variable]
sub$variable<-factor(sub$variable, levels=rev(glm_plot_data_ordering$predictor))

ggplot(sub, aes(x=rel_rank*100, y=tip, fill=tip))+
  #geom_jitter(size=0.25, width=0.25, alpha=0.25)+
  geom_boxplot(size=0.25, alpha=0.75, outlier.size = 0, outlier.color = NA, width=0.5)+
  scale_fill_manual(values=c(None="gray",PWWP="orange", Tudor="green2"), guide="none")+
  theme_bw(base_size=6)+
  #ggtitle(real_sci_format(signif(test1$p.value, digits = 2)), subtitle = real_sci_format(signif(fit1$coefficients[2,4], digits = 2)))+
  scale_x_continuous(name="Percentile")+
  theme(legend.position = "none", panel.grid = element_blank(), axis.text.y=element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())+
  facet_grid(variable~., space="free", scales="free_x")

dev.off()

pdf("Figures/organism_properties_readers_metaplots_percentile_box.pdf", width=1, height=1.35)
trait_names<-cols[cols %in% glm_plot_data[group3==F]$predictor]
sub<-trait_data[,c("tip", trait_names), with=F]
sub<-melt(sub, id.vars = "tip")[!is.na(value)]
sub[,rel_rank:=rank(value)/.N, by=variable]

sub$variable<-factor(sub$variable, levels=(c("Ne", "body_mass", "longevity","snp_u","indel_u")))
ggplot(sub, aes(x=rel_rank*100, y=tip, fill=tip))+
  #geom_jitter(size=0.25, width=0.25, alpha=0.25)+
  geom_boxplot(size=0.25, alpha=0.75, outlier.size = 0, outlier.color = NA, width=0.5)+
  scale_fill_manual(values=c(None="gray",PWWP="orange", Tudor="green2"), guide="none")+
  theme_bw(base_size=6)+
  #ggtitle(real_sci_format(signif(test1$p.value, digits = 2)), subtitle = real_sci_format(signif(fit1$coefficients[2,4], digits = 2)))+
  scale_x_continuous(name="Percentile")+
  theme(plot.background = element_blank(), panel.background = element_blank(), legend.position = "none", panel.grid = element_blank(), axis.text.y=element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())+
  facet_grid(variable~., space="free", scales="free_x")
dev.off()


pdf("Figures/genome_organism_properties_readers_metaplots_log10_box.pdf", width=7, height=4)

all_col_labels<-c(col_labels, c(body_mass= "Body mass", longevity="Max life span", Ne="Pop. size (Ne)",snp_u="SNP mut. rate", indel_u="InDel mut. rate"))
all_col_labels<-all_col_labels[names(all_col_labels)%in% colnames(trait_data)]
# trait_names<-cols[cols %in% names(col_labels)]
# sub<-trait_data[,c("tip", trait_names), with=F]
# sub<-melt(sub, id.vars = "tip")[!is.na(value)]
# sub$variable<-factor(sub$variable, levels=rev(glm_plot_data_ordering$predictor))

plotlist<-lapply(c(names(all_col_labels)), function(p){
  message(p)
  sub<-trait_data[,c("tip"), with=F]
  sub$predictor=trait_data[[p]]
  ggplot(sub, aes(x=predictor, y=tip, fill=tip))+
    #geom_jitter(size=0.25, width=0.25, alpha=0.25)+
    geom_boxplot(size=0.25, alpha=0.75, outlier.size = 0, outlier.color = NA, width=0.5)+
    scale_fill_manual(values=c(None="gray",PWWP="orange", Tudor="green2"), guide="none")+
    theme_bw(base_size=6)+
    #ggtitle(real_sci_format(signif(test1$p.value, digits = 2)), subtitle = real_sci_format(signif(fit1$coefficients[2,4], digits = 2)))+
    scale_x_log10(name=paste0("log10(", all_col_labels[p],")"))+
    theme(legend.position = "none", panel.grid = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())+
    ggtitle(all_col_labels[p])

})

plot_grid(plotlist = plotlist, ncol=3)
dev.off()

pdf("Figures/genome_organism_properties_readers_metaplots_log10_box_protostomes.pdf", width=7, height=3.5)

all_col_labels<-c(col_labels, c(body_mass= "Body mass", longevity="Max life span", Ne="Pop. size (Ne)",snp_u="SNP mut. rate", indel_u="InDel mut. rate"))
all_col_labels<-all_col_labels[names(all_col_labels)%in% colnames(trait_data)]
# trait_names<-cols[cols %in% names(col_labels)]
# sub<-trait_data[,c("tip", trait_names), with=F]
# sub<-melt(sub, id.vars = "tip")[!is.na(value)]
# sub$variable<-factor(sub$variable, levels=rev(glm_plot_data_ordering$predictor))

plotlist<-lapply(c(names(all_col_labels)), function(p){
  message(p)
  Prots<-trait_data[Protostomia==T]
  sub<-Prots[,c("tip"), with=F]
  sub$predictor=Prots[[p]]
  sub$reader=sub$tip=="None"
  sub$log10p<-log10(sub$predictor)
  mod<-summary(glm(reader~log10p, data=sub[is.finite(log10p)], family="binomial"))
  pval<-mod$coefficients[2,4]
  pvalannt<-ifelse(pval<1e-8,"***", ifelse(pval<1e-3,"**", ifelse(pval<0.05, "*","")))
  ggplot(sub, aes(x=predictor, y=tip, fill=tip))+
    #geom_jitter(size=0.25, width=0.25, alpha=0.25)+
    geom_boxplot(size=0.25, alpha=0.75, outlier.size = 0, outlier.color = NA, width=0.5)+
    scale_fill_manual(values=c(None="gray",PWWP="orange", Tudor="green2"), guide="none")+
    theme_bw(base_size=6)+
    #ggtitle(real_sci_format(signif(test1$p.value, digits = 2)), subtitle = real_sci_format(signif(fit1$coefficients[2,4], digits = 2)))+
    scale_x_log10(name=paste0("log10(", all_col_labels[p],")"))+
    theme(legend.position = "none", panel.grid = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())+
    ggtitle(paste(all_col_labels[p],pvalannt))

})

plot_grid(plotlist = plotlist, ncol=3)
dev.off()



# phlyolm each reader -----------------------------------------------------

phylolm_plotdata<-reader_specific_phylolm_pgls(trait_data = trait_data, cols=cols, tree_time = tree_time)

phylolm_plotdata$group2<-ifelse(grepl("rank", phylolm_plotdata$test),"Rank","log10")
phylolm_plotdata$group3<-!phylolm_plotdata$predictor %in% cols[1:18]
phylolm_plotdata$pval[-log10(phylolm_plotdata$pval)==Inf]<-1e-300
phylolm_plotdata$predictor<-factor(phylolm_plotdata$predictor, levels=glm_plot_data_ordering$predictor)

pdf("figures/phylolm_each_reader.pdf", width=3, height=3)
ggplot(phylolm_plotdata[group3==F & predictor %in% names(col_labels)], aes(y=predictor, x=-log10(pval), fill=reader, shape=group2))+
  geom_point(position = position_dodge(width=0.5), size=1, stroke=0.25, alpha=0.75)+
  facet_grid(~group, scales="free")+
  geom_vline(xintercept = -log10(0.05), linetype="dashed", size=0.25)+
  theme_bw(base_size = 6)+
  scale_y_discrete(labels=col_labels, name="")+
  theme(legend.position = "none", legend.key.size = unit(.25,"cm"), panel.background = element_blank(), strip.background = element_blank())+
  scale_fill_manual(values=c(None="gray",PWWP="orange", Tudor="green2"), guide="none")+
  scale_shape_manual(values=c(21, 22))

ggplot(phylolm_plotdata, aes(y=predictor, x=-log10(pval), fill=reader, shape=group2))+
  geom_point(position = position_dodge(width=0.5), size=1, stroke=0.25, alpha=0.75)+
  facet_grid(~group, scales="free")+
  geom_vline(xintercept = -log10(0.05), linetype="dashed", size=0.25)+
  theme_bw(base_size = 6)+
  #scale_y_discrete(labels=col_labels, name="")+
  theme(legend.position = "none", legend.key.size = unit(.25,"cm"), panel.background = element_blank(), strip.background = element_blank())+
  scale_fill_manual(values=c(None="gray",PWWP="orange", Tudor="green2"), guide="none")+
  scale_shape_manual(values=c(21, 22))
dev.off()


# change plot sizes -------------------------------------------------------

pdf("Figures/genome_traits_GLMs.pdf", width=5, height=3)
ggplot(glm_plot_data[group3==T & predictor %in% names(col_labels)], aes(x=-log10(pval), y=predictor, fill=group2, group=group2))+
  #geom_point(position=position_dodge(width=0.25))+
  geom_bar(stat="identity", position = "dodge", width=0.25)+
  facet_grid(~factor(group, levels=c("GLM","lm","phyloglm","binaryPGLMM", "pgls","phylolm")), scales="free")+
  geom_vline(xintercept = -log10(0.05), linetype="dashed", size=0.25)+
  theme_bw(base_size = 6)+
  scale_y_discrete(labels=col_labels, name="")+
  scale_fill_manual(values=c("black","gray"), name=NULL, labels=c("log10", "rank"))+
  theme(legend.position = "none", legend.key.size = unit(.25,"cm"), panel.background = element_blank(), strip.background = element_blank())

ggplot(glm_plot_data[group!="phyloglm" & group3==T & predictor %in% names(col_labels)], aes(x=-log10(pval), y=predictor, col=group, group=group2, shape=group2))+
  geom_point()+
  geom_vline(xintercept = -log10(0.05), linetype="dashed", size=0.25)+
  theme_bw(base_size = 6)+
  facet_grid(~group!="GLM", scales="free")+

  scale_y_discrete(labels=col_labels, name="")+
  #scale_fill_manual(values=c("black","gray"), name=NULL, labels=c("log10", "rank"))+
  scale_fill_viridis_d()+
  theme(legend.position = "top", legend.key.size = unit(.25,"cm"), panel.background = element_blank(), strip.background = element_blank())
dev.off()

pdf("Figures/organism_traits_GLMs.pdf", width=4.8, height=1.6)
ggplot(glm_plot_data[group3==F], aes(x=-log10(pval), y=factor(predictor, levels=rev(c("Ne", "body_mass", "longevity","snp_u","indel_u"))), fill=group2, group=group2))+
  #geom_point(position=position_dodge(width=0.25))+
  geom_bar(stat="identity", position = "dodge", width=0.25)+
  facet_grid(~factor(group, levels=c("GLM","lm","phyloglm","binaryPGLMM", "pgls","phylolm")), scales="free")+
  geom_vline(xintercept = -log10(0.05), linetype="dashed",size=0.25)+
  theme_bw(base_size = 6)+
  scale_fill_manual(values=c("black","gray"), name=NULL)+
  scale_y_discrete(labels=c(body_mass= "Body mass", longevity="Max life span", Ne="Pop. size (Ne)",snp_u="SNP mut. rate", indel_u="InDel mut. rate"), name="")+
  theme(legend.position = "none", strip.background = element_blank())
dev.off()




# comparisons boxplots ----------------------------------------------------

pdf("Figures/genome_properties_readers_metaplots_percentile_box.pdf", width=1, height=2.75)
# all_data_table<-rbindlist(lapply(all_results, function(x) x[["all_data_table"]]))
#
# all_data_table<-all_data_table[trait=="all_traits" & predictor %in% all_results_table$col]
# all_data_table$predictor<-factor(all_data_table$predictor, levels=rev(levels(all_results_table[col %in% all_data_table$predictor]$col)))
#
#
# all_data_table_mean<-all_data_table[trait=="all_traits" & predictor %in% all_results_table$col,.(mean=mean(value), se=sd(value)/sqrt(.N)*2), by=.(predictor, tip)]
#
# ggplot(all_data_table_mean, aes(x=mean, y=tip, fill=tip, col=tip))+
#   geom_point(size=0.25, shape=21)+
#   geom_errorbarh(aes(xmin=mean-se, xmax=mean+se), height=0)+
#   scale_color_manual(values=c(None="gray",PWWP="orange", Tudor="green2"), guide="none")+
#   scale_fill_manual(values=c(None="gray",PWWP="orange", Tudor="green2"), guide="none")+
#   theme_bw(base_size=6)+
#   theme(legend.position = "none", panel.grid = element_blank())+
#   facet_grid(predictor~., space="free", scales="free_x")


trait_names<-cols[cols %in% names(col_labels)]
sub<-trait_data[,c("tip", trait_names), with=F]
sub<-melt(sub, id.vars = "tip")[!is.na(value)]
sub[,rel_rank:=rank(value)/.N, by=variable]
sub$variable<-factor(sub$variable, levels=rev(glm_plot_data_ordering$predictor))

ggplot(sub, aes(x=rel_rank*100, y=tip, fill=tip))+
  #geom_jitter(size=0.25, width=0.25, alpha=0.25)+
  geom_boxplot(size=0.25, alpha=0.75, outlier.size = 0, outlier.color = NA, width=0.5)+
  scale_fill_manual(values=c(None="gray",PWWP="orange", Tudor="green2"), guide="none")+
  theme_bw(base_size=6)+
  #ggtitle(real_sci_format(signif(test1$p.value, digits = 2)), subtitle = real_sci_format(signif(fit1$coefficients[2,4], digits = 2)))+
  scale_x_continuous(name="Percentile")+
  theme(legend.position = "none", panel.grid = element_blank(), axis.text.y=element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())+
  facet_grid(variable~., space="free", scales="free_x")

dev.off()

pdf("Figures/organism_properties_readers_metaplots_percentile_box.pdf", width=1, height=1.35)
trait_names<-cols[cols %in% glm_plot_data[group3==F]$predictor]
sub<-trait_data[,c("tip", trait_names), with=F]
sub<-melt(sub, id.vars = "tip")[!is.na(value)]
sub[,rel_rank:=rank(value)/.N, by=variable]

sub$variable<-factor(sub$variable, levels=(c("Ne", "body_mass", "longevity","snp_u","indel_u")))
ggplot(sub, aes(x=rel_rank*100, y=tip, fill=tip))+
  #geom_jitter(size=0.25, width=0.25, alpha=0.25)+
  geom_boxplot(size=0.25, alpha=0.75, outlier.size = 0, outlier.color = NA, width=0.5)+
  scale_fill_manual(values=c(None="gray",PWWP="orange", Tudor="green2"), guide="none")+
  theme_bw(base_size=6)+
  #ggtitle(real_sci_format(signif(test1$p.value, digits = 2)), subtitle = real_sci_format(signif(fit1$coefficients[2,4], digits = 2)))+
  scale_x_continuous(name="Percentile")+
  theme(plot.background = element_blank(), panel.background = element_blank(), legend.position = "none", panel.grid = element_blank(), axis.text.y=element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())+
  facet_grid(variable~., space="free", scales="free_x")
dev.off()

