library(polymorphology2)

# Function to plot conservation scores
plot_conservation <- function(file_path, highcol) {
  amino_acids <- c("G", "A", "V", "L", "I", "P", "M", "F", "Y", "W", "S", "T", "C", "N", "Q", "K", "R", "H", "D", "E")
  properties <- c(rep("Non-polar\naliphatic", 7), rep("Aro-\nmatic", 3), rep("Polar\nuncharged", 5), rep("Pos", 3), rep("Neg", 2))

  properties_table <- data.frame(AA = amino_acids, Property = properties)
  properties_prop.table<-data.table(prop.table(table(property=properties_table$Property)))

  conservation_scores <- fread(file_path)

  # Filter out rows with "-" in the AA column if present
  conservation_scores <- conservation_scores[consensus != "-"]
  conservation_scores$POS<-1:nrow(conservation_scores)

  AA_freq<-conservation_scores[,-c("refAA","consensus","conservation","entropy","X","-"), with=F]
  AA_freq_melt<-data.table(melt(AA_freq, id.vars = c("POS")))
  AA_freq_melt$value[is.na(AA_freq_melt$value)]<-0
  AA_freq_melt$property<-properties_table$Property[match(AA_freq_melt$variable, properties_table$AA)]
  AA_freq_melt_properties<-AA_freq_melt[,.(sum=sum(value)), by=.(POS, property)]

  property_chi<-rbindlist(lapply(1:max(AA_freq_melt_properties$POS), function(p){
    pos<-AA_freq_melt_properties[POS==p]
    pos$count<-round(pos$sum*100)
    pos$prop<-prop.table(pos$count)
    pos<-merge(pos, properties_prop.table, by="property")
    chi<-chisq.test(pos$count, p=pos$N)
    pos$exp<-chi$expected
    pos$oe<-pos$count/pos$exp
    return(data.table(POS=p,
                      chi.p=chi$p.value,
                      Property=pos$property[which(pos$count==max(pos$count))],
                      Property_cons=pos$prop[which(pos$count==max(pos$count))],
                      Property_oe=pos$oe[which(pos$count==max(pos$count))]))
  }))
  conservation_scores<-merge(conservation_scores, property_chi, by="POS")

  p <- ggplot(conservation_scores, aes(x=POS, y=conservation, col=conservation, fill=conservation, label=consensus)) +
    geom_text(fontface = "bold", aes(y=-0.1), size=1.5) +
    geom_bar(stat="identity", col=NA, width=0.9) +
    scale_color_gradient(low='gray', high=highcol) +
    scale_fill_gradient(low='gray', high=highcol) +
    theme_classic(base_size = 6) +
    theme(axis.line = element_blank(),     # Remove axis lines
          axis.text.x = element_blank(),   # Remove x-axis labels
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none",
          panel.background = element_blank(),
          plot.background = element_blank())+
    scale_y_continuous(name="f(consensus A.A.)")# Remove x-axis ticks

  print(p)

  p <- ggplot(conservation_scores, aes(x=POS, y=-log10(chi.p), col=Property_cons, fill=Property_cons, label=consensus)) +
    geom_text(fontface = "bold", aes(y=-1), size=1.5) +
    geom_bar(stat="identity", col=NA, width=0.9) +
    scale_color_gradient(low='gray', high=highcol) +
    scale_fill_gradient(low='gray', high=highcol) +
    theme_classic(base_size = 6) +
    theme(axis.line = element_blank(),     # Remove axis lines
          axis.text.x = element_blank(),   # Remove x-axis labels
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none",
          panel.background = element_blank(),
          plot.background = element_blank())+
    scale_y_continuous(name="-log10(Property_oe)")# Remove x-axis ticks
  print(p)

  p <- ggplot(conservation_scores, aes(x=POS, y=Property_cons, col=Property_cons, fill=Property_cons, label=consensus)) +
    geom_text(fontface = "bold", aes(y=-0.1), size=1.5) +
    geom_bar(stat="identity", col=NA, width=0.9) +
    scale_color_gradient(low='gray', high=highcol) +
    scale_fill_gradient(low='gray', high=highcol) +
    theme_classic(base_size = 6) +
    theme(axis.line = element_blank(),     # Remove axis lines
          axis.text.x = element_blank(),   # Remove x-axis labels
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none",
          panel.background = element_blank(),
          plot.background = element_blank())+
    scale_y_continuous(name="Property_cons")# Remove x-axis ticks

  print(p)

  p <- ggplot(conservation_scores, aes(x=POS, y=Property_oe, col=Property_cons, fill=Property_oe, label=consensus)) +
    geom_text(fontface = "bold", aes(y=-0.1), size=1.5) +
    geom_bar(stat="identity", col=NA, width=0.9) +
    scale_color_gradient(low='gray', high=highcol) +
    scale_fill_gradient(low='gray', high=highcol) +
    theme_classic(base_size = 6) +
    theme(axis.line = element_blank(),     # Remove axis lines
          axis.text.x = element_blank(),   # Remove x-axis labels
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none",
          panel.background = element_blank(),
          plot.background = element_blank())+
    scale_y_continuous(name="Property_oe")# Remove x-axis ticks

  print(p)


  p <- ggplot(conservation_scores, aes(x=POS, y=entropy, col=entropy, fill=entropy, label=consensus)) +
    #geom_text(fontface = "bold", aes(y=0.1), size=2.5) +
    geom_bar(stat="identity", col=NA) +
    scale_color_gradient(high='gray', low=highcol) +
    scale_fill_gradient(high='gray', low=highcol) +
    theme_classic(base_size = 6) +
    theme(axis.line = element_blank(),     # Remove axis lines
          axis.text.x = element_blank(),   # Remove x-axis labels
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none",
          panel.background = element_blank(),
          plot.background = element_blank())+
    scale_y_reverse(name="Shannon Entropy")# Remove x-axis ticks

  print(p)

  return(conservation_scores)
}

pdf("Figures/msh6_domains_consevation.pdf", width=5, height=1)
tudor_conservation<-plot_conservation(paste0("data/", "Muscle", "_msh6_tudors_conservation.csv"), highcol="green4")
pwwp_conservation<-plot_conservation(file_path = paste0("data/", "Muscle", "_msh6_pwwps_conservation.csv"), highcol="purple")
dev.off()

fwrite(tudor_conservation,"tables/tudor_conservation.csv")
fwrite(pwwp_conservation,"tables/pwwp_conservation.csv")


### see plot_conservation() for making AA_freq_melt
pdf("Figures/msh6_domains_consevation_matrix.pdf", width=5, height=2)
ggplot(AA_freq_melt[variable!="X"], aes(x=POS, y=variable, fill=property, alpha=value))+
  #scale_fill_gradient(low='white', high=highcol)+
  geom_vline(xintercept = 1:max(AA_freq_melt$POS), linewidth=0.1, col="gray")+
  geom_tile()+
  scale_alpha_continuous(range=c(0,1))+
  scale_x_continuous(breaks=1:max(AA_freq_melt$POS), labels=conservation_scores$consensus)+
  theme_classic(base_size = 6)+
  facet_grid(property~., scales="free", space = "free")+
  theme(axis.line = element_blank(),     # Remove axis lines
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        plot.background = element_blank())
dev.off()
