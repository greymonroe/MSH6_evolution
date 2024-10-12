

### code for tudor domain pdb update to make general for PWWP too
library(bio3d)
library(data.table)
pdb <- read.pdb("files/structures/AF-O04716-F1-model_v4.pdb")
atoms<-data.table(pdb$atom)
tudor<-atoms[resno>115 & resno < 179]
dist_matrix <- dist(tudor[,.(x,y,z)])
dist_dt<-data.table(as.matrix(dist_matrix))
colnames(dist_dt)<-as.character(tudor$eleno)
dist_dt$eleno<-tudor$eleno
dist_dt_melt<-data.table(melt(dist_dt, id.vars = "eleno"))
dist_dt_melt$resid<-tudor$resid[match(dist_dt_melt$eleno, tudor$eleno)]
dist_dt_melt$resid2<-tudor$resid[match(dist_dt_melt$variable, tudor$eleno)]
dist_dt_melt$resno<-tudor$resno[match(dist_dt_melt$eleno, tudor$eleno)]
dist_dt_melt$resno2<-tudor$resno[match(dist_dt_melt$variable, tudor$eleno)]

contacts<-dist_dt_melt[abs(resno-resno2)>1]
contacts_POS<-contacts[,.(min_dist=min(value), contacts=sum(value<10)), by=.(POS_pdb=resno-117, resid, MSH6_POS=resno)]
contacts_POS_conservation<-merge(tudor_conservation, contacts_POS, by.y="POS_pdb", by.x="POS")
contacts_POS_conservation_contacts<-contacts_POS_conservation


pdb <- read.pdb("files/structures/pdb7de9.ent")

atoms<-data.table(pdb$atom)[resid!="HOH"]
dist_matrix <- dist(atoms[,.(x,y,z)])
dist_dt<-data.table(as.matrix(dist_matrix))
colnames(dist_dt)<-as.character(atoms$eleno)
dist_dt$eleno<-atoms$eleno
dist_dt_melt<-data.table(melt(dist_dt, id.vars = "eleno"))
dist_dt_melt$resid<-atoms$resid[match(dist_dt_melt$eleno, atoms$eleno)]
dist_dt_melt$resid2<-atoms$resid[match(dist_dt_melt$variable, atoms$eleno)]
dist_dt_melt$resno<-atoms$resno[match(dist_dt_melt$eleno, atoms$eleno)]
dist_dt_melt$resno2<-atoms$resno[match(dist_dt_melt$variable, atoms$eleno)]
dist_dt_melt$chain<-atoms$chain[match(dist_dt_melt$eleno, atoms$eleno)]
dist_dt_melt$chain2<-atoms$chain[match(dist_dt_melt$variable, atoms$eleno)]

contacts<-dist_dt_melt[resno!=resno2]
contacts<-contacts[chain=="A"]
contacts_POS<-contacts[chain2=="P",.(mean_dist=mean(value, na.rm=T), min_dist=min(value, na.rm=T), contacts=sum(value<10)), by=.(POS_pdb=resno-600, resid, MSH6_POS=resno)]
tudor_conservation<-fread("tables/tudor_conservation.csv")
contacts_POS_conservation<-merge(tudor_conservation, contacts_POS, by.y="POS_pdb", by.x="POS")
contacts_POS_conservation$POS_AA<-paste(contacts_POS_conservation$POS, contacts_POS_conservation$consensus)
contacts_POS_conservation$MSH6_POS_AA<-paste(contacts_POS_conservation$MSH6_POS, contacts_POS_conservation$consensus)
contacts_POS_conservation$cage<-contacts_POS_conservation$MSH6_POS_AA %in% c("616 W","641 Y","623 Y")

contacts_POS_conservation$other_contacts<-contacts_POS_conservation_contacts$contacts[match(contacts_POS_conservation$POS, contacts_POS_conservation_contacts$POS)]

contacts_POS_conservation$other_contacts_min_dist<-contacts_POS_conservation_contacts$min_dist[match(contacts_POS_conservation$POS, contacts_POS_conservation_contacts$POS)]

model<-lm(entropy~other_contacts+min_dist, contacts_POS_conservation)
anova(model)
contacts_POS_conservation$Property<-factor(contacts_POS_conservation$Property)
fwrite(contacts_POS_conservation[,.(POS, consensus, entropy, min_dist)], "tables/tudor_contacts_POS_conservation.csv")

shapes <- c("Aro-\nmatic" = 7, "Neg" = 6, "Non-polar\naliphatic" = 21, "Polar\nuncharged" = 14, "Pos" = 3)

View(contacts_POS_conservation[,.(refAA, MSH6_POS, POS)])

# prepare pdb for distance plotting ---------------------------------------
combined_pdb<-read.pdb("files/structures/pdb7de9.ent", ATOM.only=T)
contacts_POS_conservation$ID<-paste(contacts_POS_conservation$POS, "A")
combined_pdb$atom$ID<-paste(combined_pdb$atom$resno-600, combined_pdb$atom$chain)
combined_pdb$atom$b<-contacts_POS_conservation$entropy[match(combined_pdb$atom$ID, contacts_POS_conservation$ID)]
combined_pdb$atom$b[is.na(combined_pdb$atom$b)]<-0

write.pdb(combined_pdb, file="files/structures/Tudor_H3K4me1_combined_structure_annot.pdb")


combined_pdb<-read.pdb("files/structures/pdb7de9.ent", ATOM.only=T)
contacts_POS_conservation$ID<-paste(contacts_POS_conservation$POS, "A")
combined_pdb$atom$ID<-paste(combined_pdb$atom$resno-600, combined_pdb$atom$chain)
combined_pdb$atom$b<-contacts_POS_conservation$min_dist[match(combined_pdb$atom$ID, contacts_POS_conservation$ID)]
combined_pdb$atom$b[is.na(combined_pdb$atom$b)]<-0

write.pdb(combined_pdb, file="files/structures/Tudor_H3K4me1_combined_structure_annot_dist.pdb")



pdf("~/Documents/Tudor_dist_constraint.pdf", width=1.5, height=1.5)
ggplot(contacts_POS_conservation, aes(y=entropy, x=min_dist, shape=Property, group=1))+
  geom_smooth(method="lm", col="green4")+
  geom_point()+
  scale_shape_manual(values = shapes) +

  #geom_point(data=contacts_POS_conservation, col="orange2", size=3, shape=21)+
  #geom_label_repel(data=contacts_POS_conservation[cage==TRUE], aes(label=consensus), col="green4", size=2)+
  theme_classic(base_size = 6)+
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    plot.background = element_blank())

ggplot(contacts_POS_conservation, aes(y=entropy, x=min_dist, shape=Property, group=1))+
  geom_smooth(method="lm", col="green4")+
  geom_point()+
  scale_shape_manual(values = shapes) +

  #geom_point(data=contacts_POS_conservation, col="orange2", size=3, shape=21)+
  #geom_label_repel(data=contacts_POS_conservation[cage==TRUE], aes(label=consensus), col="green4", size=2)+
  theme_classic(base_size = 6)+
  theme(
    panel.background = element_blank(),
    plot.background = element_blank())
dev.off()


cor.test(contacts_POS_conservation$min_dist, contacts_POS_conservation$entropy)
pdf("~/Documents/tudor_histone_dist_AA.pdf", width=4, height=.75)
highcol="green4"
p <- ggplot(contacts_POS_conservation, aes(x=POS, y=min_dist, col=scale(Property_cons), fill=scale(Property_cons), label=consensus)) +
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
  scale_y_continuous(name="Min dist")# Remove x-axis ticks

print(p)

p <- ggplot(contacts_POS_conservation, aes(x=POS, y=entropy, col=scale(Property_cons), fill=scale(Property_cons), label=consensus)) +
  #geom_text(fontface = "bold", aes(y=-1), size=1.5) +
  geom_point( aes(shape=Property, y=-.2), col="black", size=0.75, stroke=0.1, fill=NA) +
  scale_shape_manual(values = shapes) +
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
  scale_y_reverse(name="Shannon entropy")# Remove x-axis ticks

print(p)
dev.off()

