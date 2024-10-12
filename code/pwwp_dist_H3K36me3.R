library(bio3d)

cif.6OQM.A <- read.cif("files/structures/6OQM.A.cif") #aligned MSH6 6OQM to 5CIU
pdb.5ciu <- read.pdb("files/structures/5ciu.pdb")
pdb.5ciu_D <- trim.pdb(pdb.5ciu, chain="D")
combined_pdb <- cat.pdb(cif.6OQM.A, pdb.5ciu_D)

write.pdb(combined_pdb, file="files/structures/PWWP_H3K36me3_combined_structure.pdb")

atoms<-data.table(combined_pdb$atom)[chain=="A"]

dist_matrix <- dist(atoms[,.(x,y,z)])
dist_dt<-data.table(as.matrix(dist_matrix))
colnames(dist_dt)<-as.character(atoms$eleno)
dist_dt$eleno<-atoms$eleno
dist_dt_melt<-data.table(melt(dist_dt, id.vars = "eleno"))
dist_dt_melt$resid<-atoms$resid[match(dist_dt_melt$eleno, atoms$eleno)]
dist_dt_melt$resid2<-atoms$resid[match(dist_dt_melt$variable, atoms$eleno)]
dist_dt_melt$resno<-atoms$resno[match(dist_dt_melt$eleno, atoms$eleno)]
dist_dt_melt$resno2<-atoms$resno[match(dist_dt_melt$variable, atoms$eleno)]

contacts<-dist_dt_melt[abs(resno-resno2)>1]
contacts_POS<-contacts[,.(min_dist_contacts=min(value), contacts=sum(value<10)), by=.(POS_pdb=resno, resid, MSH6_POS=resno)]
contacts_POS$POS=contacts_POS$POS_pdb-84
contacts_POS_conservation_contacts<-contacts_POS

#check
#paste(contacts_POS$POS, contacts_POS$resid)
#paste(pwwp_conservation$POS, pwwp_conservation$consensus)


atoms<-data.table(combined_pdb$atom)
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
contacts<-contacts[chain %in% c("A")]

contacts_POS<-contacts[chain2 %in% c("B"),.(mean_dist=mean(value, na.rm=T), min_dist=min(value, na.rm=T), contacts=sum(value<10)), by=.(POS_pdb=resno, resid)]
contacts_POS$POS=contacts_POS$POS_pdb-84

contacts_POS_conservation<-merge(pwwp_conservation, contacts_POS, by.y="POS", by.x="POS", all.x=T)
paste(contacts_POS_conservation$consensus, contacts_POS_conservation$resid)
contacts_POS_conservation$POS
contacts_POS_conservation$POS_AA<-paste(contacts_POS_conservation$POS, contacts_POS_conservation$consensus)
contacts_POS_conservation$MSH6_POS_AA<-paste(contacts_POS_conservation$MSH6_POS, contacts_POS_conservation$consensus)


contacts_POS_conservation$other_contacts<-contacts_POS_conservation_contacts$contacts[match(contacts_POS_conservation$POS, contacts_POS_conservation_contacts$POS)]

contacts_POS_conservation$other_contacts_min_dist<-contacts_POS_conservation_contacts$min_dist_contacts[match(contacts_POS_conservation$POS, contacts_POS_conservation_contacts$POS)]

paste(contacts_POS_conservation$resid, contacts_POS_conservation$consensus)

model<-lm(entropy~min_dist+other_contacts_min_dist, contacts_POS_conservation)
summary(model)

contacts_POS_conservation$Property<-factor(contacts_POS_conservation$Property)
levels(contacts_POS_conservation$Property)
shapes <- c("Aro-\nmatic" = 7, "Neg" = 6, "Non-polar\naliphatic" = 21, "Polar\nuncharged" = 14, "Pos" = 3)


# prepare pdb for distance plotting ---------------------------------------
combined_pdb<-read.pdb("files/structures/PWWP_H3K36me3_combined_structure.pdb", ATOM.only=T)
contacts_POS_conservation$ID<-paste(contacts_POS_conservation$POS_pdb, "A")
combined_pdb$atom$ID<-paste(combined_pdb$atom$resno, combined_pdb$atom$chain)
combined_pdb$atom$b<-contacts_POS_conservation$entropy[match(combined_pdb$atom$ID, contacts_POS_conservation$ID)]
combined_pdb$atom$b[is.na(combined_pdb$atom$b)]<-0

write.pdb(combined_pdb, file="files/structures/PWWP_H3K36me3_combined_structure_annot.pdb")

combined_pdb<-read.pdb("files/structures/PWWP_H3K36me3_combined_structure.pdb", ATOM.only=T)
contacts_POS_conservation$ID<-paste(contacts_POS_conservation$POS_pdb, "A")
combined_pdb$atom$ID<-paste(combined_pdb$atom$resno, combined_pdb$atom$chain)
combined_pdb$atom$b<-contacts_POS_conservation$min_dist[match(combined_pdb$atom$ID, contacts_POS_conservation$ID)]
combined_pdb$atom$b[is.na(combined_pdb$atom$b)]<-0

write.pdb(combined_pdb, file="files/structures/PWWP_H3K36me3_combined_structure_annot_dist.pdb")

fwrite(contacts_POS_conservation[,.(POS, consensus, entropy, min_dist)], "tables/PWWP_contacts_POS_conservation.csv")

pdf("~/Documents/shapes_legend.pdf", width=1.5, height=1.5)

dt_shapes<-data.table(shapes=names(shapes),pos=1:5)
ggplot(dt_shapes, aes(x=pos, y=1, shape=shapes))+
  geom_point()+
  scale_shape_manual(values = shapes)+
  theme_void(base_size = 6)
dev.off()
pdf("~/Documents/PWWP_dist_constraint.pdf", width=1.5, height=1.5)
ggplot(contacts_POS_conservation, aes(y=entropy, x=min_dist, shape=Property, group=1))+
  geom_smooth(method="lm", col="orange")+
  geom_point()+
  scale_shape_manual(values = shapes) +
  #geom_point(data=contacts_POS_conservation, col="orange2", size=3, shape=21)+
  #geom_label_repel(data=contacts_POS_conservation[cage==TRUE], aes(label=consensus), col="green4", size=2)+
  theme_classic(base_size = 6)+
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    plot.background = element_blank())

dev.off()

cor.test(contacts_POS_conservation$entropy, contacts_POS_conservation$min_dist)
pdf("~/Documents/pwwp_histone_dist_AA.pdf", width=4, height=.75)
highcol="orange"
p <- ggplot(contacts_POS_conservation[!is.na(min_dist)], aes(x=POS, y=min_dist, col=scale(Property_cons), fill=scale(Property_cons), label=consensus)) +
  geom_text(aes(y=-2), size=1) +
  geom_bar(stat="identity", col=NA, width=0.9, position = "identity") +
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

p <- ggplot(contacts_POS_conservation[!is.na(min_dist)], aes(x=POS, y=entropy, col=scale(Property_cons), fill=scale(Property_cons), label=consensus)) +
  geom_point( aes(shape=Property, y=-.2), col="black", size=0.75, stroke=0.1, fill=NA) +
  geom_bar(stat="identity", col=NA, width=0.9, position = "identity") +
  scale_color_gradient(low='gray', high=highcol) +
  scale_fill_gradient(low='gray', high=highcol) +
  theme_classic(base_size = 6) +
  scale_shape_manual(values = shapes) +
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

