require(lubridate)
require(readxl)
require(ggplot2)
require(ape)
require(phytools)
require(nlme)
require(geiger)
require(ggsignif)
require(scales)
require(cowplot)



require(scales)
setwd("/Users/christiandamsgaard/Dropbox/Projects/FishPhysiology_aquaculture/")
# Import data
df<-read_xlsx("./Data/EPA Ammmonia 2013 c.xlsx")
df<-as.data.frame(df)
df<-subset(x = df,select = c("Species","pH","Temp. (˚C)","Total Ammonia LC50 (mM)","Resp mode"))
colnames(df)<-c("sp","ph","temp","lc50","rm")
df$sp             <-sub("\\ ", "_", df$sp)                      # Separate genus and species by underscore

# Correct typos in marks excel sheet
df$sp[which(df$sp=="Boleophthalmus_boddaerti")]<-"Boleophthalmus_boddarti"
df$sp[which(df$sp=="Channa_Striatus")]<-"Channa_striata"
df$sp[which(df$sp=="Periophthalmus_schlosseri")]<-"Periophthalmodon_schlosseri"
df$sp[which(df$sp=="Morone_saxatilis x chrysops")]<-"Morone_saxatilis"

# Format variables
df$ph<-as.numeric(df$ph)
df$lc50<-as.numeric(df$lc50)
df$temp<-as.numeric(df$temp)

# Import tree
tree<-read.tree(file = "/Users/christiandamsgaard/Dropbox/Projects/Air-breathing review_full/mcc.nexus")

tree<-keep.tip(tree,df$sp)
tree<-ladderize(tree)

prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
prior

tree$node.label<-seq(1,length(tree$node.label),1)

inv.phylo<-MCMCglmm::inverseA(tree,nodes="TIPS",scale=TRUE)

df.w<-df[df$rm=="W",]

glmm0w<-
  MCMCglmm::MCMCglmm(
    log10(lc50)~ph,
    random=~sp,
    family="gaussian",
    ginverse=list(sp=inv.phylo$Ainv),
    prior=prior,
    data=df.w,
    nitt=5000000,
    burnin=1000,
    thin=500)
summary(glmm0w)


int<-summary(glmm0w)$solutions[1,1]
slo<-summary(glmm0w)$solutions[2,1]


df$resid<-log10(df$lc50)-(df$ph*slo+int)

resid<-aggregate(resid~sp,df,median)
resid<-setNames(resid$resid,resid$sp)

rm<-aggregate(rm~sp,df,unique)
rm<-setNames(rm$rm,rm$sp)

#phyloaov
paov<-phylANOVA(tree = tree,rm,resid,nsim = 10000)

phylosig(tree = tree,x = resid,method = "lambda",test = T)


# Set colors
col = to.matrix(rm[tree$tip.label],seq=sort(unique(rm)))[,1]
col[which(col == 0)]<-"#377eb8"
col[which(col == 1)]<-"#e41a1c"
col<-col[tree$tip.label]
df$alpha<-ifelse(df$rm=="W",0.6,1)

# Plot log10 lc50 vs. pH scatter plot
p1w<-
  ggplot(data = df, mapping = aes(x = ph, y = lc50, col = rm))+
  
  scale_y_continuous(trans = "log10")+
  theme_classic()+
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.position = c(0.75, 0.90),
    legend.spacing.y = unit(0.01,"cm"),
    legend.spacing.x = unit(0.001,"cm"),
    legend.key.height = unit(0.01,"cm"))+
  
  geom_abline(slope = summary(glmm0w)$solutions[2,1],intercept = summary(glmm0w)$solutions[1,1])+
  geom_point(shape = 16, size = 1,alpha=df$alpha)+
  scale_color_manual(values = c("#377eb8","#e41a1c"),breaks = c("W","A"),labels = c("Water-breathers","Air-breathers"))+
  
  annotation_logticks(sides = "l") +
  labs(col = "Respiratory mode",
       x = expression("Water pH"), 
       y = expression("LC"[50]*" (mM)"));p1w

## Jitter plot of residuals
rm2<-rm
rm2[which(rm2=="A")]<-"zA"

p2<-
  ggplot(
    data = data.frame(
      resid = resid,
      rm = rm2
    ),
    mapping = aes(x=rm,y=resid,color=rm)
    )+
    theme_classic()+
    theme(legend.position = "none",
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8))+
    geom_jitter(shape =16, width = 0.3)+
    scale_color_manual(values = c("#377eb8","#e41a1c"))+
    scale_x_discrete(labels = c("Water-breathers","Air-breathers"))+
    labs(x = expression("Respiratory mode"),
         y = expression("residual log"[10]*"[LC"[50]*" (mM)]"))+
    geom_signif(comparisons=list(c("zA", "W")),y_position = max(1.08*resid), vjust=0,annotation=ifelse(test = paov$Pf<0.001, no = paste("P = ",paov$Pf,sep=""), yes = paste("P < 0.001")),size = 0.2,textsize = 6/3,color="black") 
p2

pdf("./Figures/Fig X AB.pdf",width = 2.25,height = 4.5,useDingbats = F)
plot_grid(p1w,p2,nrow = 2,ncol =1,align = "v",labels = "AUTO",label_size = 12)
dev.off()

pdf("./Figures/Fig X C.pdf",width = 2.25,height = 4.5,pointsize = 6,useDingbats = F)
plotTree.wBars(
  tree = tree,
  x = resid,
  tip.labels = T,method = "plotTree",
  col = col,
  args.plotTree=list(font = 1,edge.width = 0.5), 
  args.barplot=list(),
  scale = 90)
dev.off()
plotTree.barplot
plotTree.barplot(tree = tree,
                 x = resid,
                 tip.labels = T,args.plotTree = list(cex = 0.3))

pdf("./Figures/Fig X Ca.pdf",width = 1.75,height = 4.5,pointsize = 6,useDingbats = F)
plot(tree,no.margin = T,cex =0.8762706)
dev.off()


barplot(height = resid,horiz = T)



is_tip <- tree$edge[,2] <= length(tree$tip.label)
ordered_tips <- tree$edge[is_tip, 2]

par(mar = c(0,0,0,0))

pdf("./Figures/Fig X Cb.pdf",width = 1.75,height = 6,pointsize = 6)

barplot(height = resid[tree$tip.label[ordered_tips]],
        col = col[tree$tip.label[ordered_tips]],
        names.arg = rep("",length(tree$tip.label)),
        horiz = T,space = .3)
dev.off()







### Pond data ###
df<-
  data.frame(
    depth = c(0,1,2,3,0,1,2,0,1,2,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,0,1,2,0,1,2),
    pco2 = c(5.56847885550076,8.38984147562115,8.46408786036116,27.0884218494468,3.6958037927232,5.18131110146138,5.92284294191895,3.6940579037392,4.43909153306809,5.18012862587324,0.738986611280672,0.739160758544639,0.739160758544639,2.95733610497611,0.737742640180414,1.47762316149568,2.21748227563392,2.21748227563392,0.737011873401243,1.47727132646333,2.21695983384202,2.95664303417856,22.1759171498646,29.160194103323,25.9696385144098,37.8414732638542,14.800367502495,15.5403858776197,19.9804961283682,28.8540949649426,13.3111155253128,14.7901283614586,14.0506219433857,17.0086476156774,14.7620510297083,23.6249968748933,25.8398403319145,25.0893836553658,28.7929649412762,29.5241020594166,19.1813086446908,19.1813086446908,19.185999265868),
    ph = c(6.92,6.88,6.66,6.73,6.79,6.74,6.69,6.9,6.87,6.87,7.44,7.44,7.31,7.05,7.45,7.37,7.31,7.26,7.36,7.3,7.27,7.25,6.37,6.37,6.36,6.18,6.37,6.36,6.37,6.35,6.5,6.49,6.5,6.46,6.31,6.29,6.3,6.28,6.26,6.27,6.38,6.39,6.43),
    po2 = c(59.1075468915221,54.4411616106125,43.55292928849,1.55480119930387,98.7969661478351,74.1230964302063,60.800944220269,86.6773746533367,70.67921538951,64.1840137296948,145.838268749629,120.941113732695,109.172196135147,63.0404247487364,149.147585208674,111.597119866171,99.4163828634955,88.112027802693,124.604017894391,102.904874012278,91.9617108676007,86.563486013542,24.860760630885,17.0991519875592,13.9902152625485,3.10893672501077,4.96108318683631,4.80604933724768,4.49598163807041,2.32497419044441,31.9148784847735,24.9431817283909,23.858695566287,21.2249434583202,49.7916600206546,46.5556403260104,45.4729510161032,46.9968431247922,45.7822908189338,36.4932663505419,103.707802772042,94.4343777849738,95.5396350365842),
    size = c("s","s","s","s","s","s","s","s","s","s","s","s","s","s","s","s","s","s","s","s","s","s","l","l","l","l","l","l","l","l","l","l","l","l","l","l","l","l","l","l","l","l","l")
  )

df.pco2<-
  cbind(
    aggregate(pco2 ~ depth+size,df,mean),
    aggregate(pco2 ~ depth+size,df,function(x){sd(x)/sqrt(length(x))})[,3])
colnames(df.pco2)[2]<-"Fish size"
colnames(df.pco2)[3]<-"mean"
colnames(df.pco2)[4]<-"err"
df.pco2$`Fish size`<-as.character(df.pco2$`Fish size`)

p3<-
  ggplot(data = df.pco2, mapping = aes(x = mean, y = depth, shape = `Fish size`))+
  scale_y_reverse(limits=c(3.1,-0.5))+
  theme_classic()+
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 8),
    axis.text = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8,hjust = 1),
    legend.position = c(0.85, 0.95),
    legend.spacing.y = unit(0.02,"cm"),
    legend.spacing.x = unit(0.001,"cm"),
    legend.key.height = unit(0.02,"cm"))+
  geom_point(size = 2)+
  geom_path()+
  geom_errorbarh(mapping = aes(xmin = mean - err, xmax = mean + err),height = 0.1)+
  scale_shape_manual(values = c(16,17),breaks = c("s","l"),labels = c("Small","Large"))+
  labs(col = "Fish size",
       x = expression("PCO"[2]*" (mmHg)"), 
       y = expression("Depth (m)"));p3


df.ph<-
  cbind(
    aggregate(ph ~ depth+size,df,mean),
    aggregate(ph ~ depth+size,df,function(x){sd(x)/sqrt(length(x))})[,3])
colnames(df.ph)[3]<-"mean"
colnames(df.ph)[4]<-"err"

p4<-
  ggplot(data = df.ph, mapping = aes(x = mean, y = depth, shape = size))+
  scale_y_reverse(limits=c(3.1,-0.5))+
  theme_classic()+
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 8),
    axis.text = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.position = "none",#c(0.75, 0.90),
    legend.spacing.y = unit(0.01,"cm"),
    legend.spacing.x = unit(0.001,"cm"),
    legend.key.height = unit(0.01,"cm"))+
  geom_point(size = 2)+
  geom_path()+
  geom_errorbarh(mapping = aes(xmin = mean - err, xmax = mean + err),height = 0.1)+
  scale_shape_manual(values = c(16,17),breaks = c("s","l"),labels = c("Small","Large"))+
  labs(col = "Fish size",
       x = expression("pH"), 
       y = expression("Depth (m)"));p4



df.po2<-
  cbind(
    aggregate(po2 ~ depth+size,df,mean),
    aggregate(po2 ~ depth+size,df,function(x){sd(x)/sqrt(length(x))})[,3])
colnames(df.po2)[3]<-"mean"
colnames(df.po2)[4]<-"err"

p5<-
  ggplot(data = df.po2, mapping = aes(x = mean, y = depth, shape = size))+
  scale_y_reverse(limits=c(3.1,-0.5))+
  theme_classic()+
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 8),
    axis.text = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.position = "none",#c(0.75, 0.90),
    legend.spacing.y = unit(0.01,"cm"),
    legend.spacing.x = unit(0.001,"cm"),
    legend.key.height = unit(0.01,"cm"))+
  geom_point(size = 2)+
  geom_path()+
  geom_errorbarh(mapping = aes(xmin = mean - err, xmax = mean + err),height = 0.1)+
  scale_shape_manual(values = c(16,17),breaks = c("s","l"),labels = c("Small","Large"))+
  labs(col = "Fish size",
       x = expression("PO"[2]*" (mmHg)"), 
       y = expression("Depth (m)"));p5


# Set up figure in matrix
plot<-
  plot_grid(
    p3,p4,p5,
    nrow = 1,
    ncol =3,
    align = 'h',hjust = -2,
    labels = "AUTO",
    label_size = 12
    )

# Save pond data figure PDF
pdf("./Figures/3_Pond data.pdf",width = 4.5,height = 2.25,useDingbats = F)
grid.arrange(arrangeGrob(plot, left = textGrob(label = "Depth (m)",rot = 90, gp=gpar(fontsize=8))))
dev.off()

jpeg("./Figures/3_Pond data.jpeg",width = 4.5,height = 2.25,units = "in",res = 300)
grid.arrange(arrangeGrob(plot, left = textGrob(label = "Depth (m)",rot = 90, gp=gpar(fontsize=8))))
dev.off()




### Production data ###
df<-
  data.frame(
    year = rep(seq(from = 2008, to = 2017, by = 1),3),
    group = c(rep("salmo",10),rep("ab",10),rep("panga",10)),
    prod = c(2.2,2.25,2.27,2.3,2.3,3,3,3.1,3.2,3.3,3.1,3.2,3.3,3.5,3.7,3.95,4.2,4.5,4.9,5.5,1.9,1.95,2.0,2.05,2.1,2.15,2.22,2.3,2.4,2.5)*10^6
    )
df<-
  data.frame(
    year = rep(seq(from = 2008, to = 2017, by = 1),3),
    group = c(rep("salmo",10),rep("ab",10),rep("panga",10)),
    prod = c(2278404.54,2427226.04,2380941.34,2742162.32,3187830.41,3131954.47,3365152.2,3339563.44,3266564.06,3427116.67,3089517.06,3187335.35,3618379.61,3929396.57,4489869.93,4788618.79,5084058.85,5144339.65,5462466.03,5705521.76,1810677,1698889,1931673.28,2050911.75,2364206.23,2408650.67,2493769.67,2469442.89,2645719.49,2811803.16)
    )
p6<-
  ggplot(data = df, mapping = aes(x = year, y = prod, col = group,shape = group))+
  theme_classic()+
  scale_x_continuous(breaks = seq(from = 2009, to = 2017, by = 2))+
  scale_y_continuous(label=function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))} )+
  scale_color_manual(name = "",
    values = c("#e41a1c","#377eb8","#e41a1c"),
    breaks = c("ab","salmo","panga"),
    labels = c("All air-breathers","Salmonidae","Pangasiidae")
    )+
  scale_shape_manual(name = "",
    values = c(16,17,17),
    breaks = c("ab","salmo","panga"),
    labels = c("All air-breathers","Salmonidae","Pangasiidae")
    )+
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.position = c(0.32, 0.95),
    legend.spacing.y = unit(0.01,"cm"),
    legend.spacing.x = unit(0.001,"cm"),
    legend.key.height = unit(0.01,"cm"))+
  geom_point(size = 2)+
  geom_path()+
  #scale_shape_manual(values = c(16,17),breaks = c("s","l"),labels = c("Small","Large"))+
  labs(col = "",
       x = expression("Year"), 
       y = expression("Live mass (metric tonnes)"));p6

pdf("./Figures/1_Production.pdf",width = 2.25,height = 2.25,useDingbats = F)
p6
dev.off()

jpeg("./Figures/1_Production.jpeg",width = 2.25,height = 2.25,units = "in",res = 300)
p6
dev.off()





######################## 
### Sjannies figures ###
######################## 

#### FISH 1 ####
df<-read_xlsx("./Data/fish1.xlsx")
df.o<-df.d<-as.data.frame(df)
df.o$day<-day(df.o$time_o2)
df.o$hour<-hour(df.o$time_o2)
df.o<-df.o[df.o$day>=24,]
df.o<-df.o[df.o$day<=26,]
df.o<-df.o[-which(df.o$day==26&df.o$hour>11),]
df.o<-df.o[-which(df.o$day==24&df.o$hour<12),]
df.d$day<-day(df.d$time_o2)
df.d$hour<-hour(df.d$time_o2)
df.d<-df.d[df.d$day>=24,]
df.d<-df.d[df.d$day<=26,]
df.d<-df.d[-which(df.d$day==26&df.d$hour>11),]
df.d<-df.d[-which(df.d$day==24&df.d$hour<12),]


# Plot depth data
ggplot(data = df.d,mapping = aes(x = time_d,y=depth))+
  scale_y_reverse(limits=c(2.5,0))+
  theme_classic()+
  theme_classic()+
  theme(
    axis.title = element_text(size = 8),
    plot.margin = margin(b = 0,unit = "pt"),
    axis.text.y = element_text(size = 6),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )+  
  geom_point(size = 0.5)+
  labs(x = expression(""), 
       y = expression(""))->p7;p7


# Plot PO2 data
ggplot(data = df.o,mapping = aes(x = time_o2,y=o2))+
  scale_y_continuous(limits=c(0,150))+
  theme_classic()+
  theme(
    plot.margin = margin(t = 0,unit = "pt"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6)
  )+
  geom_point(size = 0.5)+
  labs(x = expression(""), 
       y = expression(""))->p8;p8

# Set up figure in matrix
plot.1<-
  plot_grid(
    p8,p7,
    nrow = 2,
    ncol =1,
    align = 'v',hjust = -2,
    #labels = "AUTO",
    label_size = 12
  )

plot.1


#### FISH 2 ####

df<-read_xlsx("./Data/fish2.xlsx")
df.o<-df.d<-as.data.frame(df)
df.o$day<-day(df.o$time_o2)
df.o$hour<-hour(df.o$time_o2)
df.o<-df.o[df.o$day>=24,]
df.o<-df.o[df.o$day<=26,]
df.o<-df.o[-which(df.o$day==26&df.o$hour>11),]
df.o<-df.o[-which(df.o$day==24&df.o$hour<12),]
df.d$day<-day(df.d$time_o2)
df.d$hour<-hour(df.d$time_o2)
df.d<-df.d[df.d$day>=24,]
df.d<-df.d[df.d$day<=26,]
df.d<-df.d[-which(df.d$day==26&df.d$hour>11),]
df.d<-df.d[-which(df.d$day==24&df.d$hour<12),]


# Plot depth data
ggplot(data = df.d,mapping = aes(x = time_d,y=depth))+
  scale_y_reverse(limits=c(2.5,0))+
  theme_classic()+
  theme_classic()+
  theme(
    axis.title = element_text(size = 8),
    plot.margin = margin(b = 0,unit = "pt"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )+  
  geom_point(size = 0.5)+
  labs(x = expression(""), 
       y = expression(""))->p7;p7


# Plot PO2 data
ggplot(data = df.o,mapping = aes(x = time_o2,y=o2))+
  scale_y_continuous(limits=c(0,150))+
  theme_classic()+
  theme(
    plot.margin = margin(t = 0,unit = "pt"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title = element_text(size = 8)
  )+
  geom_point(size = 0.5)+
  labs(x = expression(""), 
       y = expression(""))->p8;p8


# Set up figure in matrix
plot.2<-
  plot_grid(
    p8,p7,
    nrow = 2,
    ncol =1,
    align = 'v',hjust = -2,
    #labels = "AUTO",
    label_size = 12
  )

plot.2



#### FISH 6 ####
df<-read_xlsx("./Data/fish6.xlsx")
#df<-read_xlsx("./Data/fish10.xlsx")
df.o<-df.d<-as.data.frame(df)
df.o$day<-day(df.o$time_o2)
df.o$hour<-hour(df.o$time_o2)
df.o<-df.o[df.o$day>=8,]
df.o<-df.o[df.o$day<=10,]
df.o<-df.o[-which(df.o$day==10&df.o$hour>11),]
df.o<-df.o[-which(df.o$day==8&df.o$hour<12),]

df.d$day<-day(df.d$time_o2)
df.d$hour<-hour(df.d$time_o2)
df.d<-df.d[df.d$day>=8,]
df.d<-df.d[df.d$day<=10,]
df.d<-df.d[-which(df.d$day==10&df.d$hour>11),]
df.d<-df.d[-which(df.d$day==8&df.d$hour<12),]


# Plot depth data
ggplot(data = df.d,mapping = aes(x = time_d,y=depth))+
  scale_y_reverse(limits=c(2.5,0))+
  theme_classic()+
  theme_classic()+
  theme(
    axis.title = element_text(size = 8),
    plot.margin = margin(b = 0,unit = "pt"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )+  
  geom_point(size = 0.5)+
  labs(x = expression(""), 
       y = expression(""))->p7;p7


# Plot PO2 data
ggplot(data = df.o,mapping = aes(x = time_o2,y=o2))+
  scale_y_continuous(limits=c(0,150))+
  theme_classic()+
  theme(
    plot.margin = margin(t = 0,unit = "pt"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title = element_text(size = 8)
  )+
  geom_point(size = 0.5)+
  labs(x = expression(""), 
       y = expression(""))->p8;p8


# Set up figure in matrix
plot.6<-
  plot_grid(
    p8,p7,
    nrow = 2,
    ncol =1,
    align = 'v',hjust = -2,
    #labels = "AUTO",
    label_size = 12
  )






#### FISH 10 ####
df<-read_xlsx("./Data/fish10.xlsx")
df.o<-df.d<-as.data.frame(df)
df.o$day<-day(df.o$time_o2)
df.o$hour<-hour(df.o$time_o2)
df.o<-df.o[df.o$day>=8,]
df.o<-df.o[df.o$day<=10,]
df.o<-df.o[-which(df.o$day==10&df.o$hour>11),]
df.o<-df.o[-which(df.o$day==8&df.o$hour<12),]

df.d$day<-day(df.d$time_o2)
df.d$hour<-hour(df.d$time_o2)
df.d<-df.d[df.d$day>=8,]
df.d<-df.d[df.d$day<=10,]
df.d<-df.d[-which(df.d$day==10&df.d$hour>11),]
df.d<-df.d[-which(df.d$day==8&df.d$hour<12),]


# Plot depth data
ggplot(data = df.d,mapping = aes(x = time_d,y=depth))+
  scale_y_reverse(limits=c(2.5,0),position = "right")+
  theme_classic()+
  theme_classic()+
  theme(
    axis.title = element_text(size = 8),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    plot.margin = margin(b = 0,unit = "pt"),
    axis.text.y = element_text(size = 6),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )+  
  geom_point(size = 0.5)+
  labs(x = expression(""), 
       y = expression(""))->p7;p7


# Plot PO2 data
ggplot(data = df.o,mapping = aes(x = time_o2,y=o2))+
  scale_y_continuous(limits=c(0,150),position = "right")+
  theme_classic()+
  theme(
    plot.margin = margin(t = 0,unit = "pt"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6)
  )+
  geom_point(size = 0.5)+
  labs(x = expression(""), 
       y = expression(""))->p8;p8


# Set up figure in matrix
plot.10<-
  plot_grid(
    p8,p7,
    nrow = 2,
    ncol =1,
    align = 'v',hjust = -2,
    #labels = "AUTO",
    label_size = 12
  )



pdf("./Figures/5_telemetry.pdf",width = 4,height = 2.25,useDingbats = F)
plot_grid(
  plot.1,plot.2,plot.6,plot.10,
  nrow = 1,
  ncol =4,
  vjust = 1.2,
  hjust = c(-2.25,-1,-.95,0.15),
  labels = "AUTO",
  label_size = 12
)
dev.off()


####################
#### NH3 Growth ####
####################
df<-read_xlsx("./Data/Cong growth ammonia.xlsx",sheet = "Sheet3")
df<-as.data.frame(df)
df$ph<-as.character(df$ph)

26/14

p11<-
  ggplot(data = df[df$ph=="7",], mapping = aes(x = time, y = fcr, col = treat))+
  scale_x_continuous(limits = c(0,92),breaks = c(30,60,90), labels = c("0-30","0-60","0-90"))+
  scale_y_continuous(limits = c(1,3.5))+
  theme_classic()+
  theme(
    axis.title = element_blank(),
    plot.margin = margin(l = 0,r=0,unit = "pt"),
    #axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8,hjust = 1),
    legend.background = element_rect(colour = NA),
    #legend.key = element_rect(colour = "white", fill = NA),
    legend.position = c(0.5, 0.9),
    legend.spacing.y = unit(0.02,"cm"),
    legend.spacing.x = unit(0.001,"cm"),
    legend.key.height = unit(0.02,"cm"))+
  geom_point(size = 2)+
  geom_path()+
  geom_errorbar(mapping = aes(ymin = fcr - fcr.se, ymax = fcr + fcr.se),width = 2)+
  scale_color_manual(values = c("black","#377eb8","#e41a1c"),breaks = c("c","10","26"),labels = c("Control","0.7 mmol/L","1.9 mmol/L"))+
  labs(col = "Water pH 6.5-7.0",
       x = expression("Time (days)"), 
       y = expression("FCR"));p11


p12<-
  ggplot(data = df[df$ph=="8",], mapping = aes(x = time, y = fcr, col = treat))+
  scale_x_continuous(limits = c(0,92),breaks = c(30,60,90), labels = c("0-30","0-60","0-90"))+
  scale_y_continuous(limits = c(1,3.5),position = "right")+
  theme_classic()+
  theme(
    axis.title = element_blank(),
    plot.margin = margin(l = 0,r=0,unit = "pt"),
    axis.text = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8,hjust = 1),
    legend.position = c(0.5, 0.9),
    legend.spacing.y = unit(0.02,"cm"),
    legend.spacing.x = unit(0.001,"cm"),
    legend.key.height = unit(0.02,"cm"))+
  geom_point(size = 2)+
  geom_path()+
  geom_errorbar(mapping = aes(ymin = fcr - fcr.se, ymax = fcr + fcr.se),width = 2)+
  scale_color_manual(values = c("black","#377eb8","#e41a1c"),breaks = c("c","7","10"),labels = c("Control","0.5 mmol/L","0.7 mmol/L"))+
  labs(col = "Water pH 7.5-8.0",
       x = expression("Time (days)"), 
       y = expression("FCR"));p12


p13<-
  ggplot(data = df[df$ph=="7",], mapping = aes(x = time, y = ww, col = treat))+
  scale_x_continuous(limits = c(0,92),breaks = c(0,30,60,90))+
  scale_y_continuous(limits = c(10,40))+
  theme_classic()+
  theme(
    axis.title = element_blank(),
    plot.margin = margin(l = 0,r=0,unit = "pt"),
    axis.text = element_text(size = 6),
    legend.position = "none")+
  geom_point(size = 2)+
  geom_path()+
  geom_errorbar(mapping = aes(ymin = ww - ww.se, ymax = ww + ww.se),width = 2)+
  scale_color_manual(values = c("black","#377eb8","#e41a1c"),breaks = c("c","10","26"),labels = c("Control","10 mg/L","26 mg/L"))+
  labs(col = "Ammonia concentration",
       x = expression("Time (days)"), 
       y = expression("FCR"));p13


p14<-
  ggplot(data = df[df$ph=="8",], mapping = aes(x = time, y = ww, col = treat))+
  scale_x_continuous(limits = c(0,92),breaks = c(0,30,60,90))+
  scale_y_continuous(limits = c(10,40),position = "right")+
  theme_classic()+
  theme(
    plot.margin = margin(l = 0,r=0,unit = "pt"),
    axis.title = element_blank(),
    #axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    legend.position = "none")+
  geom_point(size = 2)+
  geom_path()+
  geom_errorbar(mapping = aes(ymin = ww - ww.se, ymax = ww + ww.se),width = 2)+
  scale_color_manual(values = c("black","#377eb8","#e41a1c"),breaks = c("c","7","10"),labels = c("Control","7 mg/L","10 mg/L"))+
  labs(col = "Ammonia concentration",
       x = expression("Time (days)"), 
       y = expression("FCR"));p14



pdf("./Figures/6_NH3growth.pdf",width = 4,height = 4,useDingbats = F)
plot_grid(
    p11,p12,p13,p14,
    nrow = 2,
    ncol =2,
    align = 'hv'
    #labels = "AUTO",
    #label_size = 12
  )
dev.off()
