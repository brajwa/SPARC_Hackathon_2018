#### keywords: spatial statistics, SPARC, Powley, K-function, Kashif

#####################################################################################
#####################################################################################
####### needed libraries

require(tibble)
require(magrittr)
require(dplyr)
require(multcomp)
require(emmeans)
require(readxl)
library(httr)
require(ggfortify)


#####################################################################################
#### Get the data from Dropbox

url1 <- "https://www.dropbox.com/s/dumktnih865k88o/Reorganized%20IMA%20Data%20-%20DMJ.xlsx?dl=1"
GET(url1, write_disk(tf <- tempfile(fileext = ".xlsx")))

spatial.data.long <- read_excel(tf, sheet=1)
spatial.data.circ <- read_excel(tf, sheet=2)

##spatial.data.long <- read_excel("/Users/rajwa/Documents/Revolution/Powley/Reorganized IMA Data - DMJ.xlsx", sheet=1)
##spatial.data.circ <- read_excel("/Users/rajwa/Documents/Revolution/Powley/Reorganized IMA Data - DMJ.xlsx", sheet=2)

spatial.data.long <- add_column(spatial.data.long, `layer`=factor("lon", levels=c("lon", "circ")))
spatial.data.circ <- add_column(spatial.data.circ, `layer`=factor("circ", levels=c("lon", "circ")))
spatial.data <- dplyr::bind_rows( spatial.data.long, spatial.data.circ )

#####################################################################################
#### Take a subset of columns which are relevant for our work

spatial.data.s <- (spatial.data[,c(5,6,11:16,18:20)])
spatial.data.s <- spatial.data.s %>% filter(complete.cases(.))
spatial.data.sub <- (spatial.data[,c(11:16,18:19)])
spatial.data.sub <- add_column(spatial.data.sub, `Layer`=spatial.data$'layer', `Region`=factor(spatial.data$'Region'))
spatial.data.sub <- spatial.data.sub %>% filter(complete.cases(.))

#####################################################################################
#### Transform for variance stabilization 

spatial.data.tr <- spatial.data.sub %>% mutate_at(c(1:6), funs(log))
spatial.data.tr <- spatial.data.tr %>% mutate_at(c(7:8), funs(sqrt))

#####################################################################################
#### See if the regions overlap in the PC dimentions
plot(prcomp(as.matrix(spatial.data.tr[,-c(9,10)]))$x[,1:2], col=spatial.data.tr$'Region', pch=19)

#####################################################################################
#### A nicer plot can be done using ggplot engine

#install.packages("ggfortify")
#require(ggfortify)

df <- subset(spatial.data.tr[,-c(9,10)], spatial.data.tr$Layer=="circ")
autoplot(prcomp(df), data = subset(spatial.data.tr, Layer=="circ"), colour = 'Region', loadings = TRUE, loadings.colour = 'blue', frame.type = 'norm',
         loadings.label = TRUE, loadings.label.size = 4) + 
  theme(legend.position="bottom", text = element_text(size=20, color = "black"), axis.text.x = element_text(size = 20, color="black"), axis.text.y = element_text(size = 20, color="black"),
        axis.title=element_text(size=20, color = "black", face="bold"), aspect.ratio=1, plot.margin = margin(t = 5, r = 20, b = 5, l = 5, unit = "pt"))
#  axis.text.x = blue.bold.italic.16.text


cl.names <- c("Area.of.innervation", "Length.of.arbor", "Width.of.arbor", 
              "Total.tree.length", "Highest.branch.order", "Average.branch.length", 
              "Total.no.branches", "No.of.nodes", "Layer", 'Region')

#####################################################################################
#### Convert data from tibble to data.frame (for simplicity)

spatial.data.df <- as.data.frame(spatial.data.sub)
colnames(spatial.data.df) <- cl.names

#####################################################################################
#### Use polynomial contrasts
#### (or other contrasts of your choice inf you are know what you are doing)

contrasts(spatial.data.df$Layer) <- contr.poly(2)
contrasts(spatial.data.df$Region) <- contr.poly(3)

#####################################################################################
#### Just a bunch of different ways of formatting the axis, choose one you like. It is all about aestethics!

library(scales)

scientific_10 <-
  function(x) {
    ifelse(x == 0, "0", parse(text = gsub(
      "[+]", "", gsub("e", " %*% 10^", scientific_format()(x))
    )))
  }

scientific_10_1 <-
  function(x) {
    ifelse(x == 0, "0", parse(text = gsub("e", " %*% 10^", scientific_format()(x))))
  }

scientific_10_2 <-
  function(x) {
    ifelse(x == 0, "0", (
      ifelse(abs(log10(abs(x))) <= 2, x, parse(text = gsub(
        "[+]", "", gsub("e", " %*% 10^", scientific_format()(x))
      )))))
  }

scientific_10_3 <-
  function(x) {
    ifelse(x==0, "0", 
           ifelse(abs(log10(abs(x))) <= 2, x, sfsmisc::pretty10exp(x, sub10 = T, drop.1 = F))
    )
  }


#####################################################################################
#### Go through the features and compute relevant statistics

i=1  #### pick i to get other features (numbered)

x_lab <- sapply(colnames(spatial.data.sub)[i], function(x){gsub("\\(µm)", "[µm]", x)})
gg_plot <- ggplot(data=spatial.data.df, aes( (spatial.data.df[[i]]) )) + 
  geom_histogram(aes(y =..density..), 
                 col="black", 
                 fill="black",
                 bins=20,
                 alpha = 1)  + 
  geom_density(col=2, size=2, fill=2, alpha=.2) + 
  theme(text = element_text(size=20, color = "black"), axis.text.x = element_text(size = 20, color="black"), axis.text.y = element_text(size = 20, color="black"),
        axis.title=element_text(size=20, color = "black", face="bold"), aspect.ratio=1, plot.margin = margin(t = 5, r = 20, b = 5, l = 5, unit = "pt"))+ 
  # uncomment the appropriate formatting of axes
  #  scale_y_continuous(labels = scales::trans_format("log10", math_format(10^.x))) +
  #  scale_x_sqrt(paste("",x_lab, sep=""),
  #               breaks = trans_breaks("sqrt", function(x) x^2),
  #               labels = trans_format("sqrt", math_format(.x^2)))+
  scale_x_log10(paste("",x_lab, sep=""),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  # scale_x_continuous(paste("\n",x_lab, sep=""), labels = scientific_10) +
  scale_y_continuous("Density", labels = scientific_10_3)

gg_plot  ### print the output


#####################################################################################
#### Compute the statistics

tmp.var <- colnames(spatial.data.df)[1:8]
#tmp.var  <- sapply(tmp.var, function(x){paste("log(",x,")", sep="")})
fit <- fit.em <- fit.cont <- list()
gg_plot_f <- list()

spatial.data.df$int <- with(spatial.data.df, interaction(Region, Layer))
contrasts(spatial.data.df$int) <- contr.poly

for (i in 1:6) {
  fit[[i]] <-    lm(log(get(tmp.var[i])) ~ int - 1, data = spatial.data.df)    ##### for some feature the transformation should be "sqrt"
  fit[[i]] <-    update(ref_grid(fit[[i]]), tran = "log")                      ##### same as above
  fit.em[[i]] <- emmeans(fit[[i]], c("int"))
  summary(fit.em[[i]], type = "response")
  
  x_lab <-
    sapply(colnames(spatial.data.sub)[i], function(x) {
      gsub("\\(µm)", "[µm]", x)
    })
  
  par(mar = c(6, 6, 1, 1),
      oma = c(0, 0, 0, 0),
      pty = "s")
  gg_plot_f[[i]] <-
    plot(fit.em[[i]], comparisons = TRUE, type = "response") +
    theme(
      text = element_text(size = 20, color = "black"),
      axis.text.x = element_text(
        size = 20,
        face = "bold",
        color = "black"
      ),
      axis.title = element_text(
        size = 30,
        color = "black",
        face = "bold"
      )
    ) +
    #labs(x = paste("\nlog(", gsub("[.]", "-", tmp.var[[i]]), ")", sep=""), y="Stomach region\n") +
    labs(x = x_lab, y = "Stomach region\n") +
    scale_y_discrete(
      labels = c(
        "Antrum\nlongitudinal layer",
        "Corpus\nlongitudinal layer",
        "Forestomach\nlongitudinal layer",
        "Antrum\ncircular layer",
        "Corpus\ncircular layer",
        "Forestomach\ncircular layer"
      )
    )
}


#####################################################################################
#### Plot one of the features

gg_plot_f[[6]]


#####################################################################################
#### Compute the statistics and plot the barplots

i=1
fit.cont[[i]] <- contrast(fit.em[[i]], method="eff", type="response")

yy <- as.data.frame(fit.cont[[i]])
yy <- cbind(yy[,1], yy)

yy[,1] <- c("Antrum","Corpus", "Forestomach", "Antrum","Corpus", "Forestomach") 
yy[,2] <- c("Longitudinal","Longitudinal", "Longitudinal", 
            "Circular","Circular", "Circular")

zz <- as.data.frame(summary(fit.em[[i]], type="response"))
zz <- cbind(zz[,1], zz)
zz[,1] <- c("Antrum","Corpus", "Forestomach", "Antrum","Corpus", "Forestomach") 
zz[,2] <- c("Longitudinal","Longitudinal", "Longitudinal", 
            "Circular","Circular", "Circular")

colnames(yy)[1:3] <- c("Segment", "Layer", "ES")
colnames(zz)[1:3] <- c("Segment", "Layer", "Mean")

print(yy)
print(zz)

#####################################################################################
#### Create a barplot

i=3 ### feature number


fit <- lm(log(get(tmp.var[i])) ~ Region, data=subset(spatial.data.df, Layer="circ"))    ##### for some feature the transformation should be "sqrt"
fit <- update(ref_grid(fit, tran = "log"))                                              ##### same as above
fit.em <- emmeans(fit, "Region")
fit.em.sum <- summary(fit.em, type = "response")
fit.em.tb <- data.frame(Region=fit.em.sum$Region, means=fit.em.sum$response, lower.CL=fit.em.sum$lower.CL,upper.CL=fit.em.sum$upper.CL)

ggplot(fit.em.tb, aes(
  x = Region,
  y = means,
  fill = as.factor(Region)
)) +
  geom_bar(stat = "identity",
           alpha = 0.7,
           width = 0.7) +
  scale_fill_brewer(guide = FALSE, palette = "Set1") +
  geom_errorbar(
    aes(x = Region, ymin = lower.CL * 1, ymax = upper.CL * 1),
    width = 0.2,
    colour = "darkblue",
    alpha = 0.9,
    size = 2
  ) +
  theme(
    text = element_text(size = 20, color = "black"),
    axis.text.x = element_text(size = 20, face = "bold", color = "black"),
    axis.title = element_text(size = 30, color = "black", face = "bold")
  ) +
  labs(x = "\nStomach region", y = paste(gsub("[.]", " ", tmp.var), "\n")) #+


#####################################################################################
#####################################################################################  
#####################################################################################
##### spatial stats

library(png)
library(spatstat)


stom.win <- content(GET("https://www.dropbox.com/s/u1nzn4ztwqa41am/Stomach03.png?dl=1"))
stom.win <- readPNG(stom.win)
stom.win.mask <- matrix(ncol=ncol(stom.win[,,1]), nrow=nrow(stom.win[,,1]), data=as.logical(stom.win[,,1]), byrow=F)

stom.owin <- owin(mask=stom.win.mask)
stom.poly <- as.polygonal(stom.owin)
plot(stom.poly)
spatial.data.xy <- data.frame(X=round(spatial.data.s[,1]*930.05/100+25,0), Y=round(spatial.data.s[,2]*700.05/100+82,0))


i=3  #### feature number (a column where the feature is located)

##################################################################################### 
#### create point-patterns

sp.dat.circ <- subset(spatial.data.s, layer=="circ")
sp.xy.circ <- subset(spatial.data.xy, spatial.data.s$layer=="circ")
spat.test.circ <- ppp(x=sp.xy.circ[,1], y=sp.xy.circ[,2], marks=log(sp.dat.circ[[i]]), window=stom.poly, checkdup=F)
sp.dat.lon <- subset(spatial.data.s, layer=="lon")
sp.xy.lon <- subset(spatial.data.xy, spatial.data.s$layer=="lon")
spat.test.lon <- ppp(x=sp.xy.lon[,1], y=sp.xy.lon[,2], marks=log(sp.dat.lon[[i]]), window=stom.poly, checkdup=F)


### plot point-patterns
plot(spat.test.lon)
plot(spat.test.circ)


##################################################################################### 
#### generate discretely marked patterns

breaks.hist <-
  quantile(log(spatial.data.s[[i]]), probs = c(0, 0.33, 0.66, 1))   #### generate tertiles
spat.test.cut.lon <-
  cut.ppp(spat.test.lon,
          breaks = (breaks.hist),
          include.lowest = T)
spat.test.cut.circ <-
  cut.ppp(spat.test.circ,
          breaks = (breaks.hist),
          include.lowest = T)

levels(spat.test.cut.lon$marks) <- c("low", "medium", "high")
levels(spat.test.cut.circ$marks) <- c("low", "medium", "high")

### plot point-patterns
plot(spat.test.cut.lon)
plot(spat.test.cut.circ)

### plot point-patterns
par(mar = c(0, 0, 3, 0),
    oma = c(0, 0, 0, 0),
    pty = "m")
plot(
  spat.test.cut.circ,
  bg = c(2, 3, 4),
  pch = 21,
  main = "Circular layer",
  main.cex = 2,
  cex = 2
)
plot(
  spat.test.cut.lon,
  bg = c(4, 3, 2),
  pch = 21,
  main = "Longitudinal layer",
  main.cex = 2,
  cex = 2
)

### plot density maps
par(mar = c(0, 0, 0, 0),
    oma = c(0, 0, 0, 1),
    pty = "m")
plot(
  density(split(spat.test.cut.lon)),
  main = "",
  useRaster = T,
  log = F,
  col = colorRampPalette(c("blue", "yellow", "red"))(120),
  cex.main = 2
)
plot(
  density(split(spat.test.cut.circ)),
  main = "",
  useRaster = T,
  log = F,
  col = colorRampPalette(c("blue", "yellow", "red"))(120),
  cex.main = 2
)

