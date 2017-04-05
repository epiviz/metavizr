library("devtools")
library("msd16s")
library("matrixStats")
setwd("d:\\EpiViz\\github\\metavizr\\")
#install_github("nosson/metagenomeSeq")
load_all("../epivizr")
load_all("./")
age  = pData(msd16s)$Age
ord = order(age)
msd16s = msd16s[,ord]
msd16s = filterData(msd16s,present=100)
fData(msd16s) = fData(msd16s)[,c(3:9, 1)]
pData(msd16s) = pData(msd16s)[,c("Type","Country","AgeFactor","Age")]
levels(pData(msd16s)$AgeFactor) = c("0-5 mo.", "6-11 mo.", "12-17 mo.", "18-23 mo.", "24-59 mo.")
fData(msd16s)[,1] = "Bacteria"
df = fData(msd16s)
fData(msd16s)$OTU = paste(df$OTU, " (", df$genus, ")", sep="")
#mgr = startMetaviz(localURL="http://epiviz-dev.cbcb.umd.edu/metavis/", workspace = "Ey4CYuaTjNd",useDevel=TRUE, debug=TRUE, verbose=TRUE)
#mgr = startMetaviz(localURL="http://localhost/epiviz-dev", workspace = "YGwCd2zrYOs", useDevel=FALSE, debug=TRUE, verbose=FALSE)
mgr = startMetaviz(localURL="http://localhost/epiviz-dev/index-standalone.html", workspace = "YGwCd2zrYOs", useDevel=FALSE, debug=TRUE, verbose=FALSE)

#####
# STACK PLOT

ind = which(pData(msd16s)$Type=="Control" & pData(msd16s)$Country=="Gambia")
obj = metavizTransformSelect(msd16s[,ind],rowSums, control=metavizControl(aggregateAtDepth="genus",n=10))
gates = metavizStack(obj,mgr,control=metavizControl(aggregateFun=colSums,aggregateAtDepth=5))
mgr$visualize("icicle", datasource=gates)

# HEATMAP
obj = msd16s[,which(pData(msd16s)$Type=="Control")]
top15 = metavizTransformSelect(obj,rowVars,control=metavizControl(n=50,aggregateAtDepth="OTU"))
ok = which(!fData(top15)$species=="Other")
heat  = metavizHeatmap(top15[ok,],mgr,control=metavizControl(aggregateAtDepth="OTU"))

#####

obj = msd16s
pData(obj)$AgeStatus = interaction(pData(obj)$AgeFactor, pData(obj)$Type)
#pData(obj)$AgeStatus = interaction(pData(obj)$Type, pData(obj)$AgeFactor)
#obj = msd16s[,which(pData(msd16s)$Type=="Control")]
#obj = obj[,sort(pData(obj)$Type,pData(obj)$AgeFactor)]
obj = obj[,order(pData(obj)$AgeStatus)]
top15 = metavizTransformSelect(obj,rowVars,control=metavizControl(n=50,aggregateAtDepth="OTU"))
ok = which(!fData(top15)$species=="Other")
heat  = metavizHeatmap(top15[ok,],mgr,control=metavizControl(aggregateAtDepth="OTU"))


#####

obj = msd16s[,which(pData(msd16s)$Type=="Control" & pData(msd16s)$AgeFactor=="0-5 mo.")]
top15 = metavizTransformSelect(obj,rowVars,control=metavizControl(n=50,aggregateAtDepth="OTU"))
ok = which(!fData(top15)$species=="Other")
heat  = metavizHeatmap(top15[ok,],mgr,control=metavizControl(aggregateAtDepth="OTU"))

obj = msd16s[,which(pData(msd16s)$Type=="Control" & pData(msd16s)$AgeFactor=="6-11 mo.")]
top15 = metavizTransformSelect(obj,rowVars,control=metavizControl(n=50,aggregateAtDepth="OTU"))
ok = which(!fData(top15)$species=="Other")
heat  = metavizHeatmap(top15[ok,],mgr,control=metavizControl(aggregateAtDepth="OTU"))

obj = msd16s[,which(pData(msd16s)$Type=="Control" & pData(msd16s)$AgeFactor=="12-17 mo.")]
top15 = metavizTransformSelect(obj,rowVars,control=metavizControl(n=50,aggregateAtDepth="OTU"))
ok = which(!fData(top15)$species=="Other")
heat  = metavizHeatmap(top15[ok,],mgr,control=metavizControl(aggregateAtDepth="OTU"))

obj = msd16s[,which(pData(msd16s)$Type=="Control" & pData(msd16s)$AgeFactor=="18-23 mo.")]
top15 = metavizTransformSelect(obj,rowVars,control=metavizControl(n=50,aggregateAtDepth="OTU"))
ok = which(!fData(top15)$species=="Other")
heat  = metavizHeatmap(top15[ok,],mgr,control=metavizControl(aggregateAtDepth="OTU"))

obj = msd16s[,which(pData(msd16s)$Type=="Control" & pData(msd16s)$AgeFactor=="24-59 mo.")]
top15 = metavizTransformSelect(obj,rowVars,control=metavizControl(n=50,aggregateAtDepth="OTU"))
ok = which(!fData(top15)$species=="Other")
heat  = metavizHeatmap(top15[ok,],mgr,control=metavizControl(aggregateAtDepth="OTU"))

#####

#obj = msd16s[,which(pData(msd16s)$AgeFactor=="0-5 mo.")]
e = msd16s
topN = metavizTransformSelect(e,rowVars,control=metavizControl(n=50,aggregateAtDepth="OTU"))
pData(topN)$CountryStatus = interaction(pData(topN)$Country, pData(topN)$Type)
withoutOther = which(!fData(topN)$species=="Other")

ageGroup = which(pData(topN)$AgeFactor=="0-5 mo.")
heat  = metavizHeatmap(topN[withoutOther, ageGroup],mgr,control=metavizControl(aggregateAtDepth="OTU"))

ageGroup = which(pData(topN)$AgeFactor=="6-11 mo.")
heat  = metavizHeatmap(topN[withoutOther, ageGroup],mgr,control=metavizControl(aggregateAtDepth="OTU"))

ageGroup = which(pData(topN)$AgeFactor=="12-17 mo.")
heat  = metavizHeatmap(topN[withoutOther, ageGroup],mgr,control=metavizControl(aggregateAtDepth="OTU"))

ageGroup = which(pData(topN)$AgeFactor=="18-23 mo.")
heat  = metavizHeatmap(topN[withoutOther, ageGroup],mgr,control=metavizControl(aggregateAtDepth="OTU"))

ageGroup = which(pData(topN)$AgeFactor=="24-59 mo.")
heat  = metavizHeatmap(topN[withoutOther, ageGroup],mgr,control=metavizControl(aggregateAtDepth="OTU"))

#####

# e = msd16s
# pData(e)$CountryStatus = interaction(pData(e)$Country, pData(e)$Type)
#
# ageGroup = which(pData(e)$AgeFactor=="0-5 mo.")
# topN = metavizTransformSelect(e[,ageGroup],rowVars,control=metavizControl(n=50,aggregateAtDepth="OTU"))
#
# withoutOther = which(!fData(topN)$species=="Other")
# heat  = metavizHeatmap(topN[withoutOther, ],mgr,control=metavizControl(aggregateAtDepth="OTU"))
#
# ageGroup = which(pData(e)$AgeFactor=="24-59 mo.")
# topN = metavizTransformSelect(e[,ageGroup],rowVars,control=metavizControl(n=50,aggregateAtDepth="OTU"))
# withoutOther = which(!fData(topN)$species=="Other")
# heat  = metavizHeatmap(topN[withoutOther, ],mgr,control=metavizControl(aggregateAtDepth="OTU"))

#t = gates$taxonomyTable()
#####
# LINE DIFFERENCE PLOT
data(mouseData)
time = order(as.Date(pData(mouseData)[,2],format="%Y-%m%-%d"))
mouseData = mouseData[,time]
status = pData(mouseData)[,"status"]
pm = pData(mouseData)[,"mouseID"]
s0 = aggSamp(mouseData[,which(pData(mouseData)$status==0)],"date")
s1 = aggSamp(mouseData[,which(pData(mouseData)$status==1)],"date")
d = MRcounts(s0,norm=TRUE,log=TRUE) - MRcounts(s1,norm=TRUE,log=TRUE)
obj = newMRexperiment(d,featureData=featureData(s0))
metavizLine(obj,mgr=mgr,metavizControl(aggregateAtDepth="class",n=nrow(obj),norm=FALSE,aggregateFun=colMeans,minValue=-1,maxValue=2))

metavizLine(obj,mgr=mgr,metavizControl(aggregateAtDepth="class",n=nrow(obj),norm=FALSE,minValue=-1,maxValue=2,aggregateFun=colMeans,
    valuesAnnotationFuns=list(errMinus=function(t)colMeans(t) - colSds(t), errPlus=function(t)colMeans(t) + colSds(t))
))

metavizLine(obj,mgr=mgr,metavizControl(aggregateAtDepth="class",n=nrow(obj),norm=FALSE,minValue=-1,maxValue=2,aggregateFun=colMeans,
                                       valuesAnnotationFuns=list(errMinus=function(t){mu = colMeans(t);s = colSds(t);s=ifelse(is.na(s),yes=0,no=s);mu + s},
                                                                 errPlus=function(t){mu = colMeans(t);s = colSds(t);s=ifelse(is.na(s),yes=0,no=s);mu - s})
))

#####
# ALTERNATIVE LINE PLOT
pData(mouseData)$diet[which(pData(mouseData)$diet=="BK")] = "LF/PP"
df = pData(mouseData)
pData(mouseData)$label = paste(df$date, " (", df$diet, ")", sep="")
obj = metavizTransformSelect(mouseData[,which(pm=="PM10")],fun=rowSums,control=metavizControl(aggregateAtDepth="class",n=5))
metavizLine(obj,mgr=mgr,metavizControl(aggregateAtDepth="class",n=nrow(obj)))

#####
# MDS (Scatter Plot)
pData(mouseData) = pData(mouseData)[,c("mouseID", "diet", "date", "relativeTime")]
pData(mouseData)$diet[which(pData(mouseData)$diet=="BK")] = "LF/PP"
metavizMDS(mouseData,mgr)

#####
# Heatmap
#library(matrixStats)
obj = msd16s[,which(pData(msd16s)$Type=="Control")]
top15 = metavizTransformSelect(obj,rowVars,control=metavizControl(n=25,aggregateAtDepth="OTU"))
ok = which(!fData(top15)$species=="Other")
heat  = metavizHeatmap(top15[ok,],mgr,control=metavizControl(aggregateAtDepth="OTU"))
#  Replace `FOO` with either, `rowVars` or `rowMads` or `rowRanks` (edited)

#####
mgr$service()
mgr$stopServer()
rm(list=ls())
