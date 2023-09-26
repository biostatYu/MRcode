library(networkD3)
library(jsonlite)
library(data.table)
URL<-'E:\\Cloud_disk\\program\\Hematology\\Pleiotropy\\AD_BALL\\Figure\\Landscape\\AD_BALL.json'
flare <- fromJSON(URL, simplifyDataFrame = FALSE)
p=radialNetwork(List = flare, fontSize = 12, opacity = 0.9, margin=0)

setwd("E:\\Cloud_disk\\program\\Hematology\\Pleiotropy\\AD_BALL\\Figure\\Landscape")
require(htmlwidgets)
saveWidget(p, file="circus.html")

require(webshot2)
Sys.setenv(CHROMOTE_CHROME = "C:\\Program Files (x86)\\Microsoft\\Edge\\Application\\msedge.exe")
webshot("circus.html", "circus.pdf")