library(ggplot2)
library(data.table)

options(echo=TRUE)
options(scipen=99999)
args <- commandArgs(trailingOnly = TRUE)
`%!in%` <- Negate(`%in%`)
x <- fread(file.path(args[1]), sep = '\t', header = T)
x$query.name <- gsub(":.*","",x$query.name)
x$query.name <- gsub("#J.*","",x$query.name)
x$ref.name<- gsub(":.*","",x$query.name)

p<-ggplot(x, aes(x=query.start, xend=query.end, y=ref.start,yend=ref.end)) + 
            geom_segment(size=0.3) +
            facet_wrap(. ~ query.name, ncol=10,nrow=10, scale="free") +
    theme(
      axis.text.x = element_text(angle = 45)
    )  +
    xlab("Query start") +
    ylab("Reference start")

y<-fread(file.path(args[2]), sep = '\t', header = F)

fill_<-c("#EE2222", "#2222EE", "#46FF33", "#FFFC33", "#812778", "#A8A4A4")

for (i in c(1:nrow(y))) {
  
  p <- p +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = y$V2[i], ymax = y$V3[i], fill = fill_[i], alpha = .2, color = "#444444", size = 0.1) + 
    annotate("text", x = 1000, y =(y$V3[i]+y$V2[i])/2, label = y$V4[i], size=1)
    
}

ggsave(file.path(args[3], paste0("grch38.dotplot.pdf")), width = 25, height = 20, bg = "transparent")
