library("seqinr")
library("plotly")
library(htmlwidgets)
monkey <- read.fasta(file = "Monkeypox.fasta")
monkeyseq <- monkey[[1]] #Putting sequence in vector

#To read sequences from particular location
monkeyseq[100:200]


#Sliding window to analyse GC content
GC(monkeyseq[1:10000])
GC(monkeyseq[10001:20000])
GC(monkeyseq[20001:30000])

#Instead of typing again and again for above  sliding window we can wite a function and see below for automaticity
starts <- seq(1,length(monkeyseq)-1000, by = 1000)
starts
length(monkeyseq)

n<- length(starts)
for(i in 1:n){
  chunk <- monkeyseq[starts[i]:(starts[i]+999)]
  chunkGC <- GC(chunk)
  print(chunkGC)
}
chunk

##Sliding window plot
starts <- seq(1,length(monkeyseq)-1000, by =1000)
n <- length(starts)
chunkGCs <- numeric(n)
for(i in 1:n){
  chunk <- monkeyseq[starts[i]:(starts[i]+999)]
  chunkGC <- GC(chunk)
  print(chunkGC)
  chunkGCs[i]<- chunkGC
}
data <- data.frame(Nucleotide_Start_Position = starts, GC_Content = chunkGCs)
fig <- plot_ly(data, x = ~Nucleotide_Start_Position, y = ~GC_Content, type = 'scatter', mode = 'lines')
fig
html_file_path <- "GC_content_plot.html"
saveWidget(fig, file = html_file_path, selfcontained = TRUE)


##Sliding window plot using function
slidingwindowplot <- function(windowsize, inputseq)
{
  starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
  n <- length(starts)    # Find the length of the vector "starts"
  chunkGCs <- numeric(n) # Make a vector of the same length as vector "starts", but just containing zeroes
  for (i in 1:n) {
    chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)]
    chunkGC <- GC(chunk)
    print(chunkGC)
    chunkGCs[i] <- chunkGC
  }
  data <- data.frame(Nucleotide_Start_Position = starts, GC_Content = chunkGCs)
  fig <- plot_ly(data, x = ~Nucleotide_Start_Position, y = ~GC_Content, type = 'scatter', mode = 'lines',line=list(color='red'))
  fig
  html_file_path <- "GC_content_plot_function.html"
  saveWidget(fig, file = html_file_path, selfcontained = TRUE)
}
slidingwindowplot(500,monkeyseq)


