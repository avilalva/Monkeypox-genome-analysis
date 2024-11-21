library("seqinr")
library("plotly")
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
data<- data.frame(starts,chunkGCs)
plot(starts,chunkGCs,type="b",xlab="Nucleotide start position",ylab="GC content")
fig <- plot_ly(data, x = "Nucleotide_start_position", y = "GC content", type = 'scatter', mode = 'lines')
fig
