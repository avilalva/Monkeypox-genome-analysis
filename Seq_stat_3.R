library("seqinr")
library("Biostrings")

# Function to find potential start and stop codons
findPotentialstartandstop <- function(sequence) {
  codons <- c("atg", "taa", "tag", "tga")
  
  positions <- vector()  # Initialize positions
  types <- vector()      # Initialize types
  
  # Convert sequence to lowercase to handle case sensitivity
  sequence <- tolower(sequence)
  
  # Ensure the sequence is in the correct format (DNAString)
  sequence <- DNAString(sequence)  # Convert to DNAString
  
  # Loop through each codon and find its occurrences
  for (i in 1:4) {
    codon <- codons[i]
    occurrences <- matchPattern(codon, sequence)  # matchPattern for DNAString
    
    # Print out debugging information
    print(paste("Searching for codon:", codon))
    
    if (length(occurrences) > 0) {
      codonpositions <- start(occurrences)  # Directly get start positions
      numoccurrences <- length(codonpositions)
      
      print(paste("Number of occurrences of", codon, ":", numoccurrences))
      
      if (numoccurrences > 0) {
        positions <- append(positions, codonpositions, after = length(positions))
        types <- append(types, rep(codon, numoccurrences), after = length(types))
      }
    }
  }
  
  # Order the positions
  indices <- order(positions)
  positions <- positions[indices]
  types <- types[indices]
  
  mylist <- list(positions, types)
  return(mylist)
}


# Example usage
monkey <- read.fasta(file = "Monkeypox.fasta")
monkeyseq <- monkey[[1]]  # Putting sequence in vector
monkeystartseq <- monkeyseq[1:500]  # Extract the first 500 bases (for testing)

# Convert the sequence to string format
monkeyseqstring <- c2s(monkeystartseq)

# Run the function to find potential start and stop codons
findPotentialstartandstop(monkeyseqstring)




############################################################################################
