
# Moving average ----------------------------------------------------------
ma <- function(x, n=3, sides=2){
  if(length(x) >= n){stats::filter(x, rep(1/n, n), sides=sides)
  }else NA_real_
}


# Z normalisation ---------------------------------------------------------

znorm <- function(ts){
  ts.mean <- mean(ts)
  ts.dev <- sd(ts)
  (ts - ts.mean)/ts.dev
}

# Fix negative edge lengths -----------------------------------------------

fix_negative_edge_length <- function(tree, collapse_multi = FALSE) {
  edge_infos <- cbind(tree$edge, tree$edge.length) %>% data.table::data.table()
  colnames(edge_infos) <- c('from', 'to', 'length')
  negative_edges <- edge_infos[length < 0, sort(unique(from))]
  negative_edges
  
  for (negative_edges in negative_edges) {
    minus_length <- edge_infos[from == negative_edges, ][order(length)][1, length]
    edge_infos[from == negative_edges, length := length - minus_length]
    edge_infos[to == negative_edges, length := length + minus_length]
  }
  tree$edge.length <- edge_infos$length
  if(collapse_multi){
    tree <- di2multi(tree)
  }
  return(tree)
}

# Sliding window from alignment -------------------------------------------
alignment_sw <- function(x, width, interval = 1, maxgaps=Inf){
  alignment_width <- max(width(x))
  win <- seq(1,  alignment_width - width, by = interval) #Get all possible windows
  out <- vector("list", length=length(win))
  for(i in 1:length(win)){
    amplicon <- Biostrings::subseq(x, start=win[i], end = win[i]+width) #may need to add 1 to start
    rem <- names(amplicon)[Biostrings::letterFrequency(amplicon, "-") > maxgaps]
    amplicon <- amplicon[!names(amplicon) %in% rem]
    # remove any windows with all gaps 
    if(any(letterFrequency(amplicon, "ACGT") > 0)){
      out[[i]] <- amplicon
      attr(out[[i]], "window") <- win[i]
    }
  }
  rem <- sapply(out, is.null)
  out <- out[!rem]
  names(out) <- paste(win[!rem], win[!rem]+width, sep="-")
  return(out)
}

# Replace terminal N ------------------------------------------------------
replace_terminal_N <- function(x, new="-"){
  rep_N <- function(seq){
    splitseq <- strsplit(seq, "")[[1]]
    #Find the first real base match
    begin <- match(c("A","C","G","T"), splitseq) -1
    if(length(begin[!is.na(begin)]) == 0) {return(seq)}
    begin <- min(begin[!is.na(begin)])
    # Find the last real base match
    last <- match(c("A","C","G","T"), rev(splitseq)) - 2
    last <- min(last[!is.na(last)])  # Getting the min, ignores deletions in the middle
    if (begin > 0 ){ splitseq[1:begin] <- "-" } 
    if (last > 0){splitseq[(length(splitseq) - last ): length(splitseq)] <- new}
    return(paste(splitseq, collapse=""))
  }
  
  if (is.list(x) & is.raw(x[[1]])){
    res <- insect::char2dna(sapply(insect::dna2char(x),rep_N))
  } else if (is(x, "DNAStringSet")){
    res <- DNAStringSet(sapply(as.character(x), rep_N))
  } else if (length(x) == 1 & is(x, "character")){
    res <- rep_N(x)
  }
  return(res)
}


# Evaluate primer ---------------------------------------------------------


evaluate_primer <- function(x, primer, positions=NULL, direction="F", mm_position=NULL, mm_type=NULL, adjacent=2, count_gaps = FALSE, count_N = FALSE){
  
  # Verify inputs
  if(is.null(mm_position)){
    message("Using the default mismatch position scoring table from PrimerMiner (Position_v1)")
    mm_position <- c(242.4, 202.8, 169.8, 142.4, 119.5, 100.4, 84.5, 71.2,
                     60.2, 51, 43.3, 36.9, 31.6, 27.2, 23.5, 20.4, 17.8, 15.7,
                     13.9, 12.4, 11.2, 10.2, 9.3, 8.6, 8, 7.5, 7.1, 6.7, 6.4,
                     6.2, NA)
    names(mm_position) <- c(seq(1:(length(mm_position)-1)), NA)
    mm_position <- as.data.frame(tibble::enframe(mm_position, name = "pos", value="score"))
  }
  if(is.null(mm_type)){
    message("Using the default mismatch type scoring table from PrimerMiner (Position_v1)")
    mm_type <- as.data.frame(tibble::tribble(
      ~X, ~A, ~T, ~C, ~G,
      "A", 2, NA, 0.5, 2,
      "T", NA, 1, 1, 0.5,
      "C", 0.5, 1, 2, NA,
      "G", 2, 0.5, NA, 2,
    ))
  }
  
  # Define nucleotide table
  upac <- as.data.frame(tibble::tribble(
    ~ID, ~comment, ~A, ~T, ~C, ~G, ~comp,
    "A", "Adenine", 1, 0, 0, 0, "T",
    "C", "Cytosine", 0, 0, 1, 0, "G",
    "G", "Guanine", 0, 0, 0, 1, "C",
    "T", "Thymine", 0, 1, 0, 0, "A",
    "R", "A or G", 0.5, 0, 0, 0.5, "Y",
    "Y", "C or T", 0, 0.5, 0.5, 0, "R",
    "S", "G or C", 0, 0, 0.5, 0.5, "S",
    "W", "A or T", 0.5, 0.5, 0, 0, "W",
    "K", "G or T", 0, 0.5, 0, 0.5, "M",
    "M", "A or C", 0.5, 0, 0.5, 0, "K",
    "B", "C or G or T", 0, 1/3, 1/3, 1/3, "V",
    "D", "A or G or T", 1/3, 1/3, 0, 1/3, "H",
    "H", "A or C or T", 1/3, 1/3, 1/3, 0, "D",
    "V", "A or C or G", 1/3, 0, 1/3, 1/3, "B",
    "N", "any base", 0.25, 0.25, 0.25, 0.25, "N",
    "I", "inosine", 0.25, 0.25, 0.25, 0.25, "N",
  ))
  
  # if positons is character of lenght 2, its start and end position of primer
  if (is(positions, "numeric") & length(positions) == 2){
    start <- sort(positions)[1]
    stop <- sort(positions)[2]
  } else if (is.null(positions)){
    model <- aphid::derivePHMM(x)
    positions <- get_binding_position()
    start <-  sort(positions)[1]
    stop <- sort(positions)[2]
  } else {
    stop("Positions must be a vector of length 2 containing the start and stop position within the alignment, or NULL to automatically detect")
  }
  # check primer length
  if(!(stop+1 - start) == nchar(primer)){
    stop("length of the given region does NOT match the primer length")
  }
  
  # spit primer into nuceleotides
  primer_split <- strsplit(primer, "")[[1]]
  
  # Turn alignment into matrix
  alignment <- matrix(unlist(ape::as.character.DNAbin(x)), ncol = length(x[[1]]), byrow = TRUE)
  rownames(alignment) <- names(x)
  alignment <- apply(alignment,2,toupper)
  
  # extract primer binding region
  primer_region <- alignment[,start:stop]
  
  # If direction is R - make rev comp of alignment
  if(direction=="R"){
    primer_region <- primer_region[,ncol(primer_region):1]
    primer_region <- matrix(upac$comp[match(primer_region, upac$ID)], nrow=nrow(primer_region), ncol=ncol(primer_region), byrow=F)
  }
  
  # Process degeneracy
  primer_degen <- upac[match(primer_split, upac$ID), 3:6] # get wobble scores for primer
  row.names(primer_degen) <- 1:nrow(primer_degen)
  primer_degen <- data.frame(primer_degen > 0) # TRUE = base present in wobble
  
  # duplicate row if only one sequence! - Why is this necesary?
  if (nrow(alignment)==1){
    primer_region <- rbind(primer_region, primer_region)
    primer_region_1 <- TRUE
  } else {
    primer_region_1 <- FALSE
  }
  
  mm_scores <- matrix(ncol = length(primer_split), nrow = length(x))
  
  # Score mismatches for each base in primer
  for(i in 1:length(primer_split)){
    # Check for degenerate bases in sequ
    sequ <- upac[match(primer_region[,i], upac$ID), 3:6]
    sequ <- data.frame(sequ > 0)
    degen <- rowSums(sequ)
    
    sequ[,1] <- sequ[,1]-unlist(primer_degen[i,])[1]*2
    sequ[,2] <- sequ[,2]-unlist(primer_degen[i,])[2]*2
    sequ[,3] <- sequ[,3]-unlist(primer_degen[i,])[3]*2
    sequ[,4] <- sequ[,4]-unlist(primer_degen[i,])[4]*2
    
    # get sequences that do not match the primer
    error <- which(!rowSums(sequ==-1)>0)
    numberofmatches <- rowSums(sequ==-1)
    match_bases <- rep(0, nrow(sequ)) #
    
    # remove bases which arent a match (-2 = not a match)
    sequ[sequ==-2] <- 0
    sequ[sequ==-1] <- 1
    
    wob_sequ_adj_factor <- (rowSums(sequ)-numberofmatches)/rowSums(sequ)
    wob_sequ_adj_factor[error] <- 1
    wob_sequ_adj_factor <- as.vector(wob_sequ_adj_factor) # adjust for wobbles in sequence
    
    match_bases[error] <- mm_position[length(primer_split)+1-i,2]
    match_bases[wob_sequ_adj_factor < 1] <- mm_position[length(primer_split)+1-i,2]
    
    # adjustment of mismatch mm_type
    sequ[,unlist(primer_degen[i,])] <- 0 # keep missmatches
    
    # switch out bases, build rev comp! primer binds on complementary bases : )
    sequ <- sequ[c(2,1,4,3)]
    names(sequ) <- names(sequ)
    mm <- unlist(primer_degen[i,]) # bases present in primer
    
    # Score mismatches
    type_scores <- rep(NA, nrow(sequ)) # empty table
    
    for (m in c(1:4)[mm]){ # cycle trough bases present in primer
      type_temp <- 0 # calculate mm_scores for individual mm of primer base in sequ
      for(n in 2:5){ # acount for wobbles in sequence
        type_temp <- cbind(type_temp, sequ[n-1]*mm_type[m,n])
      }
      type_temp[type_temp==0] <- NA
      n_degen <- rowSums(!is.na(type_temp), na.rm=T) # count number of nucleotides (wobles) in sequ
      n_degen <- 1/n_degen
      n_degen[is.infinite(n_degen)] <- NA
      type_temp2 <- type_temp*n_degen
      type_temp2 <- rowSums(type_temp2, na.rm=T)
      type_scores <- cbind(type_scores, type_temp2)
    }
    
    type_scores <- cbind(type_scores, NA)
    f <- rowSums(!is.na(type_scores), na.rm=T)
    f <- 1/f
    f[is.infinite(f)] <- NA
    type_scores2 <- type_scores*f
    type_scores2 <- rowSums(type_scores2, na.rm=T)
    type_scores2[type_scores2==0] <- 1
    
    match_bases <- match_bases * type_scores2
    
    # Mark gaps as NA?
    if(!count_gaps){
      match_bases[primer_region[,i]=="-"] <- NA
    } else {
      match_bases[primer_region[,i]=="-"] <- mm_position[length(primer_split)+1-i,2]
    }
    
    # Mark N's as NA?
    if(!count_N){
      match_bases[primer_region[,i]=="N"] <- NA
    } else {
      match_bases[primer_region[,i]=="N"] <- mm_position[length(primer_split)+1-i,2]
    }
    
    match_bases <- round(match_bases* wob_sequ_adj_factor, digits=2)
    mm_scores[,i] <- match_bases
  }
  
  # calculate increased error score for adjacent bases, factor has to be bigger than 1!
  if(adjacent > 1){
    for (i in 1:(ncol(mm_scores)-1)){
      adjTT <- paste(mm_scores[,i]>0, mm_scores[,(i+1)]>0) =="TRUE TRUE" # find adjacent values
      mm_scores[,i][adjTT] <- mm_scores[,i][adjTT]* adjacent
      mm_scores[,(i+1)][adjTT] <- mm_scores[,(i+1)][adjTT]* adjacent
    }
    
  }
  # Return primer mismatch mm_scores
  mm_scores <- data.frame(mm_scores)
  names(mm_scores) <- paste("V", length(primer_split):1, sep="")
  
  exp <- 1:nrow(primer_region) # save revcomp of primer region
  
  for(k in 1:nrow(primer_region)){
    exp[k] <- paste(primer_region[k,], collapse="")
  }
  
  mm_scores <- data.frame(template = names(x), "template"= exp, mm_scores, "sum"=rowSums(mm_scores), direction=direction)
  #row.names(mm_scores) <- names(x)
  
  if(primer_region_1){mm_scores <- mm_scores[1,]}# remove duplicated row, as only one sequence
  return(mm_scores)
}
