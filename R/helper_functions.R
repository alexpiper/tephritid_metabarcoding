
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
alignment_sw <- function(x, width, interval = 1, maxgaps=0){
  alignment_width <- max(width(x))
  win <- seq(1,  alignment_width - width, by = interval) #Get all possible windows
  out <- vector("list", length=length(win))
  for(i in 1:length(win)){
    amplicon <- Biostrings::subseq(x, start=win[i], end = win[i]+width) #may need to add 1 to start
    rem <- names(amplicon)[Biostrings::letterFrequency(amplicon, "-") > maxgaps]
    amplicon <- amplicon[!names(amplicon) %in% rem]
    out[[i]] <- amplicon
  }
  names(out) <- paste(win, win+width, sep="-")
  return(out)
}