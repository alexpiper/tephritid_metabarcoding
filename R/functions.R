
# Quality control ---------------------------------------------------------

step_seq_qc <- function(fcid, quiet=FALSE){
  
  seq_dir <- paste0("data/", fcid, "/")
  qc_dir <- paste0("output/logs/", fcid,"/" )
  
  # Check that required files exist
  if(!dir.exists(seq_dir)) {
    stop("seq_dir doesnt exist, check that the correct path was provided")
  }
  if(!dir.exists(paste0(seq_dir,"/InterOp"))){
    stop("InterOp folder must be present to run quality checks")
  }
  
  # Create qc_dir if it doesnt exist
  if(!dir.exists(qc_dir)) {dir.create(qc_dir, recursive = TRUE)}
  
  ## Sequencing run quality check using savR
  fc <- savR::savR(seq_dir)
  
  readr::write_csv(savR::correctedIntensities(fc), normalizePath(paste0(qc_dir, "/correctedIntensities.csv")))
  readr::write_csv(savR::errorMetrics(fc), normalizePath(paste0(qc_dir, "/errorMetrics.csv")))
  readr::write_csv(savR::extractionMetrics(fc), normalizePath(paste0(qc_dir, "/extractionMetrics.csv")))
  readr::write_csv(savR::qualityMetrics(fc), normalizePath(paste0(qc_dir, "/qualityMetrics.csv")))
  readr::write_csv(savR::tileMetrics(fc), normalizePath(paste0(qc_dir, "/tileMetrics.csv")))
  
  gg.avg_intensity <- fc@parsedData[["savCorrectedIntensityFormat"]]@data %>%
    dplyr::group_by(tile, lane) %>%
    dplyr::summarise(Average_intensity = mean(avg_intensity), .groups="drop") %>% 
    dplyr::mutate(side = case_when(
      str_detect(tile, "^11") ~ "Top",
      str_detect(tile, "^21") ~ "Bottom"
    ))%>%
    ggplot(aes(x=lane, y=as.factor(tile), fill=Average_intensity)) +
    geom_tile() +
    facet_wrap(~side, scales="free") +
    scale_fill_viridis_c()
  
  pdf(file=normalizePath(paste0(qc_dir, "/gg.avg_intensity")), width = 11, height = 8 , paper="a4r")
  plot(gg.avg_intensity)
  try(dev.off(), silent=TRUE)
  
  pdf(file=normalizePath(paste0(qc_dir, "/PFclusters.pdf")), width = 11, height = 8 , paper="a4r")
  savR::pfBoxplot(fc)
  try(dev.off(), silent=TRUE)
  
  for (lane in 1:fc@layout@lanecount) {
    pdf(file=normalizePath(paste0(qc_dir, "/QScore_L", lane, ".pdf")), width = 11, height = 8 , paper="a4r")
    savR::qualityHeatmap(fc, lane, 1:fc@directions)
    try(dev.off(), silent=TRUE)
  } 
  if(!quiet){message("Flow cell quality metrics written to: ", qc_dir)}
  
  # Return the total number of reads and passing filter
  out <- fc@parsedData[["savTileFormat"]]@data %>%
    dplyr::filter(code %in% c(100,101)) %>%
    dplyr::mutate(code = case_when(
      code == 100 ~ "reads_total",
      code == 101 ~ "reads_pf"
    ),fcid = stringr::str_remove(fc@flowcell, "^.*-")) %>% 
    group_by(fcid, code) %>%
    dplyr::summarise(reads = sum(value), .groups="drop") %>%
    tidyr::pivot_wider(names_from = code,
                       values_from = reads)
  return(out)
}

step_sample_qc <- function(sample_id, fcid, multithread=FALSE, quiet=FALSE){
  
  seq_dir <- paste0("data/", fcid, "/")
  qc_dir <- paste0("output/logs/", fcid,"/" )
  
  print(qc_dir)
  print(seq_dir)
  
  fq <- paste0(seq_dir, sample_id, "*.fastq.gz")
  # Check that required files exist
  if(!dir.exists(seq_dir)) {
    stop("seq_dir doesnt exist, check that the correct path was provided")
  }
  # Create qc_dir if it doesnt exist
  if(!dir.exists(qc_dir)) {dir.create(qc_dir, recursive = TRUE)}
  
  # Handle multithreading
  cores <- setup_multithread(multithread)
  
  # Define fastqc function
  fastqc <- function (fq, qc.dir = NULL, threads = 1, fastqc.path = "bin/FastQC/fastqc") {
    if (is.null(qc.dir))  qc.dir <- file.path(fq.dir, "FASTQC")

    if (.Platform$OS.type == "unix") {
      cmd <- paste0(fastqc.path, " ", fq, "  --threads ", 
                    threads, " --outdir ", qc.dir)
      result <- system(cmd)
    }
    else {
      .threads <- paste0("--threads ", threads)
      .qc.dir <- paste0("--outdir ", qc.dir)
      result <- processx::run(command = "perl", args = c(fastqc.path, 
                                                         fq, .threads, .qc.dir), echo = TRUE, echo_cmd = TRUE, 
                              spinner = TRUE, windows_verbatim_args = TRUE, error_on_status = FALSE, 
                              cleanup_tree = TRUE)
    }
    return(result)
  }
  
  # Run fastqc
  qc_out <- fastqc(fq = fq, qc.dir= qc_dir, fastqc.path = "bin/FastQC/fastqc", threads=cores)

}

step_multiqc <- function(fcid, quiet=FALSE){
  qc_dir <- paste0("output/logs/", fcid,"/" )
  # Check that required files exist
  if(!dir.exists(qc_dir)) {
    stop("qc_dir doesnt exist, check that the correct path was provided")
  }
  # Write out fastqc report
  if(!quiet){
    ngsReports::writeHtmlReport(qc_dir, overwrite = TRUE, gcType ="Genome",  quiet=quiet)
    message("Sample quality metrics written to: ", qc_dir)
  } else {
    suppressMessages(ngsReports::writeHtmlReport(qc_dir, overwrite = TRUE, gcType ="Genome",  quiet=quiet))
  }
  
}




step_switching_calc <- function(fcid, multithread=FALSE, quiet=FALSE){
  
  seq_dir <- paste0("data/", fcid, "/")
  qc_dir <- paste0("output/logs/", fcid,"/" )
  
  # Check that required files exist
  if(!dir.exists(seq_dir)) {
    stop("seq_dir doesnt exist, check that the correct path was provided")
  }
  if(!any(str_detect(list.files(seq_dir, pattern="_R1_", full.names = TRUE), "Undetermined"))){
    stop("Error, an Undetermined reads fastq must be present to calculate index switching")
  }
  # Create qc_dir if it doesnt exist
  if(!dir.exists(qc_dir)) {dir.create(qc_dir, recursive = TRUE)}
  
  # Handle multithreading
  cores <- setup_multithread(multithread)
  
  # Summarise indices
  indices <- sort(list.files(seq_dir, pattern="_R1_", full.names = TRUE)) %>%
    purrr::set_names() %>%
    furrr::future_map(seqateurs::summarise_index) %>%
    dplyr::bind_rows(.id="Sample_Name")%>%
    dplyr::arrange(desc(Freq)) %>% 
    dplyr::mutate(Sample_Name = Sample_Name %>% 
                    stringr::str_remove(pattern = "^(.*)\\/") %>%
                    stringr::str_remove(pattern = "(?:.(?!_S))+$"))
  
  combos <- indices %>% 
    dplyr::filter(!str_detect(Sample_Name, "Undetermined")) %>%
    dplyr::select(index, index2) %>%
    tidyr::expand(index, index2)
  
  #get unused combinations resulting from index switching
  switched <- combos %>%
    dplyr::left_join(indices, by=c("index", "index2")) %>%
    tidyr::drop_na()
  
  other_reads <- anti_join(indices,combos, by=c("index", "index2")) %>%
    dplyr::summarise(sum = sum(Freq)) %>%
    dplyr::pull(sum)
  
  #Summary of index switching rate
  
  res <- switched %>%
    dplyr::mutate(type = case_when(
      !stringr::str_detect(Sample_Name, "Undetermined") ~ "expected",
      stringr::str_detect(Sample_Name, "Undetermined") ~ "observed"
    )) %>%
    dplyr::group_by(type) %>%
    dplyr::summarise(switch_rate = sum(Freq), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = type,
                       values_from = switch_rate) %>%
    dplyr::mutate(switch_rate =  observed / expected )
  
  readr::write_csv(res, normalizePath(paste0(qc_dir, "index_switching.csv")))
  
  if(!quiet){message("Index switching rate calculated as: ", res$switch_rate)}
  
  #Plot switching
  # Need to find a way to collapse all indexes by hamming distance for this plot
  gg.switch <- switched %>%
    ggplot(aes(x = index, y = index2), stat="identity") +
    geom_tile(aes(fill = Freq),alpha=0.8)  + 
    scale_fill_viridis_c(name="log10 Reads", begin=0.1, trans="log10")+
    theme(axis.text.x = element_text(angle=90, hjust=1), 
          plot.title=element_text(hjust = 0.5),
          plot.subtitle =element_text(hjust = 0.5)#,
          #legend.position = "none"
    ) +
    labs(title= fcid, subtitle = paste0(
      "Total Reads: ", sum(indices$Freq),
      ", Switch rate: ", sprintf("%1.4f%%", res$switch_rate*100),
      ", other Reads: ", other_reads)) 
  pdf(file=normalizePath(paste0(qc_dir, "/index_switching.pdf")), width = 11, height = 8 , paper="a4r")
  plot(gg.switch)
  try(dev.off(), silent=TRUE)
  return(res)
}


# Plots -------------------------------------------------------------------

plot_read_quals <- function(sample_id, input_dir, truncLen = NULL, quiet=FALSE){
  input_dir <- normalizePath(input_dir)
  # Seq dir might need ot be changed to trimmed or other sub folders
  fastqFs <- fs::dir_ls(path = input_dir, glob = paste0("*",sample_id, "*_R1_*"))
  fastqRs <- fs::dir_ls(path = input_dir, glob = paste0("*",sample_id, "*_R2_*"))
  if(!file.size(fastqFs) > 28) {
    message(paste0("Sample ", sample_id, " Has no reads"))
    return(NULL)
  }

  #Plot qualities
  gg.Fqual <- plot_quality(fastqFs) +
    ggtitle(paste0(sample_id, " Forward Reads")) +
    scale_x_continuous(breaks=seq(0,300,25))

  gg.Fee <- plot_maxEE(fastqFs) + 
    ggtitle(paste0(sample_id, " Forward Reads")) +
    scale_x_continuous(breaks=seq(0,300,25)) +
    theme(legend.position = "bottom")
  
  gg.Rqual <- plot_quality(fastqRs) + 
    ggtitle(paste0(sample_id, " Reverse Reads")) +
    scale_x_continuous(breaks=seq(0,300,25)) +
    theme(legend.position = "bottom")
  
  gg.Ree <- plot_maxEE(fastqRs) +
    ggtitle(paste0(sample_id, " Reverse Reads")) +
    scale_x_continuous(breaks=seq(0,300,25)) +
    theme(legend.position = "bottom")
  
  if(!is.null(truncLen)){
    gg.Fqual <- gg.Fqual +
      geom_vline(aes(xintercept=truncLen[1]), colour="blue") +
      annotate("text", x = truncLen[1]-10, y =2, label = paste0("Suggested truncLen = ", truncLen[1]), colour="blue")
    gg.Fee <-gg.Fee +
      geom_vline(aes(xintercept=truncLen[1]), colour="blue")+
      annotate("text", x = truncLen[1]-10, y =-3, label = paste0("Suggested truncLen = ", truncLen[1]), colour="blue")
    gg.Rqual <- gg.Rqual +
      geom_vline(aes(xintercept=truncLen[2]), colour="blue")+
      annotate("text", x = truncLen[1]-10, y =2, label = paste0("Suggested truncLen = ", truncLen[2]), colour="blue")
    gg.Ree <- gg.Ree + 
      geom_vline(aes(xintercept=truncLen[2]), colour="blue")+
      annotate("text", x = truncLen[1]-10, y =-3, label = paste0("Suggested truncLen = ", truncLen[2]), colour="blue")
  }
  
  Qualplots <- (gg.Fqual + gg.Rqual) / (gg.Fee + gg.Ree)
  return(Qualplots)
}


# Primer trimming ---------------------------------------------------------

step_primer_trim <- function(sample_id, seq_dir, qc_dir, for_primer_seq, rev_primer_seq, quiet=FALSE){
  # Check inputs
  if(!is.character(for_primer_seq) | !length(for_primer_seq)==1){
    stop("for_primer_seq must be a character vector of length 2, containing forward and reverse primers")
  }
  if(!is.character(rev_primer_seq) | !length(rev_primer_seq)==1){
    stop("rev_primer_seq must be a character vector of length 2, containing forward and reverse primers")
  }
  
  # Check that required files exist
  if(!dir.exists(seq_dir)) {
    stop("seq_dir doesnt exist, check that the correct path was provided")
  }
  # Create qc_dir if it doesnt exist
  if(!dir.exists(qc_dir)) {dir.create(qc_dir, recursive = TRUE)}
  
  primers <- c(for_primer_seq, rev_primer_seq)
  primers <- str_replace_all(primers, "I", "N")
  
  fastqFs <- fs::dir_ls(path = seq_dir, glob = paste0("*",sample_id, "*_R1_*"))
  fastqRs <- fs::dir_ls(path = seq_dir, glob = paste0("*",sample_id, "*_R2_*"))
  if(length(fastqFs) != length(fastqRs)) stop(paste0("Forward and reverse files for ",sample_id," do not match."))
  
  #Check if there were more than 1 primer per sample
  multi_primer <- any(stringr::str_detect(primers, ";"))
  
  if (multi_primer) {
    # Create output directory
    demuxpath <- file.path(seq_dir, "00_demux")
    dir.create(demuxpath)
    
    Fprimers <- unlist(str_split(unique(primers[1]), ";"))
    Rprimers <- unlist(str_split(unique(primers[2]), ";"))
    primer_names <- unlist(str_split(unique(primers), ";"))
    
    # Demultiplex reads
    demux <- bbdemux2(install="bin/bbmap", fwd=fastqFs, rev=fastqRs, Fbarcodes = Fprimers, 
                      Rbarcodes = Rprimers, names=primer_names, degenerate=TRUE, out.dir=demuxpath,
                      threads=1, mem=4,  hdist=0, force=TRUE)
    
    # Rename output files to match the required pipeline format
    old_names <- fs::dir_ls(path = demuxpath, glob = paste0("*",sample_id, "*_R1R2_*"))
    new_names <- old_names %>%
      basename() %>%
      stringr::str_remove(".fastq.gz") %>%
      stringr::str_split("_", n=Inf) %>%
      purrr::map_chr(function(x){
        paste0(paste0(c(x[1:3], x[length(x)], x[4:(length(x)-1)]), collapse = "_"),".fastq.gz")
      }) 
    if(file.exists(file.path(dirname(old_names), new_names))){
      file.remove(file.path(dirname(old_names), new_names))
    }
    file.rename(old_names, file.path(dirname(old_names), new_names))
    
    # Trim primers from demultiplexed fastq
    demux_fastqs <- fs::dir_ls(path = demuxpath, glob = paste0("*",sample_id, "*_R1R2_*"))
    res <- bbtrim2(install="bin/bbmap", fwd = demux_fastqs,
                   primers = c(Fprimers, Rprimers), checkpairs = FALSE,
                   degenerate = TRUE, out.dir=file.path(seq_dir, "01_trimmed"), trim.end = "left", 
                   kmer=NULL, tpe=TRUE, tbo=TRUE,
                   ordered = TRUE, mink = FALSE, hdist = 2,
                   maxlength =(max(step_data$for_read_length, step_data$rev_read_length
                   ) - sort(nchar(c(Fprimers, Rprimers)), decreasing = FALSE)[1]) +5, force = TRUE, quiet=FALSE)
    
    # Re-split interleaved fastq's
    trimmedpath <- file.path(seq_dir, "01_trimmed") 
    trimmed_fastqs<- fs::dir_ls(path = trimmedpath, glob = paste0("*",sample_id, "*_R1R2_*"))
    bbsplit2(install="bin/bbmap", file=trimmed_fastqs, force=TRUE)
    
  } else if (!multi_primer) {
    Fprimers <- primers[1]
    Rprimers <- primers[2]
    
    res <- bbtrim2(install="bin/bbmap", fwd = fastqFs, rev = fastqRs,
                   primers = c(Fprimers, Rprimers), checkpairs = FALSE,
                   degenerate = TRUE, out.dir=file.path(seq_dir, "01_trimmed"), trim.end = "left", 
                   kmer=NULL, tpe=TRUE, tbo=TRUE,
                   ordered = TRUE, mink = FALSE, hdist = 2,
                   maxlength =(max(step_data$for_read_length, step_data$rev_read_length
                   ) - sort(nchar(c(Fprimers, Rprimers)), decreasing = FALSE)[1]) +5,
                   force = TRUE, quiet=FALSE)
  }
  return(res)
}




step_primer_trim2 <- function(sample_id, input_dir, output_dir, qc_dir, for_primer_seq, rev_primer_seq, pcr_primers,
                             n = 1e6, qualityType = "Auto", check_paired = TRUE, compress =TRUE, quiet=FALSE){
  input_dir <- normalizePath(input_dir)
  output_dir <- normalizePath(output_dir)
  qc_dir <- normalizePath(qc_dir)
  # Check inputs
  if(!is.character(for_primer_seq) | !length(for_primer_seq)==1){
    stop("for_primer_seq must be a character vector of length 2, containing forward and reverse primers")
  }
  if(!is.character(rev_primer_seq) | !length(rev_primer_seq)==1){
    stop("rev_primer_seq must be a character vector of length 2, containing forward and reverse primers")
  }
  # Check that required files exist
  if(!dir.exists(input_dir)) {
    stop("input_dir doesnt exist, check that the correct path was provided")
  }
  
  # Create output directory if it doesnt exist
  if(!dir.exists(output_dir)) {dir.create(output_dir)}

  # Create qc_dir if it doesnt exist
  if(!dir.exists(qc_dir)) {dir.create(qc_dir, recursive = TRUE)}
  
  fastqFs <- normalizePath(fs::dir_ls(path = input_dir, glob = paste0("*",sample_id, "*_R1_*")))
  fastqRs <- normalizePath(fs::dir_ls(path = input_dir, glob = paste0("*",sample_id, "*_R2_*")))
  if(length(fastqFs) != length(fastqRs)) stop(paste0("Forward and reverse files for ",sample_id," do not match."))
  
  #Check if there were more than 1 primer per sample
  multi_primer <- any(stringr::str_detect(for_primer_seq, ";"), stringr::str_detect(rev_primer_seq, ";"))
  
  if (multi_primer) {
    for_primer_seq <- unlist(stringr::str_split(for_primer_seq, ";"))
    rev_primer_seq <- unlist(stringr::str_split(rev_primer_seq, ";"))
    
    primer_names <- unlist(stringr::str_split(pcr_primers, ";"))
    
    # Define outfiles
    fwd_out <- purrr::map(primer_names, ~{
      normalizePath(paste0(output_dir,"/", str_replace(basename(fastqFs), ".fastq", paste0("_",.x, ".fastq"))))
    }) %>%
      unlist()
    rev_out <- purrr::map(primer_names, ~{
      normalizePath(paste0(output_dir,"/", str_replace(basename(fastqRs), ".fastq", paste0("_",.x, ".fastq"))))
    }) %>%
      unlist()
    
  } else if (!multi_primer) {
    fwd_out <- normalizePath(paste0(output_dir,"/", basename(fastqFs)))
    rev_out <- normalizePath(paste0(output_dir,"/", basename(fastqRs)))
  } 
  
  # Demultiplex reads and trim primers
  res <- trim_primers(fwd = fastqFs,
                     rev = fastqRs,
                     fwd_out = fwd_out,
                     rev_out = rev_out,
                     for_primer_seq = str_replace_all(for_primer_seq, "I", "N"),
                     rev_primer_seq = str_replace_all(rev_primer_seq, "I", "N"),
                     n = n, qualityType = qualityType, check_paired = check_paired,
                     compress =compress, quiet=quiet
  )
  return(res)
}


step_read_filter <- function(sample_id, input_dir, output_dir, quiet=FALSE, ...){
  input_dir <- normalizePath(input_dir)
  output_dir <- normalizePath(output_dir)
  
  # Check that required files exist
  if(!dir.exists(input_dir)) {
    stop("input_dir doesnt exist, check that the correct path was provided")
  }
  # Create output directory if it doesn't exist
  if(!dir.exists(output_dir)) {dir.create(output_dir)}
  
  
  fastqFs <- fs::dir_ls(path = input_dir, glob = paste0("*",sample_id, "*_R1_*"))
  fastqRs <- fs::dir_ls(path = input_dir, glob = paste0("*",sample_id, "*_R2_*"))
  
  if(length(fastqFs) != length(fastqRs)) stop(paste0("Forward and reverse files for ",sample_id," do not match."))

  # Run read filter
  res <- dada2::filterAndTrim(
    fwd = fastqFs, filt = file.path(output_dir, basename(fastqFs)),
    rev = fastqRs, filt.rev = file.path(output_dir, basename(fastqRs)),
    rm.phix = TRUE, multithread = FALSE, compress = TRUE, verbose = !quiet, ...) %>%
    as_tibble() %>%
    dplyr::rename(filter_input = reads.in,
                  filter_output = reads.out)
  
  filtered_summary <- res %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    dplyr::mutate(reads_remaining = signif(((filter_output / filter_input) * 100), 2)) %>%
    dplyr::filter(!is.na(reads_remaining))
  
  if (filtered_summary$reads_remaining < 10) {
    message(paste0("WARNING: Less than 10% reads remaining for ",
                   filtered_summary$sample), "Check filtering parameters ")
  }
  return(res)
}

trim_primers <- function(fwd, rev, fwd_out, rev_out, for_primer_seq, rev_primer_seq, 
            n = 1e6, qualityType = "Auto", check_paired = TRUE, id.field = NULL, id.sep = "\\s", compress =TRUE, quiet=FALSE){
  
  ## iterating through forward and reverse files using fastq streaming
  fF <- ShortRead::FastqStreamer(file.path(fwd), n = n)
  on.exit(close(fF))
  fR <- ShortRead::FastqStreamer(file.path(rev), n = n)
  on.exit(close(fR), add=TRUE)
  
  #Check if there were more than 1 primer per sample
  if(any(length(fwd_out), length(rev_out), length(for_primer_seq), length(rev_primer_seq) >0)){
    multi_primer <- TRUE
  }
  
  # Check if number of F and R primers match the number of outfiles
  if (!all.equal(length(fwd_out), length(rev_out), length(for_primer_seq), length(rev_primer_seq))) {
      stop("fwd_out and rev_out must be the same length as for_primer_seq and rev_primer_seq")
  }
  
  #if(!C_isACGT(primer)) stop("Non-ACGT characters detected in primers")
  
  # Delete output files if they already exist
  c(fwd_out, rev_out) %>%
    purrr::walk(remove_if_exists, quiet=quiet)
  
  first=TRUE
  append <- vector("logical", length= length(fwd_out)) 
  remainderF <- ShortRead::ShortReadQ(); remainderR <- ShortRead::ShortReadQ()
  casava <- "Undetermined" #ID field format
  # Setup read tracking
  inseqs <- 0
  outseqs <- vector("numeric", length= length(fwd_out))
  
  while( TRUE ) {
    suppressWarnings(fqF <- ShortRead::yield(fF, qualityType = qualityType))
    suppressWarnings(fqR <- ShortRead::yield(fR, qualityType = qualityType))
    if(length(fqF) == 0 && length(fqR) == 0) { break } # Stop loop if theres no reads left to process
    
    inseqs <- inseqs + length(fqF)
    
    # Make sure that forward and reverse reads are correctly paired
    # Determine the sequence identifier field. Looks for a single 6-colon field (CASAVA 1.8+ id format)
    if(check_paired) {
      if(first) { 
        if(is.null(id.field)) {
          # or a single 4-colon field (earlier format). Fails if it doesn't find such a field.
          id1 <- as.character(ShortRead::id(fqF)[[1]])
          id.fields <- strsplit(id1, id.sep)[[1]]
          ncolon <- sapply(gregexpr(":", id.fields), length)
          ncoltab <- table(ncolon)
          if(max(ncolon) == 6 && ncoltab["6"] == 1) { # CASAVA 1.8+ format
            casava <- "Current"
            id.field <- which(ncolon == 6)
          } else if (max(ncolon) == 4 && ncoltab["4"] == 1) { # CASAVA <=1.7 format
            casava <- "Old"
            id.field <- which(ncolon == 4)
          } else { # Couldn't unambiguously find the seq id field
            stop("Couldn't automatically detect the sequence identifier field in the fastq id string.")
          }
        }
      } else { # !first
        # Prepend the unmatched sequences from the end of previous chunks
        # Need ShortRead::append or the method is not dispatched properly
        fqF <- ShortRead::append(remainderF, fqF)
        fqR <- ShortRead::append(remainderR, fqR)
      }
    } else { # !check_paired
      if(length(fqF) != length(fqR)) stop("Mismatched forward and reverse sequence files: ", length(fqF), ", ", length(fqR), ".")
    }
    
    # Enforce id matching (ASSUMES SAME ORDERING IN F/R, BUT ALLOWS DIFFERENTIAL MEMBERSHIP)
    # Keep the tail of unmatched sequences (could match next chunk)
    if(check_paired) {
      idsF <- sapply(strsplit(as.character(ShortRead::id(fqF)), id.sep), `[`, id.field)
      idsR <- sapply(strsplit(as.character(ShortRead::id(fqR)), id.sep), `[`, id.field)
      if(casava == "Old") { # Drop the index number/pair identifier (i.e. 1=F, 2=R)
        idsF <- sapply(strsplit(idsF, "#"), `[`, 1)
      }
      lastF <- max(c(0,which(idsF %in% idsR)))
      lastR <- max(c(0,which(idsR %in% idsF)))
      if(lastF < length(fqF)) {
        remainderF <- fqF[(lastF+1):length(fqF)]
      } else {
        remainderF <- ShortRead::ShortReadQ() 
      }
      if(lastR < length(fqR)) {
        remainderR <- fqR[(lastR+1):length(fqR)]
      } else {
        remainderR <- ShortRead::ShortReadQ() 
      }
      fqF <- fqF[idsF %in% idsR]
      fqR <- fqR[idsR %in% idsF]
    }
    
    # Demultiplex and trim each primer
    # Only keep reads where primer is detected
    
    for (p in 1:length(for_primer_seq)){
      barlenF <- nchar(for_primer_seq[p])
      barlenR <- nchar(rev_primer_seq[p])
      
      keepF <- Biostrings::neditStartingAt(
        pattern=Biostrings::DNAString(for_primer_seq[p]),
        subject= IRanges::narrow(sread(fqF), 1, barlenF),
        starting.at=1,
        with.indels=FALSE,
        fixed=FALSE )==0
      
      keepR <- Biostrings::neditStartingAt(
        pattern= Biostrings::DNAString(rev_primer_seq[p]),
        subject= IRanges::narrow(sread(fqR), 1, barlenR),
        starting.at=1,
        with.indels=FALSE,
        fixed=FALSE )==0
      
      # Only keep reads where forward primer is present in F, and reverse in R
      keep <- keepF & keepF
      
      fqF_primer <- suppressWarnings(ShortRead::ShortReadQ(sread=ShortRead::sread(fqF[keep]), 
                       quality=Biostrings::quality(Biostrings::quality(fqF[keep])),
                       id=ShortRead::id(fqF[keep])))
      
      fqR_primer <- suppressWarnings(ShortRead::ShortReadQ(sread=ShortRead::sread(fqR[keep]), 
                        quality=Biostrings::quality(Biostrings::quality(fqR[keep])),
                        id=ShortRead::id(fqR[keep])))
  
      # Trim primers from left side
      startF <- max(1, barlenF + 1, na.rm=TRUE)
      startR <- max(1, barlenR + 1, na.rm=TRUE)
      
      # Make sure reads are longer than the primer sequence
      keep <- (width(fqF_primer) >= startF & width(fqR_primer) >= startR)
      fqF_primer <- fqF_primer[keep]
      fqF_primer <- narrow(fqF_primer, start = startF, end = NA)
      fqR_primer <- fqR_primer[keep]
      fqR_primer <- narrow(fqR_primer, start = startR, end = NA)
      
      outseqs[p] <- outseqs[p] + length(fqF_primer)    
      
      if(!append[p]) {
        ShortRead::writeFastq(fqF_primer, fwd_out[p], "w", compress = compress)
        ShortRead::writeFastq(fqR_primer, rev_out[p], "w", compress = compress)
        append[p]=TRUE
        first=FALSE
      } else {
        ShortRead::writeFastq(fqF_primer, fwd_out[p], "a", compress = compress)
        ShortRead::writeFastq(fqR_primer, rev_out[p], "a", compress = compress)
      }
    }
  }
  
  if(!quiet) {
    outperc <- purrr::map(outseqs, ~{
      round(.x * 100 / inseqs, 1)
    }) %>%
      unlist()
    outperc <- paste(" (", outperc, "%),", sep="")
    message("Read in ", inseqs, " paired-sequences, output ", paste(" ", outseqs, ",", sep=""), " ", outperc, " primer-trimmed paired-sequences.", sep="")
  }
  
  if(sum(outseqs)==0) {
    message(paste("No reads remaining for:", fwd, "and", rev))
    file.remove(fwd_out)
    file.remove(rev_out)
  }
  out <- data.frame(for_primer_seq = for_primer_seq,
                    rev_primer_seq = rev_primer_seq,
                    trimmed_input = inseqs,
                    trimmed_output = outseqs)
  return(out)
}


#  DADA2 ------------------------------------------------------------------


step_dada2 <- function(fcid, input_dir, output, qc_dir, nbases=1e+08, randomize=FALSE, multithread=FALSE, pool="pseudo",   quiet=FALSE){
  
  input_dir <- normalizePath(input_dir)
  output <- normalizePath(output)
  qc_dir <- normalizePath(qc_dir)
  filtFs <- list.files(input_dir, pattern="R1_001.*", full.names = TRUE)
  filtRs <- list.files(input_dir, pattern="R2_001.*", full.names = TRUE)
  
  # Learn error rates from a subset of the samples and reads (rather than running self-consist with full dataset)
  errF <- learnErrors(filtFs, multithread = multithread, nbases = nbases,
                      randomize = randomize, qualityType = "FastqQuality", verbose=TRUE)
  errR <- learnErrors(filtRs, multithread = multithread, nbases = nbases,
                      randomize = randomize, qualityType = "FastqQuality", verbose=TRUE)
  
  #write out errors for diagnostics
  write_csv(as.data.frame(errF$trans), paste0(qc_dir, "/", fcid, "_errF_observed_transitions.csv"))
  write_csv(as.data.frame(errF$err_out), paste0(qc_dir, "/", fcid, "_errF_inferred_errors.csv"))
  write_csv(as.data.frame(errR$trans), paste0(qc_dir, "/", fcid, "_errR_observed_transitions.csv"))
  write_csv(as.data.frame(errR$err_out), paste0(qc_dir, "/", fcid, "_errR_inferred_errors.csv"))
  
  ##output error plots to see how well the algorithm modelled the errors in the different runs
  p1 <- plotErrors(errF, nominalQ = TRUE) + ggtitle(paste0(fcid, " Forward Reads"))
  p2 <- plotErrors(errR, nominalQ = TRUE) + ggtitle(paste0(fcid, " Reverse Reads"))
  pdf(paste0(qc_dir,"/",fcid,"_errormodel.pdf"), width = 11, height = 8 , paper="a4r")
  plot(p1)
  plot(p2)
  try(dev.off(), silent=TRUE)
  
  #Error inference and merger of reads
  dadaFs <- dada(filtFs, err = errF, multithread = multithread, pool = pool, verbose = TRUE)
  dadaRs <- dada(filtRs, err = errR, multithread = multithread, pool = pool, verbose = TRUE)
  
  # merge reads
  mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE, minOverlap = 12, trimOverhang = TRUE) 
  mergers <- mergers[sapply(mergers, nrow) > 0]
  bind_rows(mergers, .id="Sample") %>%
    mutate(Sample = str_replace(Sample, pattern="_S.*$", replacement="")) %>%
    write_csv(paste0(qc_dir, "/", fcid, "_mergers.csv"))
  
  #Construct sequence table
  seqtab <- makeSequenceTable(mergers)
  saveRDS(seqtab, output)
  
  # Track reads
  getN <- function(x) sum(getUniques(x))
  res <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN)) %>%
    magrittr::set_colnames(c("dadaFs", "dadaRs", "merged")) %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    mutate(sample_id = str_replace(basename(sample_id), pattern="_S.*$", replacement=""))
  return(res)
}


# ASV filtering -----------------------------------------------------------

# Group by target loci and apply
step_filter_asvs <- function(seqtab, output, qc_dir, min_length = NULL, max_length = NULL, 
                            check_frame=FALSE, genetic_code="SGC4", phmm=NULL, primers=NULL, multithread=FALSE, quiet=FALSE){
  # Print arguments for testing
  #print(as.list(match.call()))
  
  #argg <- c(as.list(environment()), list(...))
  # print(argg)
  
  if(is.matrix(seqtab) | is.data.frame(seqtab)){
    if(!quiet){message("Input is a matrix or data frame")}
  } else if (is.character(seqtab) & str_detect(seqtab, ".rds")){
    seqtab <- readRDS(seqtab)
  } else {
    stop("seqtab must be a matrix/data frame or .rds file")
  }
  reads_starting <- rowSums(seqtab)
  
  # Remove chimeras  
  seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=multithread, verbose=!quiet)
  if(!quiet){message(paste(sum(seqtab_nochim)/sum(seqtab),"of starting abundance remaining after chimera removal"))}
  reads_chimerafilt <- rowSums(seqtab_nochim)
  
  #cut to expected size allowing for some codon indels - do separately for each loci
  if(any(!is.null(c(min_length, max_length)))){
    if(!is.null(min_length) & !is.null(max_length)){
      seqtab_cut <- seqtab_nochim[,nchar(colnames(seqtab_nochim)) %in% min_length:max_length]
    } else if(is.null(min_length) & !is.null(max_length)){
      seqtab_cut <- seqtab_nochim[,nchar(colnames(seqtab_nochim)) < max_length]
    } else if(!is.null(min_length) & is.null(max_length)){
      seqtab_cut <- seqtab_nochim[,nchar(colnames(seqtab_nochim)) > min_length]
    }
    if(!quiet){
      seqs_rem <- length(colnames(seqtab_cut))/length(colnames(seqtab))
      abund_rem <- sum(seqtab_cut)/sum(seqtab)
      message(paste(seqs_rem, " of initial sequences and ",abund_rem, " of initial abundance remaining after length filtering"))
    }
    reads_lengthfilt <- rowSums(seqtab_cut)
  } else {
    seqtab_cut <- seqtab_nochim
    reads_lengthfilt <- rep(NA_integer_,nrow(seqtab)) 
  }
  
  # Load in profile hidden markov model if provided
  if(is.character(phmm) && str_detect(phmm, ".rds")){
    model <- readRDS(phmm)
  } else if (is(phmm, "PHMM")){
    model <- phmm
  }
  
  # subset PHMM if primers were provided
  if (is(model, "PHMM") && !is.null(primers)){
    model <- taxreturn::subset_model(model, primers = primers)
  }
  
  # Align against phmm
  if (is(model, "PHMM")){
    seqs <- DNAStringSet(colnames(seqtab_cut))
    names(seqs) <- colnames(seqtab_cut)
    phmm_filt <- taxreturn::map_to_model(
      seqs, model = model, min_score = 100, min_length = 100,
      shave = FALSE, check_frame = check_frame, kmer_threshold = 0.5, k=5, extra = "fill")
    seqtab_phmm <- seqtab_cut[,colnames(seqtab_cut) %in% names(phmm_filt)]
    if(!quiet){
      seqs_rem <- length(colnames(seqtab_phmm))/length(colnames(seqtab))
      abund_rem <- sum(seqtab_phmm)/sum(seqtab)
      message(paste(seqs_rem, " of initial sequences and ",abund_rem, " of initial abundance remaining after PHMM filtering"))
    }
    reads_phmmfilt <- rowSums(seqtab_phmm)
  } else {
    seqtab_phmm <- seqtab_cut
    reads_phmmfilt <- rep(NA_integer_,nrow(seqtab)) 
  }
  
  #Filter sequences containing stop codons
  if(check_frame){
    seqs <- DNAStringSet(colnames(seqtab_phmm))
    names(seqs) <- colnames(seqtab_phmm)
    codon_filt <- taxreturn::codon_filter(seqs, genetic_code = genetic_code) 
    seqtab_final <- seqtab_phmm[,colnames(seqtab_phmm) %in% names(codon_filt)]
    if(!quiet){
      seqs_rem <- length(colnames(seqtab_final))/length(colnames(seqtab))
      abund_rem <- sum(seqtab_final)/sum(seqtab)
      message(paste(seqs_rem, " of initial sequences and ",abund_rem, " of initial abundance remaining after checking reading frame"))
    }
    reads_framefilt <- rowSums(seqtab_final)
  } else {
    seqtab_final <- seqtab_phmm
    reads_framefilt <- rep(NA_integer_,nrow(seqtab)) 
  }
  reads_final <- rowSums(seqtab_final)
  
  saveRDS(seqtab_final, output)
  
  # Output a cleanup summary
  cleanup <- seqtab %>%
    as.data.frame() %>%
    pivot_longer( everything(),
                  names_to = "OTU",
                  values_to = "Abundance") %>%
    group_by(OTU) %>%
    summarise(Abundance = sum(Abundance)) %>%
    mutate(length  = nchar(OTU)) %>%
    mutate(type = case_when(
      !OTU %in% getSequences(seqtab_nochim) ~ "Chimera",
      !OTU %in% getSequences(seqtab_cut) ~ "Incorrect size",
      !OTU %in% getSequences(seqtab_phmm) ~ "PHMM",
      !OTU %in% getSequences(seqtab_final) ~ "Stop codons",
      TRUE ~ "Real"
    )) 
  write_csv(cleanup, paste0(qc_dir,"/ASV_cleanup_summary.csv"))
  
  # Output length distribution plots
  gg.abundance <- ggplot(cleanup, aes(x=length, y=Abundance, fill=type))+
    geom_bar(stat="identity") + 
    labs(title = "Abundance of sequences")
  
  gg.unique <- ggplot(cleanup, aes(x=length, fill=type))+
    geom_histogram() + 
    labs(title = "Number of unique sequences")
  
  pdf(paste0(qc_dir,"/seqtab_length_dist.pdf"), width = 11, height = 8 , paper="a4r")
  plot(gg.abundance / gg.unique)
  try(dev.off(), silent=TRUE)
  
  # Create output
  res <- tibble(
    sample_id = rownames(seqtab) %>% str_remove(pattern="_S.*$"),
    reads_starting = reads_starting,
    reads_chimerafilt = reads_chimerafilt,
    reads_lengthfilt = reads_lengthfilt,
    reads_phmmfilt = reads_phmmfilt,
    reads_framefilt = reads_framefilt,
    reads_final = reads_final
  )
  return(list(filtered_seqtab = seqtab_final, 
              filtered_asvs = res))
}


# Taxonomic assignment ----------------------------------------------------

# Group by target loci and apply
step_assign_taxonomy <- function(seqtab, output=NULL, qc_dir, database, threshold = 60, multithread=FALSE, quiet=FALSE,
                                 ranks = c("Root","Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species") ){
  # Print arguments for testing
  # Load the relevent db
  trainingSet <- readRDS(database)
  
  # get the sequences from the seqtab
  dna <- DNAStringSet(getSequences(seqtab)) # Create a DNAStringSet from the ASVs
  
  # Classify 
  ids <- DECIPHER::IdTaxa(dna, trainingSet, processors=1, threshold = threshold, verbose=!quiet, strand = "top") 
  
  # Get the filename of that db that we can use to name the output files
  db_name <- basename(database) %>% str_remove("_.*$")
  
  # Output plot of ids for each db
  pdf(paste0(qc_dir,"/", db_name,"idtaxa.pdf"), width = 11, height = 8 , paper="a4r")
  plot(ids)
  try(dev.off(), silent=TRUE)
  
  #Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
  tax <- t(sapply(ids, function(x) {
    taxa <- paste0(x$taxon,"_", x$confidence)
    taxa[startsWith(taxa, "unclassified_")] <- NA
    taxa
  })) %>%
    purrr::map(unlist) %>%
    stri_list2matrix(byrow=TRUE, fill=NA) %>%
    magrittr::set_colnames(ranks) %>%
    as.data.frame() %>%
    magrittr::set_rownames(getSequences(seqtab)) %T>%
    write.csv(paste0(qc_dir,"/", db_name,"idtaxa_results.csv")) %>%  #Write out logfile with confidence levels
    mutate_all(str_replace,pattern="(?:.(?!_))+$", replacement="") %>%
    magrittr::set_rownames(getSequences(seqtab)) 

  if(!is.null(output)){saveRDS(tax, output)}
  
  return(tax)
}

step_blast_tophit <- function(seqtab, output=NULL, qc_dir, database, identity = 97,  coverage=95, evalue=1e06,
                              max_target_seqs=5, max_hsp=5, 
                              ranks = c("Root","Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species"),
                              multithread=FALSE, quiet=FALSE){
  
  seqs <- taxreturn::char2DNAbin(getSequences(seqtab))
  names(seqs) <- getSequences(seqtab)
  
  # Get the filename of that db that we can use to name the output files
  db_name <- basename(database) %>% str_remove("_.*$")
  
  blast_spp <- blast_assign_species(query=seqs,db=database, identity=97, coverage=95, evalue=1e06,
                                    max_target_seqs=5, max_hsp=5, ranks=c("Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species") , delim=";") %>%
    dplyr::rename(blast_genus = Genus, blast_spp = Species) %>%
    dplyr::filter(!is.na(blast_spp))
  
  if(!is.null(output)){saveRDS(blast_spp, output)}
  return(blast_spp)
}

step_join_tax <- function(tax, blast_spp, output=NULL, propagate_tax=FALSE,
                          ranks = c("Root","Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species")){
  #Join together
  tax_blast <- tax %>%
    as_tibble(rownames = "OTU") %>%
    left_join(blast_spp , by="OTU") %>%
    mutate(Species = Species %>% str_replace_all("_", " ")) %>%
    dplyr::mutate(Species = case_when(
      is.na(Species) & Genus == blast_genus ~ blast_spp,
      !is.na(Species) ~ Species
    )) %>%
    dplyr::select(OTU, ranks) %>%
    column_to_rownames("OTU")
  
  if(propagate_tax){
    tax_blast <- tax_blast %>%
      seqateurs::na_to_unclassified() #Propagate high order ranks to unassigned ASVs
  }
  # Write taxonomy table for that db to disk
  if(!is.null(output)){saveRDS(tax_blast, output)}
  return(tax_blast)
}


# Utilities ---------------------------------------------------------------
setup_multithread <- function(multithread, create_future = FALSE, quiet = FALSE) {
  ncores <- future::availableCores()
  if (isTRUE(multithread)) {
    cores <- ncores - 1
    if (!quiet) {
      message("Multithreading with ", cores, " cores")
    }
  } else if (is.numeric(multithread) & multithread > 1) {
    cores <- multithread
  } else if (isFALSE(multithread) | multithread == 1) {
    cores <- 1
  } else (stop("Multithread must be a logical or numeric vector of the numbers of cores to use"))
  
  if (cores > ncores) {
    cores <- ncores
    warning("The value provided to multithread is higher than the number of cores, using ", cores, " cores instead")
  } else {
    if (!quiet & cores > 1) {message("Multithreading with ", cores, " cores")}
  }
  
  # Handle multithread types
  if(create_future & cores >1){
    future::plan(future::multiprocess, workers = cores)
  } else if(create_future & cores == 1){
    future::plan(future::sequential)
  } 
  return(cores)
}

remove_if_exists <- function(file, quiet=FALSE){
  if(file.exists(file)) {
    if(file.remove(file)) {
      if(!quiet) message("Overwriting file:", file)
    } else {
      stop("Failed to overwrite file:", file)
    }
  }
}
