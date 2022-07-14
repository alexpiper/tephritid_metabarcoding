library(targets)
library(tarchetypes)
source("R/functions.R")
options(tidyverse.quiet = TRUE)

# Load packages -----------------------------------------------------------
tar_option_set(packages = c(
  "devtools",
  "ggplot2",
  "gridExtra",
  "data.table",
  "tidyverse", 
  "stringdist",
  "patchwork",
  "vegan",
  "seqinr",
  "patchwork",
  "stringi",
  "phangorn",
  "magrittr",
  "galah",
  "phyloseq",
  "DECIPHER",
  "Biostrings",
  "ShortRead",
  "ggtree",
  "savR",
  "dada2",
  "ngsReports",
  "taxreturn",
  "seqateurs",
  "speedyseq"
), workspace_on_error = TRUE)


# Targets pipeline
list(
  # Track files
  tar_file(samdf_file, "sample_data/Sample_info.csv"),
  tar_file(params_file, "sample_data/loci_params.csv"),
  
  # Load sample tracking files
  tar_target(samdf, read_csv(samdf_file)),
  tar_target(params,read_csv(params_file)),
  
  # Create directories
  tar_target(create_dirs,  step_validate_folders(getwd())),
  
  # Look for sequencing reads
  tar_files(
    fastq_paths,
    purrr::map(list.dirs("data", recursive=FALSE),
                          list.files, pattern="_R1_", full.names = TRUE) %>%
      unlist() 
  ),

  # Check sequencing reads match sample sheet & Create temporary samdf file
  tar_file(
    make_temp_samdf1,
    {
      step_check_files(samdf, fastq_paths) %>%
      write_csv("temp/temp_samdf1.csv")
      return("temp/temp_samdf1.csv")
      }
  ),
  
  # Track temporary samdf
  tar_target(temp_samdf1, read_csv(make_temp_samdf1)),

# Per-Sample QC -----------------------------------------------------------
  tar_target(sample_qc, 
             temp_samdf1 %>% mutate(sample_qc = purrr::pmap(list(sample_id, fcid),
              .f = ~step_sample_qc(sample_id = ..1, fcid=..2, multithread=FALSE, quiet=TRUE)
              )),
            pattern = map(temp_samdf1), iteration = "vector"),
  
  tar_group_by(sample_qc_grouped, sample_qc, fcid),
  
  tar_target(multiqc, 
             sample_qc_grouped %>%
               group_by(fcid) %>%
               nest() %>%
               mutate(multi_qc = purrr::map(fcid, step_multiqc, quiet=FALSE)),
             pattern = map(sample_qc_grouped), iteration = "vector"),
  

# Sequencing run QC -------------------------------------------------------
  tar_target(seq_qc,
             sample_qc_grouped %>%
                group_by(fcid) %>%
                nest() %>%
                mutate(seq_qc = purrr::map(fcid, step_seq_qc)),
             pattern = map(sample_qc_grouped), iteration = "vector"),
  
  tar_target(switching_qc,
             sample_qc_grouped %>%
               group_by(fcid) %>%
               nest() %>%
               mutate(switching_qc = purrr::map(fcid, step_switching_calc, multithread=FALSE, quiet=TRUE)),
             pattern = map(sample_qc_grouped), iteration = "vector"),
             

# Demultiplex and trim primers --------------------------------------------
  tar_target(primer_trim,
            {
            process <- temp_samdf1 %>%
             mutate(primer_trim = purrr::pmap(dplyr::select(., sample_id, for_primer_seq, rev_primer_seq, pcr_primers, fcid),
                                     .f = ~step_primer_trim(sample_id = ..1, for_primer_seq=..2, rev_primer_seq=..3, pcr_primers = ..4,
                                                             input_dir = paste0("data/",..5), output_dir =  paste0("data/",..5,"/01_trimmed"),
                                                             qc_dir=paste0("output/logs/",..5), quiet = FALSE))) 
            #write_csv(process, "temp/primer_trim_res.csv") # This is overwriting - need a merge function
            # Return filepaths for tracking
            outF <- process$primer_trim %>%  
              bind_rows()%>%
              pull(fwd_out) 
            outR <- process$primer_trim %>%
              bind_rows()%>%
              pull(rev_out) 
            # Check for empty files
            outF <- outF[file.exists(outF)]
            outR <- outR[file.exists(outR)]
            # Return list of completed paths
            return(c(outF,outR))
            },
    pattern = map(temp_samdf1), format="file", iteration = "vector"),

 # tar_file(primer_trim_results, "temp/primer_trim_res.csv"),
  # Check sequencing reads match sample sheet & Create temporary samdf file
  tar_file(
    make_temp_samdf2,
    {
      temp_samdf1 %>%
      select(-where(is.list)) %>%
      step_demux_samdf() %>%
      step_add_params(params) %>%
      step_check_files(primer_trim) %>%
      write_csv("temp/temp_samdf2.csv")
      return("temp/temp_samdf2.csv")
    }
  ),

  # Track temporary samdf
  tar_target(temp_samdf2, read_csv(make_temp_samdf2)),
  
# Filter reads ------------------------------------------------------------
  tar_target(read_filter,
             {
             process <- temp_samdf2 %>%
             #  dplyr::slice(1:5)%>%
             mutate(read_filter = purrr::pmap(dplyr::select(., sample_id, fcid),
                   .f = ~step_read_filter(sample_id = ..1,
                   input_dir = paste0("data/",..2,"/01_trimmed/"), output_dir = paste0("data/",..2,"/02_filtered"),
                   maxEE = 1, truncLen = 150, rm.lowcomplex = 0,
                   quiet = FALSE)))
             #write_csv(process, "temp/filter_res.csv")
             # Return filepaths for tracking
             outF <- process$read_filter %>%  
               bind_rows()%>%
               pull(fwd_out) 
             outR <- process$read_filter %>%
               bind_rows()%>%
               pull(rev_out) 
             # Check for empty files
             outF <- outF[file.exists(outF)]
             outR <- outR[file.exists(outR)]
             # Return list of completed paths
             return(c(outF,outR))
             },                  
             pattern = map(temp_samdf2), format="file", iteration = "vector"),

  tar_file(
    make_temp_samdf3,
    {
      temp_samdf2 %>%
        select(-where(is.list)) %>%
        step_check_files(read_filter) %>%
        write_csv("temp/temp_samdf3.csv")
      return("temp/temp_samdf3.csv")
    }
  ),
  
  # Track temporary samdf
  tar_target(temp_samdf3, read_csv(make_temp_samdf3)),

## Pre-filtering quality plots ---------------------------------------------
#  #tar_target(prefilt_qualplots,
#  #          read_filter %>%
#  #          mutate(prefilt_qualplots = purrr::pmap(list(sample_id, fcid),
#  #                 .f = ~plot_read_quals(sample_id = ..1,
#  #                 input_dir = paste0("data/",..2,"/01_trimmed/"), truncLen=NULL, quiet = FALSE)
#  #            )),
#  #            pattern = map(read_filter), iteration = "vector"),
# #
#  ## Write out prefilt qualplots
#  #tar_target(write_prefilt_qualplots,
#  #          prefilt_qualplots %>% 
#  #            group_by(fcid) %>%
#  #            nest() %>%
#  #            purrr::pwalk(list(fcid, data),
#  #                         .f = ~{
#  #                           pdf(file=paste0("output/logs/",..1,"/prefilt_qualplots.pdf"), width = 11, height = 8 , paper="a4r")
#  #                           print(..2$prefilt_qualplots)
#  #                           try(dev.off(), silent=TRUE)
#  #                           
#  #                         })),
#
## Post-filtering quality plots --------------------------------------------
#  #tar_target(postfilt_qualplots,
#  #          read_filter %>%
#  #            mutate(postfilt_qualplots = purrr::pmap(list(sample_id, fcid),
#  #                 .f = ~plot_read_quals(sample_id = ..1,
#  #                 input_dir = paste0("data/",..2,"/02_filtered/"), truncLen=NULL, quiet = FALSE)
#  #           )),
#  #            pattern = map(read_filter), iteration = "vector"),
#  # 
#  ## Write out postfilt qualplots
#  #tar_target(write_postfilt_qualplots,
#  #            postfilt_qualplots %>% 
#  #              group_by(fcid) %>%
#  #              nest() %>%
#  #              purrr::pwalk(list(fcid, data),
#  #                           .f = ~{
#  #                             pdf(file=paste0("output/logs/",..1,"/postfilt_qualplots.pdf"), width = 11, height = 8 , paper="a4r")
#  #                             print(..2$postfilt_qualplots)
#  #                             try(dev.off(), silent=TRUE)
#  #                             
#  #                           })),
# 
#
# Infer sequence variants with DADA2 --------------------------------------

  # Create new samdf from output of previous step? Edit check function to take the 
  tar_group_by(temp_samdf3_grouped, temp_samdf3, fcid),
 
  # How to make it just redo one of the dadas if only one runs filtered file changed changed?
  tar_target(dada,{
             process <- temp_samdf3_grouped %>%
             group_by(fcid) %>%
             nest() %>%
             mutate(dada2 = purrr::map(fcid,
                                        .f = ~step_dada2(fcid = .x,
                                                         input_dir = paste0("data/",.x,"/02_filtered"),
                                                         output = paste0("output/rds/",.x,"_seqtab.rds"),
                                                         qc_dir = paste0("output/logs/",.x),
                                                         quiet = FALSE)
             )) %>%
               unnest(data)
             #write_csv(process, "temp/filter_res.csv")
             out <- paste0("output/rds/",unique(temp_samdf3_grouped$fcid),"_seqtab.rds")
             return(out)
             },
             pattern = map(temp_samdf3_grouped), format="file", iteration = "vector"),
  

  tar_file(
    make_temp_samdf4,
    {
      temp_samdf3 %>%
        select(-where(is.list)) %>%
        left_join(enframe(dada, name=NULL, value="output") %>%
                  mutate(fcid = basename(output) %>% str_remove("_seqtab.rds"))) %>%
        write_csv("temp/temp_samdf4.csv")
      return("temp/temp_samdf4.csv")
    }
  ),
  
  # Track temporary samdf
  tar_target(temp_samdf4, read_csv(make_temp_samdf4)),

##  Merge infered variants from each run and subset to target loci ---------
#
 tar_target(subset_seqtab,
            process <- temp_samdf4 %>%
             ungroup() %>%
             group_by(pcr_primers) %>%
             nest() %>%
             mutate(subset_seqtab = purrr::map(pcr_primers, 
                   .f = ~{
                   #seqtabs <- list.files("output/rds/", pattern="seqtab.rds", full.names = TRUE)
                   if(length(dada) > 1){
                   st.all <- mergeSequenceTables(tables=dada)
                   } else if(length(dada) == 1) {
                   st.all <- readRDS(dada)
                   }
                   st.all <- st.all[str_detect(rownames(st.all), .x),]
                   st.all <- st.all[,colSums(st.all) > 0]
                   return(st.all)
                   }))%>%
              unnest(data)),
 

##  Filter ASV's per locus -------------------------------------------------
#
  tar_target(filtered_seqtab,
             subset_seqtab %>%
              group_by(pcr_primers, subset_seqtab, exp_length, phmm, coding, genetic_code, for_primer_seq, rev_primer_seq) %>%
              nest() %>%
              ungroup()%>%
              mutate(filtered_seqtab = purrr::pmap(dplyr::select(.,pcr_primers, subset_seqtab, exp_length, phmm, coding, genetic_code, for_primer_seq, rev_primer_seq),
                  .f = ~step_filter_asvs(
                  seqtab = ..2,
                  output = paste0("output/rds/",..1,"_seqtab.cleaned.rds"),
                  qc_dir = "output/logs/",
                  min_length = ..3-10,
                  max_length = ..3+10,
                  phmm = ..4,
                  check_frame = ..5,
                  genetic_code = ..6,
                  primers = c(..7, ..8),
                  multithread = FALSE, 
                  quiet = FALSE)
    )) %>%
    unnest(data)),
  

## Merge all loci into final seqtab ----------------------------------------
  tar_target(merged_seqtab,
             filtered_seqtab %>%
              nest(data=everything()) %>%
              mutate(final_seqtab = purrr::map(data, ~{
                seqtabs <- fs::dir_ls("output/rds/", glob="*.cleaned.rds")
                seqtabs <- seqtabs[seqtabs %>%
                                     purrr::map_lgl(function(y){
                                       any(str_detect(y, unique(.x$pcr_primers)))
                                     })]
                if(length(seqtabs) > 1){
                  st.all <- mergeSequenceTables(tables=seqtabs)
                } else if(length(seqtabs) == 1) {
                  st.all <- readRDS(seqtabs)
                }
                saveRDS(st.all, "output/rds/seqtab_final.rds")
                return(TRUE)
              })) %>%
    unnest(data)),
  
  # Create grouped seqtab
  tar_group_by(merged_seqtab_grouped, merged_seqtab, target_gene),

# Assign taxonomy ---------------------------------------------------------

 # # Track the taxonomy files
 # tar_file(blast_db,
 #          params  %>%
 #              pull(blast_db) %>%
 #              unique(),
 # ),
 # tar_file(ref_db,
 #          params  %>%
 #            pull(ref_db) %>%
 #            unique(),
 # ),

# Shouldnt need to add the parames earlier, should just be able to join here

# GOT UP TO HERE WITH TRACKING EDITS

 tar_target(tax_idtaxa,
            {
           merged_seqtab_grouped %>%
           #dplyr::select(-ref_db) %>% # Can just join with the params here
           #left_join(params %>% dplyr::select(pcr_primers, ref_db)) %>%
           group_by(target_gene, filtered_seqtab, ref_db) %>%
           nest() %>%
           mutate(idtaxa = purrr::pmap(list(target_gene, filtered_seqtab, ref_db),
                                       .f = ~step_assign_taxonomy(
                                         seqtab = ..2$filtered_seqtab,
                                         database = ..3,
                                         ranks = c("Root","Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species") ,
                                         output = paste0("output/rds/",..1,"_taxtab.rds"),
                                         qc_dir = "output/logs/",
                                         threshold = 60,
                                         multithread = FALSE, 
                                         quiet = FALSE)
           ))%>%
           unnest(data)
              },
     pattern = map(merged_seqtab_grouped), iteration = "vector"),
# NEed to return character vector instead - Or can i just track the output vectors?
 
 tar_target(tax_blast,
            tax_idtaxa %>%
           group_by(target_gene, filtered_seqtab, blast_db) %>%
           nest() %>% 
           mutate(blast = purrr::pmap(list(target_gene, filtered_seqtab, blast_db),
                                      .f = ~step_blast_tophit(
                                        seqtab = ..2$filtered_seqtab,
                                        database = ..3,
                                        ranks = c("Root","Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species") ,
                                        output = paste0("output/rds/",..1,"_blast.rds"),
                                        qc_dir = "output/logs/",
                                        identity = 97,  
                                        coverage=95,
                                        evalue=1e06,
                                        max_target_seqs=5,
                                        max_hsp=5, 
                                        multithread = FALSE, 
                                        quiet = FALSE)
           ))%>%
             unnest(data),
 pattern = map(tax_idtaxa), iteration = "vector"),
 
#
## Aggregate joint tax -----------------------------------------------------
  tar_target(tax_blast_combined, tax_blast %>% ungroup()),
  
  tar_target(joint_tax,
             tax_blast_combined %>%
             group_by(target_gene, idtaxa, blast) %>%
             nest() %>% 
             mutate(joint_tax = purrr::pmap(list(target_gene, idtaxa, blast),
                                            .f = ~step_join_tax(
                                              tax = ..2,
                                              blast_spp = ..3,
                                              ranks = c("Root","Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species") ,
                                              output = paste0("output/rds/",..1,"_taxblast.rds"),
                                              propagate_tax = TRUE)
             )) %>%
             unnest(data)),
  
  
  # Look for taxtabs
  #tar_files(taxtabs, list.files("output/rds/", pattern="_taxblast.rds", full.names = TRUE)),
  
  tar_target(merged_tax,
             joint_tax %>%
             nest(data=everything())%>%
             mutate(merged_tax = purrr::map(data,
                                           .f = ~{
                                             taxtabs <- list.files("output/rds/", pattern="_taxblast.rds", full.names = TRUE)
                                             taxtabs <- taxtabs[taxtabs %>%
                                                                  purrr::map_lgl(function(y){
                                                                    any(str_detect(y, unique(.x$target_gene)))
                                                                  })]
                                             tax <- taxtabs %>%
                                               purrr::map(readRDS) %>%
                                               bind_rows() %>%
                                               as.matrix()
                                             saveRDS(tax, "output/rds/final_tax.rds")
                                             return(tax)
                                           })) %>%
            unnest(data)),

# Create phyloseq object --------------------------------------------------
  tar_file(seqtab_path, "output/rds/seqtab_final.rds"),
  tar_file(taxtab_path, "output/rds/final_tax.rds"),
  tar_file(samdf_path, "sample_data/Sample_info_updated.csv"),
  
  tar_target(ps,
             step_phyloseq(
              seqtab_path = seqtab_path,
              taxtab_path = taxtab_path,
              samdf_path = samdf_path,
              seqs_path=NULL,
              phy_path=NULL)
  ),

## Output unfiltered results -----------------------------------------------
  tar_file(ps_path, "output/rds/ps.rds"),
  tar_file(ps_summary_path, "output/results/unfiltered"), # Do i need to have all the outputs as files to track?
  
  tar_target(
    ps_output, {
      ps %>%
        saveRDS(ps_path)
    }
  ),
  tar_target(ps_summary,
             ps %>%
              step_output_summary(out_dir=ps_summary_path, type="unfiltered")
  ),
  
  # Subset taxa in phyloseq object
  tar_target(ps0,
             ps %>%
             subset_taxa(
               Phylum == "Arthropoda"
             ) %>%
            filter_taxa(function(x) mean(x) > 0, TRUE)
  ),
  
  # Filter low abundance samples
  tar_target(ps1,
             ps0 %>%
             step_filter_phyloseq(
               min_reads=1000, 
               quiet=FALSE
             )
  ),

## Output filtered results -------------------------------------------------
  tar_file(ps_filt_path, "output/rds/ps_filtered.rds"),
  tar_file(ps_filt_summary_path, "output/results/filtered"),
  
  tar_target(
   ps_filt_output, {
     ps1 %>%
       saveRDS(ps_filt_path)
   }
  ),
  tar_target(ps_filt_summary,
             ps1  %>%
             step_output_summary(out_dir=ps_filt_summary_path, type="filtered")
  ),
  
  tar_file(ps_imap_path, "output/results/final"),

## Output filtered results in imappests format -----------------------------
  tar_target(ps_imap_output,
             ps1  %>%
             step_output_imap(out_dir=ps_imap_path)
  )
)

