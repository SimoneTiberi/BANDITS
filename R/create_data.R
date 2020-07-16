#' Create a 'BANDITS_data' object
#'
#' \code{create_data} imports the equivalence classes and create a 'BANDITS_data' object.
#' 
#' 
#' @param salmon_or_kallisto a character string indicating the input data: 'salmon' or 'kallisto'.
#' @param gene_to_transcript a matrix or data.frame with a list of gene-to-transcript correspondances.
#' The first column represents the gene id, while the second one contains the transcript id.
#' @param salmon_path_to_eq_classes (for salmon input only) a vector of length equals to the number of samples: 
#' each element indicates the path to the equivalence classes of the respective sample (computed by salmon).
#' @param kallisto_equiv_classes (for kallisto input only) a vector of length equals to the number of samples: 
#' each element indicates the path to the equivalence classes ('.ec' files) of the respective sample (computed by kallisto).
#' @param kallisto_equiv_counts (for kallisto input only) a vector of length equals to the number of samples: 
#' each element indicates the path to the counts of the equivalence classes ('.tsv' files) of the respective sample (computed by kallisto).
#' @param kallisto_counts (for kallisto input only) a matrix or data.frame, 
#' with 1 column per sample and 1 row per transcript, 
#' containing the estimated abundances for each transcript in each sample, computed by kallisto.
#' The matrix must be unfiltered and the order or rows must be unchanged.
#' @param eff_len a vector containing the effective length of transcripts; the vector names indicate the transcript ids.
#' Ideally, created via \code{\link{eff_len_compute}}.
#' @param n_cores the number of cores to parallelize the tasks on.
#' It is highly suggested to use at least one core per sample (default if not specificied by the user).
#' @param transcripts_to_keep a vector containing the list of transcripts to keep.
#' Ideally, created via \code{\link{filter_transcripts}}.
#' @param max_genes_per_group an integer number specifying the maximum number of genes that each group can contain.
#' When equivalence classes contain transcripts from distinct genes, these genes are analyzed together.
#' For computational reasons, 'max_genes_per_group' sets a limit to the number of genes that each group can contain.
#' 
#' @return A \code{\linkS4class{BANDITS_data}} object.
#' @examples
#' # specify the directory of the internal data:
#' data_dir = system.file("extdata", package = "BANDITS")
#' 
#' # load gene_to_transcript matching:
#' data("gene_tr_id", package = "BANDITS")
#' 
#' # Specify the directory of the transcript level estimated counts.
#' sample_names = paste0("sample", seq_len(4))
#' quant_files = file.path(data_dir, "STAR-salmon", sample_names, "quant.sf")
#' 
#' # Load the transcript level estimated counts via tximport:
#' library(tximport)
#' txi = tximport(files = quant_files, type = "salmon", txOut = TRUE)
#' counts = txi$counts
#' 
#' # Optional (recommended): transcript pre-filtering
#' transcripts_to_keep = filter_transcripts(gene_to_transcript = gene_tr_id,
#'                                          transcript_counts = counts,
#'                                          min_transcript_proportion = 0.01,
#'                                          min_transcript_counts = 10,
#'                                          min_gene_counts = 20)
#' 
#' # compute the Median estimated effective length for each transcript:
#' eff_len = eff_len_compute(x_eff_len = txi$length)
#' 
#' # specify the path to the equivalence classes:
#' equiv_classes_files = file.path(data_dir, "STAR-salmon", sample_names, "aux_info", "eq_classes.txt")
#' 
#' # create data from 'salmon' and filter internally lowly abundant transcripts:
#' input_data = create_data(salmon_or_kallisto = "salmon",
#'                          gene_to_transcript = gene_tr_id,
#'                          salmon_path_to_eq_classes = equiv_classes_files,
#'                          eff_len = eff_len, 
#'                          n_cores = 2,
#'                          transcripts_to_keep = transcripts_to_keep)
#' input_data
#' 
#' # create data from 'kallisto' and filter internally lowly abundant transcripts:
#' kallisto_equiv_classes = file.path(data_dir, "kallisto", sample_names, "pseudoalignments.ec")
#' kallisto_equiv_counts  = file.path(data_dir, "kallisto", sample_names, "pseudoalignments.tsv")
#' 
#' input_data_2 = create_data(salmon_or_kallisto = "kallisto",
#'                           gene_to_transcript = gene_tr_id,
#'                           kallisto_equiv_classes = kallisto_equiv_classes,
#'                           kallisto_equiv_counts = kallisto_equiv_counts,
#'                           kallisto_counts = counts,
#'                           eff_len = eff_len, n_cores = 2,
#'                           transcripts_to_keep = transcripts_to_keep)
#' input_data_2
#' 
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#' 
#' @seealso \code{\link{eff_len_compute}}, \code{\link{filter_transcripts}}, \code{\link{filter_genes}}, \code{\linkS4class{BANDITS_data}}
#' 
#' @export
create_data = function(salmon_or_kallisto,
                       gene_to_transcript,
                       salmon_path_to_eq_classes = NULL,
                       kallisto_equiv_classes = NULL,
                       kallisto_equiv_counts = NULL,
                       kallisto_counts = NULL,
                       eff_len,
                       n_cores = NULL, 
                       transcripts_to_keep = NULL,
                       max_genes_per_group = 50){
  # wrapper to call the correct function: salmon_create_data (identical to the one already created) or kallisto_create_data (the one below)
  cond_salmon_or_kallisto = (length(salmon_or_kallisto) == 1) & is.character(salmon_or_kallisto) & (salmon_or_kallisto %in% c("salmon", "kallisto"))
  
  if( !cond_salmon_or_kallisto ){
    message("'salmon_or_kallisto' must be a character string indicating the input data as: 'salmon' or 'kallisto'")
    return(NULL)
  }
  
  # check that gene_to_transcript is a matrix or data.frame object
  if( !is.data.frame(gene_to_transcript) & !is.matrix(gene_to_transcript)  ){
    message("'gene_to_transcript' must be a matrix or data.frame")
    return(NULL)
  }
  
  if( ncol(gene_to_transcript) != 2 ){
    message("'gene_to_transcript' must be a 2 column matrix or data.frame")
    return(NULL)
  }
  
  if( !all( names(eff_len) %in% gene_to_transcript[,2])  ){
    message("All transcript names in 'names(eff_len)' must be in 'gene_to_transcript[,2]'")
    return(NULL)
  }
  
  if(max_genes_per_group > 100){
    message("'max_genes_per_group' can be at most 100")
    return(NULL)
  }
  
  # define the number of parellel cores (if left unspecified):
  if(is.null(n_cores)){ 
    n_cores = ifelse(salmon_or_kallisto == "salmon", length(salmon_path_to_eq_classes), length(kallisto_equiv_classes))
  }
  
  # initialize parallel cores (if n_cores > 1)
  if( n_cores > 1.5){ # if n_cores > 1, I use parallel computing tools
    suppressWarnings({
      cl = makeCluster(n_cores, setup_strategy = "sequential")
    })
  }
  
  # calculate the separator (1 or more _)
  if(is.null(transcripts_to_keep)){
    sep = "_"
    while(TRUE){
      Tr_id = as.character(gene_to_transcript[,2]) 
      cond = sum( vapply(Tr_id, function(x){  grepl(sep, x, fixed=TRUE) }, FUN.VALUE = logical(1)) ) == 0
      # sum(sapply(Tr_id, function(x){  grepl(sep, x, fixed=TRUE) })) == 0
      if(cond){
        break
      }
      sep = paste0(sep, "_")
    }
  }else{
    sep = "_"
    while(TRUE){
      cond = sum( vapply(transcripts_to_keep, function(x){  grepl(sep, x, fixed=TRUE) }, FUN.VALUE = logical(1)) ) == 0
      # sum(sapply(transcripts_to_keep, function(x){  grepl(sep, x, fixed=TRUE) })) == 0
      if(cond){
        break
      }
      sep = paste0(sep, "_")
    }
  }
  
  # load the data:
  if(salmon_or_kallisto == "salmon"){ # load salmon equivalence classes:
    x = load_salmon_data(sep = sep,
                         path_to_eq_classes = salmon_path_to_eq_classes,
                         n_cores = n_cores, cl = cl,
                         transcripts_to_keep = transcripts_to_keep)
  }else{ # load kallisto equivalence classes:
    x = load_kallisto_data(sep = sep,
                           kallisto_equiv_classes = kallisto_equiv_classes,
                           kallisto_equiv_counts = kallisto_equiv_counts,
                           kallisto_counts = kallisto_counts,
                           n_cores = n_cores, cl = cl, 
                           transcripts_to_keep = transcripts_to_keep)
  }
  
  if(is.null(x)){ # if the object is NULL, return it
    return(x)
  }
  
  message("Data has been loaded")
  
  N = length(x)
  # From the truth table I take the Transcript_id and associated Gene_id.
  Gene_id = as.character(gene_to_transcript[,1]); Tr_id = as.character(gene_to_transcript[,2]) 
  
  # merge the ids of all classes here:
  all_classes = unique( unlist( lapply(x, function(y){y$class_ids}) ) )
  # unique( unlist( sapply(x, function(y){y$class_ids}) ) )
  
  # I turn the transcript id from a single character into a vector
  all_classes_vector = vapply( all_classes, strsplit, split = sep, fixed = TRUE, FUN.VALUE = list(1) )
  # sapply( all_classes, strsplit, split = sep, fixed = TRUE )
  
  # match the class id of each sample to the long vector containing all classes ids.
  match_classes = lapply(x, function(u) match(u$class_ids, all_classes) )
  # match_classes is a list, every element of the list refers to a sample.
  # match_classes may now have some NAs.
  
  # I find the counts of each equiv class in all the samples.
  all_counts = matrix(0, nrow = length(all_classes), ncol = N)
  # I set the matrix to 0: in case of no matching (class not present in a sample)
  for(i in seq_len(N) ){ # super-fast:
    match_NA = is.na(match_classes[[i]])
    all_counts[match_classes[[i]][match_NA==FALSE],i] = x[[i]]$counts[match_NA==FALSE]
    # here I filter the classes which were filtered out when filtering the transcripts.
  }
  
  # nr of transcripts per class:
  n = vapply(all_classes_vector, length, FUN.VALUE = integer(1))
  # sapply(all_classes_vector, length)
  # create a data.frame structure keeping the info of what transcritps belong to what classes:
  df_all_classes = data.frame( Class_id=rep(seq_along(n), n), Tr_id = unlist(all_classes_vector) )
  rownames(df_all_classes) = c()
  # then match equivalence classes to their genes (see readDGE)
  # all_classes_vector is a list: every element of the list corresponds to a class.
  df_all_classes$Gene_id = Gene_id[match(df_all_classes$Tr_id, Tr_id)]
  
  genes_in_classes = split(df_all_classes$Gene_id, df_all_classes$Class_id)
  genes_in_classes = lapply(genes_in_classes, unique)
  
  ############################################################################################################################################################
  # 1) Consider ALL genes: first separate genes to be modelled alone from genes to be modelled together.
  ############################################################################################################################################################
  genes_SELECTED = unique( df_all_classes$Gene_id ) #  All genes
  N_genes = length(genes_SELECTED)
  
  n = vapply(genes_in_classes, length, FUN.VALUE = integer(1))
  # sapply(genes_in_classes, length)
  
  ## I need to consider classes with 1 gene only...easy: "cond = sapply(genes_in_classes, length) == 1"
  cond_1_gene_per_class = ( n == 1 )
  # sapply(genes_in_classes, length) == 1 # But I also need to make sure that the genes respecting "cond" does not happear in other classes!
  More_genes_in_classes = unique(unlist(genes_in_classes[cond_1_gene_per_class ==FALSE])) # list the genes happearing together in at least 1 class
  
  # 1 - mean(cond_1_gene_per_class) # mean of equiv classes with tr from > 1 gene.
  # 16 % before filterign, 17% after filtering.
  # sum(all_counts[!cond_1_gene_per_class,])/sum(all_counts)
  # mean counts from these classes
  # 11 % before filtering, 10% after filtering.
  
  # check what classes have at least 1 gene happearing with other genes in at least 1 class:
  # in other words, check whether it's on the "More_genes_in_classes" list or not.
  
  df_tmp = data.frame( Class_id=rep(seq_along(n), n), Gene_id = unlist(genes_in_classes))
  df_tmp$More_genes_in_classes = df_tmp$Gene_id %in% More_genes_in_classes
  cond_1_gene_per_class_FINAL = split(df_tmp$More_genes_in_classes, df_tmp$Class_id)
  cond_1_gene_per_class_FINAL = !vapply(cond_1_gene_per_class_FINAL, any, FUN.VALUE = logical(1))
  # !sapply(cond_1_gene_per_class_FINAL, any)
  # 0.3 seconds
  
  genes_names_Unique = as.character( genes_in_classes[cond_1_gene_per_class_FINAL] )
  # genes to be modelled alone:
  genes_SELECTED_Unique   = unique( genes_names_Unique )
  # genes to be modelled together (not in the previous list!):
  genes_SELECTED_Together = genes_SELECTED[genes_SELECTED %in% genes_SELECTED_Unique ==FALSE]
  
  N_genes_Unique = length(genes_SELECTED_Unique)
  # genes to be modelled together (not in the previous list!):
  genes_SELECTED_Together = unique( df_all_classes$Gene_id[df_all_classes$Gene_id %in% genes_SELECTED_Unique ==FALSE] );
  N_genes_Together = length(genes_SELECTED_Together)
  
  if(N_genes_Unique + N_genes_Together == 0){
    message("0 genes to be analyzed")
    return(NULL)
  }
  
  ############################################################################################################################################################
  # 2) Prepare data for the Unique genes:
  ############################################################################################################################################################
  # ADD CONSTRAINT: if(N_genes_Unique > 0)
  
  # Collect, for each gene, the classes associated to it.
  classes_split_per_gene_Unique = split(all_classes_vector[cond_1_gene_per_class_FINAL], genes_names_Unique)
  
  # Collect the counts associated to each class in every gene.
  counts_split_per_gene_Unique  = split(data.frame(all_counts[cond_1_gene_per_class_FINAL,]), genes_names_Unique)
  
  # Finally, for each gene, consider the full list of transcripts (to define ... and ...)
  # and make a matrix of 0, 1 indicating, for each class, what transcripts they have.
  Transcripts_per_gene_Unique   = lapply(classes_split_per_gene_Unique, function(x){ unique(unlist(x)) })
  N_transcripts_per_gene_Unique = vapply(Transcripts_per_gene_Unique, length, FUN.VALUE = integer(1))
  # sapply(Transcripts_per_gene_Unique, length)
  
  N_genes_Unique = length(counts_split_per_gene_Unique);
  
  #Error in classes_split_per_gene_Unique[[x]] : subscript out of bounds
  # length(classes_split_per_gene_Unique) == length(counts_split_per_gene_Unique) # TRUE
  
  classes_Unique = lapply(seq_len(N_genes_Unique), function(x){
    m = vapply(classes_split_per_gene_Unique[[x]], function(y){ Transcripts_per_gene_Unique[[x]] %in% y }, FUN.VALUE = logical( length( Transcripts_per_gene_Unique[[x]] ) ))
    # sapply(classes_split_per_gene_Unique[[x]], function(y){ Transcripts_per_gene_Unique[[x]] %in% y })
    ifelse(m, 1, 0)
  })
  # 0.7 seconds (Mac)
  
  ################################################################################################
  # Associate the Median eff length of transcripts:
  ################################################################################################
  n = vapply(Transcripts_per_gene_Unique, length, FUN.VALUE = integer(1))
  # sapply(Transcripts_per_gene_Unique, length)
  df_eff_len_Unique = data.frame( Gene_id=rep(seq_along(n), n), Tr_id = unlist(Transcripts_per_gene_Unique))
  df_eff_len_Unique$eff_len = eff_len[ match(df_eff_len_Unique$Tr_id, names(eff_len))  ]
  
  eff_len_tr_Unique = split(df_eff_len_Unique$eff_len, df_eff_len_Unique$Gene_id)
  
  genes_SELECTED_Unique = names(Transcripts_per_gene_Unique)
  
  ######################################################################################################
  # TOGETHER Genes:
  ######################################################################################################
  # ADD CONSTRAINT: if(N_genes_Together > 0)
  
  all_counts_Together         = all_counts[cond_1_gene_per_class_FINAL==FALSE,]
  all_classes_vector_Together = all_classes_vector[cond_1_gene_per_class_FINAL==FALSE]
  genes_in_classes_Together   = genes_in_classes[cond_1_gene_per_class_FINAL==FALSE]
  n_genes_per_class           = vapply(genes_in_classes_Together, length, FUN.VALUE = integer(1))
  # sapply(genes_in_classes_Together, length)
  
  # First, I need to separate classes corresponding to 1 gene only to classes corresponding to >1 gene...DONE
  # Now I need to group together the information about genes modelled together...they could be highly mixed!
  # maybe I can use a similar approach to the one used to build the classes from the transcripts, although a 2,000 * 17,000 matrix would be quite big...
  
  ######################################################################################################
  genes_in_classes_vector_Together = vapply( genes_in_classes_Together, paste, collapse = sep, FUN.VALUE = character(1) )
  # sapply( genes_in_classes_Together, paste, collapse = sep )
  names( genes_in_classes_vector_Together ) = genes_in_classes_vector_Together
  
  classes_split_per_gene_Together = split( all_classes_vector_Together,
                                           genes_in_classes_vector_Together )
  
  counts_split_per_gene_Together =  split( data.frame(all_counts_Together),
                                           genes_in_classes_vector_Together )
  
  # I look for what classes each gene appers in and record whether it happears uniquely or not.
  genes_in_classes_split_per_gene_Together = strsplit(names(classes_split_per_gene_Together), split = sep, fixed = TRUE )
  n_genes = vapply(genes_in_classes_split_per_gene_Together, length, FUN.VALUE = integer(1))
  # sapply(genes_in_classes_split_per_gene_Together, length)
  
  ######################################################################################################
  # I make group of genes to be modelled together and make a correspondance with the classes in classes_split_per_gene_Together
  ######################################################################################################
  GROUPs_of_genes = list()
  classes_associated_to_GROUPs = list()
  genes_included = c()
  g_id = 0
  # it must be a loop (g_id does not always increase).
  for(i in seq_len(N_genes_Together) ){
    if(genes_SELECTED_Together[i] %in% genes_included == FALSE){
      g_id = g_id + 1 # I create a new group of genes.
      classes_associated_to_GROUPs[[g_id]] = which( vapply(genes_in_classes_split_per_gene_Together, function(x){ genes_SELECTED_Together[i] %in% x}, FUN.VALUE = logical(1)) )
      # which( sapply(genes_in_classes_split_per_gene_Together, function(x){ genes_SELECTED_Together[i] %in% x}) )
      GROUPs_of_genes[[g_id]] = unique( unlist(genes_in_classes_split_per_gene_Together[classes_associated_to_GROUPs[[g_id]]]) ) # Genes associated to i-th gene
      j = 1
      genes_included = c(genes_included, genes_SELECTED_Together[i])
      while(j <= length(GROUPs_of_genes[[g_id]]) ){ # I loop on the other genes and repeat the operation
        gene = GROUPs_of_genes[[g_id]][j]
        if( !(gene %in% genes_included) ){
          classes_associated_to_GROUPs[[g_id]] = unique( c(classes_associated_to_GROUPs[[g_id]],
                                                           which( vapply(genes_in_classes_split_per_gene_Together, 
                                                                         function(x){ gene %in% x}, FUN.VALUE = logical(1)) ) ))
          GROUPs_of_genes[[g_id]] = unique( c(GROUPs_of_genes[[g_id]], 
                                              unlist(genes_in_classes_split_per_gene_Together[classes_associated_to_GROUPs[[g_id]]]) ) )
          # Genes associated to i-th gene
          genes_included = c(genes_included, gene)
        }
        j = j+1
      }
    }
  }
  
  ######################################################################################################
  # check if a Group has too many genes and split it START:
  ######################################################################################################
  n_genes_per_group = vapply(GROUPs_of_genes, length, FUN.VALUE = integer(1))
  # sapply(GROUPs_of_genes, length)
  bigGroup = which( n_genes_per_group > max_genes_per_group)
  if(length(bigGroup) > 0){ # if at least 1 group to be split
    for(p in bigGroup){
      message("One group of genes has ", n_genes_per_group[p], " genes")
      message("Splitting the group in ",  n_genes_per_group[p], " groups of individual genes")
      
      classes_bigGroup = classes_associated_to_GROUPs[[p]]
      ################################################################################################
      # Split genes and define their classes:
      ################################################################################################
      classes_split_per_gene_Together = split( all_classes_vector_Together,
                                               genes_in_classes_vector_Together )
      counts_split_per_gene_Together =  split( data.frame(all_counts_Together),
                                               genes_in_classes_vector_Together )
      Gene   = unlist(genes_in_classes_split_per_gene_Together[classes_bigGroup])
      
      classes_tmp = classes_split_per_gene_Together[classes_bigGroup]
      Ngenes = vapply(genes_in_classes_split_per_gene_Together[classes_bigGroup], length, FUN.VALUE = integer(1))
      # sapply(genes_in_classes_split_per_gene_Together[classes_bigGroup], length)
      Class  = rep(classes_tmp, Ngenes)
      Class = unname(Class, force = TRUE) # remove names of top lists of Class
      Class = lapply(Class, unname, force = TRUE) # remove names of each list in Class
      # REMOVE NAMES OF OBJECTS: otherwise it'll crash (names too long!).
      
      Counts = rep(counts_split_per_gene_Together[classes_bigGroup], vapply(genes_in_classes_split_per_gene_Together[classes_bigGroup], length, FUN.VALUE = integer(1)))
      # sapply(genes_in_classes_split_per_gene_Together[classes_bigGroup], length))
      Counts = unname(Counts, force = TRUE)
      # REMOVE NAMES OF OBJECTS: otherwise it'll crash (names too long!).
      
      class = do.call(c, Class)
      gene = rep(Gene, vapply(Class, length, FUN.VALUE = integer(1)) )
      # sapply(Class, length) )
      counts = do.call(rbind, Counts)
      
      # Split counts and classes by their gene name.
      class_byGene = split(class, gene)
      counts_byGene = split(counts, gene)
      genes_bigGroup = names(class_byGene)
      ################################################################################################
      # REMOVE THE TRANSCRIPTS THAT DO NOT BELONG THE GENE WE CONSIDER:
      ################################################################################################
      df = data.frame(Gene_id = Gene_id[Gene_id %in% genes_bigGroup], Tr_id = Tr_id[Gene_id %in% genes_bigGroup] )
      Tr_perGene_bigGroup = split(df$Tr_id, df$Gene_id)
      Tr_perGene_bigGroup = lapply(Tr_perGene_bigGroup, as.character)
      # it contains, for each gene, the list of transcripts associated to it.
      
      # the order of the gene names is identical because split orders them alphabetically in both cases.
      
      classes_split_bigGroup = lapply(seq_along(class_byGene), function(x){
        X = class_byGene[[x]]
        lapply(X, function(y){ 
          y = unlist(y)
          y[y %in% Tr_perGene_bigGroup[[x]] ] 
        })
      })
      ################################################################################################
      # Merge Identical classes and add up corresponding counts:
      ################################################################################################
      # put transcript names of classes together as tr1_tr2
      classes_split_bigGroup_num = lapply(classes_split_bigGroup, function(y){
        # as.numeric(factor (sapply(y, function(yy) paste(sort(yy), collapse=sep) )) ) 
        as.numeric(factor (vapply(y, function(yy) paste(sort(yy), collapse=sep), FUN.VALUE = character(1) ) ) ) 
      })
      
      DUPS = lapply(classes_split_bigGroup_num,  duplicated)
      classes_split_bigGroup_Unique = lapply(seq_along(classes_split_bigGroup), function(id){
        classes_split_bigGroup[[id]][ DUPS[[id]] == FALSE ]
      })
      
      classes_split_bigGroup_Unique_num = lapply(seq_along(classes_split_bigGroup_num), function(id){
        classes_split_bigGroup_num[[id]][ DUPS[[id]] == FALSE ]
      })
      
      # Select the counts for the unique classes:
      counts_byGene_unique = lapply(seq_along(counts_byGene), function(id){
        counts_byGene[[id]][ DUPS[[id]] == FALSE, ]
      })
      
      # Select the counts for the duplicated classes:
      counts_byGene_DUPS = lapply(seq_along(counts_byGene), function(id){
        counts_byGene[[id]][ DUPS[[id]], ]
      })
      
      counts_byGene_final = lapply(seq_along(counts_byGene), function(id){
        X = counts_byGene_unique[[id]]
        for(i in seq_len(sum(DUPS[[id]]==FALSE)) ){
          match_dups = classes_split_bigGroup_num[[id]][DUPS[[id]]] == classes_split_bigGroup_Unique_num[[id]][i] # I look for the matching between the unique classes and duplicatd ones (eliminated)
          if(sum(match_dups) > 0){
            X[i,] = X[i,] + colSums(counts_byGene_DUPS[[id]][match_dups,]) # I add the counts of the corresponding duplicated classes
          }
        }
        X
      })
      ################################################################################################
      # Finally, for each gene, consider the full list of transcripts (to define ... and ...)
      # and make a matrix of 0, 1 indicating, for each class, what transcripts they have.
      ################################################################################################
      Transcripts_per_gene_bigGroup   = lapply(classes_split_bigGroup_Unique, function(x){ unique(unlist(x)) })
      N_transcripts_per_gene_bigGroup = vapply(Transcripts_per_gene_bigGroup, length, FUN.VALUE = integer(1))
      # sapply(Transcripts_per_gene_bigGroup, length)
      
      N_genes_Unique_bigGroup = length(classes_split_bigGroup_Unique);
      
      classes_Unique_bigGroup = lapply(seq_len(N_genes_Unique_bigGroup), function(x){
        m = vapply(classes_split_bigGroup_Unique[[x]], function(y){ Transcripts_per_gene_bigGroup[[x]] %in% y }, FUN.VALUE = logical( length( Transcripts_per_gene_bigGroup[[x]] ) ))
        # sapply(classes_split_bigGroup_Unique[[x]], function(y){ Transcripts_per_gene_bigGroup[[x]] %in% y })
        ifelse(m, 1, 0)
      })
      # TO DO: REMOVE CLASSES from transcripts which are not present!
      # AND merge classes which are identical.
      
      ################################################################################################
      # Associate the Median eff length of transcripts:
      ################################################################################################
      df_eff_len_bigGroup = data.frame( Gene_id=rep(seq_len(N_genes_Unique_bigGroup), N_transcripts_per_gene_bigGroup), 
                                        Tr_id = unlist(Transcripts_per_gene_bigGroup))
      df_eff_len_bigGroup$eff_len = eff_len[ match(df_eff_len_bigGroup$Tr_id, names(eff_len))  ]
      
      eff_len_tr_bigGroup = split(df_eff_len_bigGroup$eff_len, df_eff_len_bigGroup$Gene_id)
      
      ################################################################################################
      # Store the info at the end of the Unique genes:
      ################################################################################################
      genes_SELECTED_Unique        = c( genes_SELECTED_Unique, genes_bigGroup)
      Transcripts_per_gene_Unique  = c( Transcripts_per_gene_Unique, Transcripts_per_gene_bigGroup)
      eff_len_tr_Unique            = c( eff_len_tr_Unique, eff_len_tr_bigGroup)
      classes_Unique               = c( classes_Unique, classes_Unique_bigGroup)
      counts_split_per_gene_Unique = c( counts_split_per_gene_Unique, counts_byGene_final)
    }
    # Remove the big groups from the list of groups to be modelled together:
    classes_associated_to_GROUPs = classes_associated_to_GROUPs[-bigGroup]
  }
  
  ######################################################################################################
  # Collect the counts associated to each class
  ######################################################################################################
  # I now need to select the classes and counts associated to each "group" of genes
  # Then "unlist" the classes and counts of each group: turn them into a single class and cont matrix.
  counts_ALL_together_per_GROUP = lapply(classes_associated_to_GROUPs, function(x){
    do.call(rbind, counts_split_per_gene_Together[x])
  })
  #counts_ALL_together_per_GROUP = lapply(classes_associated_to_GROUPs, function(x){
  #res = counts_split_per_gene_Together[[x[1]]]
  # loop only IF length(x) > 1 !!!
  # otherwise it'll create a new record, identical to the previous one
  #if(length(x) > 1){
  #  for(i in 2:length(x)){
  #    res = rbind(res, counts_split_per_gene_Together[[x[i]]])
  #  }
  #  # try to replace the for loop with: do.call(rbind, counts_split_per_gene_Together[x])
  #}
  #res
  #})
  
  classes_ALL_together_per_GROUP = lapply(classes_associated_to_GROUPs, function(x){
    res = list()
    for(i in seq_along(x)){
      res[[i]] = classes_split_per_gene_Together[x[i]]
    }
    res
  })
  
  ######################################################################################################
  # Finally, for each gene, consider the full list of transcripts (to define ... and ...)
  # and make a matrix of 0, 1 indicating, for each class, what transcripts they have.
  Transcripts_per_GROUP = lapply(classes_ALL_together_per_GROUP, function(x){ unique(unlist(x)) })
  N_transcripts_per_GROUP = vapply(Transcripts_per_GROUP, length, FUN.VALUE = integer(1))
  # sapply(Transcripts_per_GROUP, length)
  
  # since classes_ALL_together_per_GROUP[[i]] has a list per gene, 
  # and then each classes_ALL_together_per_GROUP[[i]][[j]][[1]] has a list per equiv class
  # I need to unlist manually classes_ALL_together_per_GROUP
  # otherwise unlist messes up the structure
  N_GROUPS = length(classes_ALL_together_per_GROUP)
  classes_ALL_together_per_GROUP_unlisted = list()
  for(i in seq_len(N_GROUPS) ){
    classes_ALL_together_per_GROUP_unlisted[[i]] = list(classes_ALL_together_per_GROUP[[i]][[1]][[1]][[1]])
    J = length(classes_ALL_together_per_GROUP[[i]])
    for(j in seq_len(J) ){
      K = length(classes_ALL_together_per_GROUP[[i]][[j]][[1]])
      for(k in seq_len(K) ){
        if( sum(c(j == 1, k == 1)) < 2 ){ # the first case (j and k both = 1), is considered above to initialize the list.
          classes_ALL_together_per_GROUP_unlisted[[i]] = c(classes_ALL_together_per_GROUP_unlisted[[i]],
                                                           list(classes_ALL_together_per_GROUP[[i]][[j]][[1]][[k]]))
        }
      }
    }
  }
  
  # From the transcripts, I create the classes
  classes_Together = lapply(seq_along(classes_ALL_together_per_GROUP_unlisted), function(x){
    # m = sapply(classes_ALL_together_per_GROUP_unlisted[[x]], function(y){ Transcripts_per_GROUP[[x]] %in% unlist(y) })
    m = vapply(classes_ALL_together_per_GROUP_unlisted[[x]], function(y){ Transcripts_per_GROUP[[x]] %in% unlist(y) }, FUN.VALUE = logical( length( Transcripts_per_GROUP[[x]] ) ))
    ifelse(m, 1, 0)
  })
  
  ################################################################################################
  # Associate the Median eff length of transcripts:
  ################################################################################################
  n = vapply(Transcripts_per_GROUP, length, FUN.VALUE = integer(1))
  # sapply(Transcripts_per_GROUP, length)
  df_eff_len_Unique = data.frame( Class_id=rep(seq_along(n), n), Tr_id = unlist(Transcripts_per_GROUP))
  df_eff_len_Unique$eff_len = eff_len[ match(df_eff_len_Unique$Tr_id, names(eff_len))  ]
  
  eff_len_tr_Together = split(df_eff_len_Unique$eff_len, df_eff_len_Unique$Class_id)
  
  ######################################################################################################
  # I need an "id" telling me what transcripts belong to what gene in each group.
  ######################################################################################################
  # For each Group, I match (on the truth matrix), the tr_id with the gene_to_transcript
  # 1): list genes in each Group;
  # 2): associate each tr to a gene;
  # 3): count nr of transcripts per gene.
  
  # use the gene_to_transcript matrix
  
  n = vapply(Transcripts_per_GROUP, length, FUN.VALUE = integer(1))
  # sapply(Transcripts_per_GROUP, length)
  df_tmp = data.frame( Class_id=rep(seq_along(n), n), Tr_id = unlist(Transcripts_per_GROUP))
  df_tmp$Gene_id = Gene_id[ match(df_tmp$Tr_id, Tr_id)  ]
  
  genes_per_GROUP = split(df_tmp$Gene_id, df_tmp$Class_id)
  
  genes_per_GROUP_unique_together = lapply(genes_per_GROUP, unique)
  
  for(i in seq_along(Transcripts_per_GROUP)){ # I associate each transcript to the corresponding gene
    names(Transcripts_per_GROUP[[i]]) = genes_per_GROUP[[i]]
  }
  
  ######################################################################################################
  # Filter elements with < 2 transcripts per gene:
  ######################################################################################################
  # Filter Unique genes:
  K = vapply(eff_len_tr_Unique, length, FUN.VALUE = integer(1))
  # sapply(eff_len_tr_Unique, length)
  SEL_tr = K > 1
  genes_SELECTED_Unique        = genes_SELECTED_Unique[SEL_tr]
  Transcripts_per_gene_Unique  = Transcripts_per_gene_Unique[SEL_tr]
  eff_len_tr_Unique            = eff_len_tr_Unique[SEL_tr]
  classes_Unique               = classes_Unique[SEL_tr]
  counts_split_per_gene_Unique = counts_split_per_gene_Unique[SEL_tr]
  
  # store all genes and transcripts that will be analyzed (i.e., in genes with > 2 transcripts):
  K = lapply( seq_len(N_GROUPS),
              function(id){
                # sapply(genes_per_GROUP_unique_together[[id]], function(x) sum(names(Transcripts_per_GROUP[[id]]) == x) )
                vapply(genes_per_GROUP_unique_together[[id]], function(x) sum(names(Transcripts_per_GROUP[[id]]) == x), FUN.VALUE = integer(1))
              }) # compute the number of transcrips belonging to 1 each gene.
  
  all_genes = c(genes_SELECTED_Unique, unlist( lapply(seq_along(K), function(id){ genes_per_GROUP_unique_together[[id]][ K[[id]] > 1] }) ) )
  
  # Filter Together genes (only if all genes in a group have < 2 transcripts)
  SEL_tr = vapply(K, function(x) any(x > 1), FUN.VALUE = logical(1) )
  # sapply(K, function(x) any(x > 1) ) # at least 1 gene with > 1 transcript
  genes_per_GROUP_unique_together = genes_per_GROUP_unique_together[SEL_tr]
  Transcripts_per_GROUP           = Transcripts_per_GROUP[SEL_tr]
  eff_len_tr_Together             = eff_len_tr_Together[SEL_tr]
  classes_Together                = classes_Together[SEL_tr]
  counts_ALL_together_per_GROUP   = counts_ALL_together_per_GROUP[SEL_tr]
  
  ######################################################################################################
  # automatically check the coherence of the matrixed in UNIQUE:
  ######################################################################################################
  if(length(Transcripts_per_gene_Unique) > 0){
    cond = vapply(seq_along(Transcripts_per_gene_Unique), function(i){
      K = length(Transcripts_per_gene_Unique[[i]]); J = ncol(data.frame(classes_Unique[[i]]))
      
      {length(eff_len_tr_Unique[[i]]) == K} & {nrow(data.frame(classes_Unique[[i]])) == K} &  {nrow(counts_split_per_gene_Unique[[i]]) == J} & {ncol(counts_split_per_gene_Unique[[i]]) == N} 
    }, FUN.VALUE = logical(1))
    
    if( mean(cond) != 1){ # Error!
      message("something wrong in Unique classes")
    }
  }
  ######################################################################################################
  # automatically check the coherence of the matrixed in TOGETHER:
  ######################################################################################################
  if(length(Transcripts_per_GROUP) > 0){
    
    cond = vapply(seq_along(Transcripts_per_GROUP), function(i){
      K = nrow(classes_Together[[i]]); J = ncol(classes_Together[[i]])
      
      {length(eff_len_tr_Together[[i]] ) == K} & {length(Transcripts_per_GROUP[[i]]) == K} & {nrow(counts_ALL_together_per_GROUP[[i]]) == J} & {ncol(counts_ALL_together_per_GROUP[[i]]) == N} 
    }, FUN.VALUE = logical(1))
    
    if( mean(cond) != 1){ # Error!
      message("something wrong in Together classes")
    }
  }
  
  # stop clusters:
  if( n_cores > 1.5){ # if n_cores > 1, I use parallel computing tools
    stopCluster(cl) 
  }
  
  ######################################################################################################
  # Make a class out of the results:
  ######################################################################################################
  # Max number of genes and transcripts per GROUP:
  message(paste("Max ", max(c( vapply(eff_len_tr_Unique, length, FUN.VALUE = integer(1)), vapply(eff_len_tr_Together, length, FUN.VALUE = integer(1)) )), " transcripts per group"))
  message(paste("Max ", max( vapply(genes_per_GROUP_unique_together, length, FUN.VALUE = integer(1)) ), " genes per group"))
  
  data = new("BANDITS_data",
             genes       = c(genes_SELECTED_Unique,                 genes_per_GROUP_unique_together), 
             transcripts = c(Transcripts_per_gene_Unique,           Transcripts_per_GROUP),
             effLen      = c(eff_len_tr_Unique,                     eff_len_tr_Together),
             classes     = c(classes_Unique,                        classes_Together),
             counts      = c(counts_split_per_gene_Unique,          counts_ALL_together_per_GROUP), 
             uniqueId    = c( rep(TRUE, length(eff_len_tr_Unique)), rep(FALSE, length(eff_len_tr_Together)) ),
             all_genes   = all_genes )
  
  # return results into a BANDITS_data object:
  return(data)
}
