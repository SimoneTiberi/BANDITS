load_kallisto_data = function(sep,
                              kallisto_equiv_classes,
                              kallisto_equiv_counts,
                              kallisto_counts,
                              n_cores, cl = NULL,
                              transcripts_to_keep){
  # check that kallisto_counts is a matrix or data.frame object
  if( !is.data.frame(kallisto_counts) & !is.matrix(kallisto_counts) ){
    message("'kallisto_counts' must be a matrix or data.frame")
    return(NULL)
  }
  
  if( !all(dim(kallisto_counts) > 0.5) ){
    message("'kallisto_counts' must have at least 1 row (transcripts) and 1 column (samples)")
    return(NULL)
  }
  
  if( !all(file.exists(kallisto_equiv_classes)) ){ # if at least 1 file not found
    message("'kallisto_equiv_classes' files not found")
    return(NULL)
  }
  
  if( !all(file.exists(kallisto_equiv_counts)) ){ # if at least 1 file not found
    message("'kallisto_equiv_counts' files not found")
    return(NULL)
  }
  
  N = length(kallisto_equiv_classes)
  
  kallisto_transcript_names = rownames(kallisto_counts)
  
  if( n_cores > 1.5){ # if n_cores > 1, I use parallel computing tools
    if(is.null(transcripts_to_keep)){
      x = parLapply(cl = cl, X = seq_len(N), fun = kallisto_read_eq_classes, sep = sep,
                    kallisto_equiv_classes = kallisto_equiv_classes, kallisto_equiv_counts = kallisto_equiv_counts,
                    kallisto_transcript_names = kallisto_transcript_names)
    }else{ # if transcripts_to_keep have been specified, I filter them out from the equivalence classes:
      x = parLapply(cl = cl, X = seq_len(N), fun = kallisto_read_eq_classes_filteringTranscripts, transcripts_to_keep = transcripts_to_keep, sep = sep,
                    kallisto_equiv_classes = kallisto_equiv_classes, kallisto_equiv_counts = kallisto_equiv_counts,
                    kallisto_transcript_names = kallisto_transcript_names)
    }
    
  }else{
    if(is.null(transcripts_to_keep)){
      x = lapply(X = seq_len(N), FUN = kallisto_read_eq_classes, sep = sep,
                 kallisto_equiv_classes = kallisto_equiv_classes, kallisto_equiv_counts = kallisto_equiv_counts,
                 kallisto_transcript_names = kallisto_transcript_names)
    }else{ # if transcripts_to_keep have been specified, I filter them out from the equivalence classes:
      x = lapply(X = seq_len(N), FUN = kallisto_read_eq_classes_filteringTranscripts, transcripts_to_keep = transcripts_to_keep, sep = sep,
                 kallisto_equiv_classes = kallisto_equiv_classes, kallisto_equiv_counts = kallisto_equiv_counts,
                 kallisto_transcript_names = kallisto_transcript_names)
    }
  }
  
  return(x)
}


kallisto_read_eq_classes = function(X, sep, kallisto_equiv_classes, kallisto_equiv_counts, kallisto_transcript_names){
  ec = fread(kallisto_equiv_classes[[X]], sep = "\t", quote = "", header = FALSE)
  ecs = ec$V2
  ecs.s = strsplit(ecs,",",fixed=TRUE)
  
  tsv = fread(kallisto_equiv_counts[[X]], sep = "\t", quote = "", header = FALSE)
  cnt = as.integer(tsv$V2)
  
  trans = lapply(ecs.s, function(u) 1+as.integer(u))
  # 1) I turn ecs.s from character to numeric and 2) I add 1 because indexes in ecs.s start from 0.
  
  class_ids = vapply(trans, function(u) paste(sort(kallisto_transcript_names[unique(u)]),collapse=sep), FUN.VALUE = "id")
  # I create the id for the class, made of all transcripts of the class separated by '_' (1 or more '_').
  
  # check if there are any duplicated classes:  
  if( sum(duplicated(class_ids)) > 0.5 ){ # use the other method, with transcripts_to_keep = all transcripts
    return( read_eq_classes_filteringTranscripts(fn = fn, transcripts_to_keep = unique(kallisto_transcript_names),  sep = sep) )
  }
  
  list(counts=cnt, class_ids=class_ids)
}


kallisto_read_eq_classes_filteringTranscripts = function(X, transcripts_to_keep, sep, 
                                                         kallisto_equiv_classes, kallisto_equiv_counts,
                                                         kallisto_transcript_names) {
  ec = fread(kallisto_equiv_classes[[X]], sep = "\t", quote = "", header = FALSE)
  ecs = ec$V2
  ecs.s = strsplit(ecs,",",fixed=TRUE)
  
  tsv = fread(kallisto_equiv_counts[[X]], sep = "\t", quote = "", header = FALSE)
  cnt = as.integer(tsv$V2)
  
  trans = lapply(ecs.s, function(u) 1+as.integer(u))
  trans = lapply(trans, as.integer)
  trans = lapply(trans, unique)
  
  # sort maybe not needed, double-check at the end!!!
  SEL_transcripts  = which(kallisto_transcript_names %in% transcripts_to_keep)
  
  #############  #############  #############  #############  #############  #############
  n = vapply(trans, length, FUN.VALUE = integer(1))
  # sapply(trans, length)
  df = data.frame(Class_id=rep(seq_along(n), n), Tr_id = unlist(trans))
  
  #df = df[df$Tr_id %in% SEL_transcripts, ]
  #trans_sel = split(df$Tr_id, df$Class_id)
  
  # ALTERNATIVE trans_sel computation (much faster):
  df$Tr_id_SEL = ifelse(df$Tr_id %in% SEL_transcripts, df$Tr_id, NA)
  trans_sel = split(df$Tr_id_SEL, df$Class_id)
  
  trans_sel = lapply(trans_sel, function(x) x[!is.na(x)])
  # sapply(trans_sel, function(x) x[!is.na(x)])
  
  #############  #############  #############  #############  #############  #############
  # Now I filter out classes with NO transcripts (where all transcripts were filtered out):
  SEL_classes = {vapply(trans_sel, length, FUN.VALUE = integer(1)) > 0.5} # if length == 0, I have 0 transcripts in a class
  
  class_ids_sel = vapply(trans_sel[SEL_classes], function(u) paste(sort(kallisto_transcript_names[u]),collapse=sep), FUN.VALUE = "id") # I create the id for the class, made of all transcripts of the class separated by _
  
  cnt_sel = cnt[SEL_classes]
  
  # do UNIQUE of classes:
  # I CAN TURN class_ids_sel A NUMERIC ID (factor):
  class_ids_sel_num = as.numeric(factor( class_ids_sel ))
  DUPS = duplicated(class_ids_sel_num)
  class_ids_sel_unique_num = class_ids_sel_num[DUPS == FALSE]
  
  # SUM up counts of repeated classes:
  cnt_sel_unique = cnt_sel[DUPS == FALSE] # I filter the counts of the unique classes
  cnt_sel_DUPS   = cnt_sel[DUPS] # I filter the counts of the unique classes
  
  cnt_sel_unique = vapply(seq_along(class_ids_sel_unique_num), function(i){
    match_dups = class_ids_sel_num[DUPS] == class_ids_sel_unique_num[i] # I look for the matching between the unique classes and duplicatd ones (eliminated)
    if(sum(match_dups) > 0){
      cnt_sel_unique[i] = cnt_sel_unique[i] + sum(cnt_sel_DUPS[match_dups]) # I add the counts of the corresponding duplicated classes
    }
    cnt_sel_unique[i]
  }, FUN.VALUE = integer(1))
  #  2.5 times faster than previous cnt_sel_unique computation!
  
  list(counts=cnt_sel_unique, class_ids=class_ids_sel[DUPS == FALSE])
}