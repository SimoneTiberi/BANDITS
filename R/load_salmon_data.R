load_salmon_data = function(sep,
                            path_to_eq_classes,
                            n_cores, cl = NULL,
                            transcripts_to_keep){
  if( !all(file.exists(path_to_eq_classes)) ){ # if at least 1 file not found
    message("'path_to_eq_classes' files not found")
    return(NULL)
  }
  
  N = length(path_to_eq_classes)
  
  # MAKE SURE THAT "_" is not present in the list of transcript ids:
  # otherwise crease a longer list of "_".
  
  if( n_cores > 1.5){ # if n_cores > 1, I use parallel computing tools
    if(is.null(transcripts_to_keep)){
      x = parLapply(cl = cl, X = path_to_eq_classes, fun = salmon_read_eq_classes, sep = sep)
    }else{ # if transcripts_to_keep have been specified, I filter them out from the equivalence classes:
      x = parLapply(cl = cl, X = path_to_eq_classes, fun = salmon_read_eq_classes_filteringTranscripts,
                    transcripts_to_keep = transcripts_to_keep, sep = sep)
    }
  }else{
    if(is.null(transcripts_to_keep)){
      x = lapply(X = path_to_eq_classes, FUN = salmon_read_eq_classes, sep = sep)
    }else{ # if transcripts_to_keep have been specified, I filter them out from the equivalence classes:
      x = lapply(X = path_to_eq_classes, FUN = salmon_read_eq_classes_filteringTranscripts, 
                 transcripts_to_keep = transcripts_to_keep, sep = sep)
    }
  }
  
  return(x)
}


salmon_read_eq_classes = function(fn, sep){
  fr = fread(fn, sep = " ", quote = "", header = FALSE)
  ids = fr$V1[3:(as.integer(fr$V1[1])+2)] # vector with all transcript ids
  ecs = fr$V1[(as.integer(fr$V1[1])+3):nrow(fr)]
  ecs.s = strsplit(ecs,"\t",fixed=TRUE)
  
  cnt = as.integer(vapply(ecs.s, last, FUN.VALUE = character(1)))
  # as.integer(sapply(ecs.s, last)) # counts for each equiv class
  trans = lapply(ecs.s, function(u) 1+as.integer(u[2:(length(u)-1)]))
  # sapply(ecs.s, function(u) 1+as.integer(u[2:(length(u)-1)]))
  
  class_ids = vapply(trans, function(u) paste(sort(ids[unique(u)]),collapse=sep), FUN.VALUE = "id") # I create the id for the class, made of all transcripts of the class separated by _
  
  # check if there are any duplicated classes:  
  if( sum(duplicated(class_ids)) > 0.5 ){ # use the other method, with transcripts_to_keep = all transcripts
    return( read_eq_classes_filteringTranscripts(fn = fn, transcripts_to_keep = unique(ids),  sep = sep) )
  }
  
  list(counts=cnt, class_ids=class_ids)
}


salmon_read_eq_classes_filteringTranscripts = function(fn, transcripts_to_keep, sep) {
  fr = fread(fn, sep = " ", quote = "", header = FALSE)
  ids = fr$V1[3:(as.integer(fr$V1[1])+2)] # vector with all transcript ids
  ecs = fr$V1[(as.integer(fr$V1[1])+3):nrow(fr)]
  ecs.s = strsplit(ecs,"\t",fixed=TRUE)
  
  cnt = as.integer(vapply(ecs.s, last, FUN.VALUE = character(1)))
  # as.integer(sapply(ecs.s, last)) # counts for each equiv class
  trans = lapply(ecs.s, function(u) 1+as.integer(u[2:(length(u)-1)]))
  # sapply(ecs.s, function(u) 1+as.integer(u[2:(length(u)-1)]))
  trans = lapply(trans, as.integer)
  trans = lapply(trans, unique)
  
  # sort maybe not needed, double-check at the end!!!
  SEL_transcripts  = which(ids %in% transcripts_to_keep)
  
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
  
  class_ids_sel = vapply(trans_sel[SEL_classes], function(u) paste(sort(ids[u]),collapse=sep), FUN.VALUE = "id") # I create the id for the class, made of all transcripts of the class separated by _
  
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
