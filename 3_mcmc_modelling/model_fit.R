library(ape)
library(adephylo)
library(MCMCglmm)
library(reshape2)

path_measure<-function(t){
  opt<-lapply(c("patristic"), function(e) distRoot(t,tips = t$tip.label, method=e))
  names(opt)<- c("pl")
  return(opt)
}

extract_dn_ds_trees<-function(df){
  ds_tree<-paste0(df[(grep("Synonymous tree",df) + 1)],";")
  dn_tree<-paste0(df[(grep("Non-synonymous tree",df) + 1)],";")
  ds_t<-read.tree(text = ds_tree)
  dn_t<-read.tree(text = dn_tree)
  return(list("ds"=ds_t,"dn"=dn_t))
}

tree_rooting<-function(input_t,out_g){
  input_t<-lapply(input_t,root,outgroup=out_g)
  input_t<-lapply(input_t,drop.tip,tip=out_g)
  return(input_t)
}

tree_inverse<-function(tree_input){
  tree_input<-di2multi(tree_input,0.00001)
  return(inverseA(tree_input)$Ainv)
}

find_amphib_outgroup<-function(tree_input,c_data){
  tmp<-c_data[c_data$Taxon %in% tree_input[[1]]$tip.label,]
  tmp<-tmp[tmp$Group == "Amphibia",]
  return(tmp$Taxon[1])
}

no_sites_no_species_extraction<-function(input_data){
  d<-strsplit(input_data[grep("Loaded a multiple sequence",input_data)]," ")[[1]]
  no_sequences<-as.numeric(gsub("\\*\\*","",d[grep("sequences",d)-1]))
  no_sites<-as.numeric(gsub("\\*\\*","",d[grep("codons",d)-1]))*3
  return(c(no_sites,no_sequences))
}

extract_numeric_values <- function(input_string) {
  numeric_values <- as.numeric(unlist(regmatches(input_string, gregexpr("-?\\d+\\.?\\d*", input_string))))
  
  return(numeric_values)
}

meta_data<-function(input_data){
  opt<-no_sites_no_species_extraction(input_data)
  d<-strsplit(input_data[grep("Fitting the baseline model",input_data)+1]," ")[[1]]
  df<-as.data.frame(t(extract_numeric_values(d)))
  df<-cbind(df,t(opt),"Baseline model")
  colnames(df)<-"NULL"
  
  d1<-strsplit(input_data[grep("Improving parameter estimates of the adaptive",input_data)+1]," ")[[1]]
  df1<-as.data.frame(t(extract_numeric_values(d1)))
  df1<-cbind(df1,t(opt),"Adaptive model")
  colnames(df1)<-"NULL"
  
  df<-rbind(df,df1)
  colnames(df)<-c("Log_L","AICc","No_est_para","No_sites","No_species","Model type")
  df[,c(1:5)]<-lapply(df[,c(1:5)],as.numeric)
  return(df)
}

standarize<-function(data){
  return((data-mean(data))/sd(data))
}

pred_standarize<-function(full_data,trimmed_data){
  return((full_data-mean(trimmed_data))/sd(trimmed_data))
}


single_class_model_formatter<-function(model,c_data,input_f){
  tmp<-as.data.frame(model$Sol[,c(1:10)])
  class_vec<-rep(unique(c_data$Class),each=nrow(tmp)*2)
  growth<-rep(c("ben","mal"),each=nrow(tmp))
  ben<-cbind(tmp[,1],tmp[,3],tmp[,5],tmp[,7],tmp[,9])
  mal<-cbind(tmp[,1]+tmp[,2],tmp[,4],tmp[,6],tmp[,8],tmp[,10])
  
  output<-as.data.frame(cbind(class_vec,growth,rbind(ben,mal)))
  colnames(output)<-c("class","growth","intercept","records","bm_slope","dn_pl_slope","ds_pl_slope")
  output[3:7]<-lapply(output[3:7],as.numeric)
  output$gene<-strsplit(input_f,"/")[[1]][1]
  output$no_species<-nrow(c_data)
  return(output)
}

two_class_model_formatter<-function(model,c_data,input_f){
  tmp<-as.data.frame(model$Sol[,c(1:12)])
  class_vec<-rep(c("Aves","Mammalia"),each=nrow(tmp)*2)
  growth<-rep(rep(c("ben","mal"),each=nrow(tmp)),2)
  ben_aves<-cbind(tmp[,1],tmp[,5],tmp[,7],tmp[,9],tmp[,11])
  mal_aves<-cbind(tmp[,1]+tmp[,2],tmp[,6],tmp[,8],tmp[,10],tmp[,12])
  ben_mammals<-cbind(tmp[,1]+tmp[,3],tmp[,5],tmp[,7],tmp[,9],tmp[,11])
  mal_mammals<-cbind(tmp[,1]+tmp[,2]+tmp[,4],tmp[,6],tmp[,8],tmp[,10],tmp[,12])
  
  output<-as.data.frame(cbind(class_vec,growth,rbind(ben_aves,mal_aves,ben_mammals,mal_mammals)))
  colnames(output)<-c("class","growth","intercept","records","bm_slope","dn_pl_slope","ds_pl_slope")
  output[3:7]<-lapply(output[3:7],as.numeric)
  output$gene<-strsplit(input_f,"/")[[1]][1]
  output$no_species<-NA
  output[output$class == "Aves",]$no_species<-nrow(c_data[c_data$Class == "Aves",])
  output[output$class == "Mammalia",]$no_species<-nrow(c_data[c_data$Class == "Mammalia",])
  return(output)
}

parse_pvalues <- function(input) {
  pattern <- "\\* (.*), p-value =\\s*(\\d+\\.\\d+)"
  
  matches <- regmatches(input, regexec(pattern, input))
  
  locations <- sapply(matches, function(x) x[2])
  p_values <- sapply(matches, function(x) as.numeric(x[3]))
  
  df <- data.frame(location = locations, p_value = p_values)
  
  return(df)
}


find_episodic_species<-function(input_file,input_data){
  episodic_evo<-strsplit(input_data[grep("Likelihood ratio test for episodic diversifying",input_data)]," ")[[1]]
  no_episodic_branches<-as.numeric(gsub("([-])|[[:punct:]]","\\1",episodic_evo[18]))
  if (no_episodic_branches > 0){
    start_line<-grep("Likelihood ratio test for episodic diversifying",input_data)+2
    episodic_species<-parse_pvalues(input_data[start_line:length(input_data)])
    write.csv(episodic_species,file=paste0(dirname(input_file),"/episodic_species.csv"),row.names=FALSE)
  }
  if (no_episodic_branches == 0){
    file.create(paste0(dirname(input_file),"/NO_EPISODIC_SPECIES"))
  }
}

format_cancer_data<-function(data,species_removal){
  if ((sum(data$dn.pl) == 0) | (sum(data$ds.pl) == 0)){
    return(list(c(),c()))
  }
  if ((sum(data$dn.pl) != 0) & (sum(data$ds.pl) != 0)){
    if (nrow(data[data$ds.pl == 0,]) > 0){
      data[data$ds.pl == 0,]$ds.pl<-min(data[data$ds.pl > 0,]$ds.pl)
    }
    if (nrow(data[data$dn.pl == 0,] > 0)){
      data[data$dn.pl == 0,]$dn.pl<-min(data[data$dn.pl > 0,]$dn.pl)
    }
    data$log_ds_pl<-log(data$ds.pl)
    data$log_dn_pl<-log(data$dn.pl)
    
    for (i in unique(data$Class)){
      if (nrow(data[(data$Class == i) & !(data$Species %in% species_removal),]) < 5){
        data<-data[data$Class != i,]
      }
    }
    for (i in unique(data$Class)){
      if ((sd(data[(data$Class == i) & !(data$Species %in% species_removal),]$log_ds_pl) == 0) | (sd(data[(data$Class == i) & !(data$Species %in% species_removal),]$log_dn_pl) == 0)){
        data<-data[data$Class != i,]
      }
    }
    if (nrow(data) == 0){
      return(list(c(),c()))
    }
    if (nrow(data) != 0){
      full_data<-data
      data<-data[!(data$Species %in% species_removal),]
      
      data$std_log_rec<-NA
      data$std_log_bm<-NA
      data$std_log_ds_pl<-NA
      data$std_log_dn_pl<-NA
      
      full_data$std_log_rec<-NA
      full_data$std_log_bm<-NA
      full_data$std_log_ds_pl<-NA
      full_data$std_log_dn_pl<-NA
      
      for (i in unique(data$Class)){
        data[data$Class == i,]$std_log_rec<-standarize(data[data$Class == i,]$log_rec)
        data[data$Class == i,]$std_log_bm<-standarize(data[data$Class == i,]$log_bm)
        data[data$Class == i,]$std_log_ds_pl<-standarize(data[data$Class == i,]$log_ds_pl)
        data[data$Class == i,]$std_log_dn_pl<-standarize(data[data$Class == i,]$log_dn_pl)
        
        full_data[full_data$Class == i,]$std_log_rec<-pred_standarize(full_data[full_data$Class == i,]$log_rec,data[data$Class == i,]$log_rec)
        full_data[full_data$Class == i,]$std_log_bm<-pred_standarize(full_data[full_data$Class == i,]$log_bm,data[data$Class == i,]$log_bm)
        full_data[full_data$Class == i,]$std_log_ds_pl<-pred_standarize(full_data[full_data$Class == i,]$log_ds_pl,data[data$Class == i,]$log_ds_pl)
        full_data[full_data$Class == i,]$std_log_dn_pl<-pred_standarize(full_data[full_data$Class == i,]$log_dn_pl,data[data$Class == i,]$log_dn_pl)
      }
      data$mammal = ifelse(data$Class == "Mammalia", 1, 0)
      full_data$mammal = ifelse(full_data$Class == "Mammalia", 1, 0)
      return(list(data,full_data))
    } 
  }
}


phylo_prediction<-function(model,data,file_location){
  pred<-as.data.frame(predict.MCMCglmm(model,data,interval="none",type="response",marginal=NULL))
  
  pred_df<-cbind(pred[1:nrow(data),1],pred[(nrow(data)+1):(2*nrow(data)),1])
  colnames(pred_df)<-c("ben_fit","mal_fit")
  
  output<-cbind(data,pred_df)
  
  output$ben_res<-log(abs(output$ben - output$ben_fit))
  output$std_ben_res<-(output$ben_res - mean(output$ben_res))/sd(output$ben_res)
  output$abs_std_ben_res<-abs(output$std_ben_res)
  
  output$mal_res<-log(abs(output$mal - output$mal_fit))
  output$std_mal_res<-(output$mal_res - mean(output$mal_res))/sd(output$mal_res)
  output$abs_std_mal_res<-abs(output$std_mal_res)
  
  output$class<-NA
  if (length(unique(data$Class)) == 1){
    output$hyphy_intercept<-unique(data$Class)
  }
  if (length(unique(data$Class)) == 2){
    output$hyphy_intercept<-"both"
  }
  
  write.csv(output[,c("Species","ben_res","mal_res","std_ben_res","std_mal_res","abs_std_ben_res","abs_std_mal_res","hyphy_intercept")],paste0(dirname(file_location),"/mcmcglmm_output/prediction_results.csv"),row.names=FALSE)
}


eff_sample_covariance_formatter<-function(model,file_location){
  mod_sum<-summary(model)
  eff_sample<-c(mod_sum$Gcovariances[c(1,2,4),4],mod_sum$Rcovariances[c(1,2,4),4])
  type<-rep(c("ben_var","ben_mal_cov","mal_var"),2)
  var_structure<-rep(c("G","R"),each=3)
  output<-as.data.frame(cbind(eff_sample,type,var_structure))
  row.names(output)<-NULL
  write.csv(output,paste0(dirname(file_location),"/mcmcglmm_output/covariance_effective_sample.csv"),row.names=FALSE)
}

eff_sample_covariate_formatter<-function(model,file_location){
  mod_sum<-summary(model)
  eff_sample<-mod_sum$solutions[,4]
  idx<-row.names(mod_sum$solutions)
  output<-as.data.frame(cbind(idx,eff_sample))
  row.names(output)<-NULL
  write.csv(output,paste0(dirname(file_location),"/mcmcglmm_output/covariate_effective_sample.csv"),row.names=FALSE)
}

neg_branch_check<-function(tt){
  t1<-any(tt[[1]]$edge.length < 0)
  t2<-any(tt[[2]]$edge.length < 0)
  if ((t1 == FALSE) & (t2 == FALSE)){
    return(FALSE)
  }
  else{
    return(TRUE)
  }
  return(any(tt$edge.length < 0))
}



main<-function(input_file){
  input_data<-readLines(input_file)
  
  dir.create(file.path(dirname(input_file), "mcmcglmm_output"), showWarnings = FALSE)
  
  model_meta_data<-meta_data(input_data)
  write.csv(model_meta_data,paste0(dirname(input_file),"/absrel_meta_data.csv"),row.names=FALSE)
  
  hyphy_trees<-extract_dn_ds_trees(input_data)
  write.tree(hyphy_trees,paste0(dirname(input_file),"/full_NOT_ROOTED_dn_ds_trees.trees"))
  species_data<-read.csv("./Full_taxa_list.csv")
  outgroup<-find_amphib_outgroup(hyphy_trees,species_data)
  
  
  if (is.na(outgroup) == TRUE){
    file.create(paste0(dirname(input_file),"/NO_ROOT_IN_ALIGNMENT"))
  }
  
  if (is.na(outgroup) == FALSE){
    cancer_data_full<-read.csv("./species_data.csv")
    
    cancer_data_full$ben<-cancer_data_full$neo-cancer_data_full$mal
    
    
    cancer_data_full<-cancer_data_full[cancer_data_full$Class %in% c("Aves","Mammalia"),]
    
    hyphy_trees<-tree_rooting(hyphy_trees,outgroup)
    species_removal<-hyphy_trees$ds$tip.label[!(hyphy_trees$ds$tip.label %in% cancer_data_full$Species)]
    hyphy_trees1<-lapply(hyphy_trees,drop.tip,species_removal)
    write.tree(hyphy_trees,paste0(dirname(input_file),"/cancer_matched_dn_ds_trees.trees"))
    
    if (neg_branch_check(hyphy_trees) == TRUE){
      file.create(paste0(dirname(input_file),"/mcmcglmm_output/NEGATIVE_BRANCH"))
    }
    
    if (neg_branch_check(hyphy_trees) == FALSE){
      pl_output<-as.data.frame(sapply(hyphy_trees,path_measure))
      pl_output$species<-row.names(pl_output)
      
      cancer_data_full<-merge(cancer_data_full,pl_output,by.x="Species",by.y="species")
      species_to_remove<-c("Didelphis_marsupialis","Acrobates_pygmaeus","Petaurus_breviceps","Dendrolagus_goodfellowi","Phascolarctos_cinereus")
      
      
      opt<-format_cancer_data(cancer_data_full,species_to_remove)
      cancer_data<-opt[[1]]
      cancer_data_full<-opt[[2]]
      
      if (length(cancer_data) == 0){
        file.create(paste0(dirname(input_file),"/mcmcglmm_output/NO_PL_VARIATION"))
      }
      if (length(cancer_data) > 0){
        tree<-read.nexus("./terrestrial_vertebrate_tree.nexus.trees")
        tree<-keep.tip(tree,cancer_data_full$Species)
        
        treeAinv<-tree_inverse(tree)
        
        number_iterations<-1e7
        sample_interval<-1e3
        burnin_iterations<-9e6
        
        
        no_classes<-length(unique(cancer_data$Class))
        write.csv(cancer_data,paste0(dirname(input_file),"/mcmcglmm_output/model_cancer_data.csv"),row.names=FALSE)
        
        if (no_classes == 1){
          pr1<- list(R=list(V=diag(2), nu=2, fix=2),B=list(mu=(rep(0,10)),V=diag(10)*1e8),
                     G=list(G1=list(V=diag(2), nu=2, alpha.mu=cbind(0,0), alpha.V=diag(2)*25^2)))
          
          std_mod1<-MCMCglmm(cbind(ben,mal) ~ trait+trait:std_log_rec + trait:std_log_bm + trait:std_log_dn_pl + trait:std_log_ds_pl, family=cbind("poisson","poisson"), pr=TRUE,prior = pr1, data=cancer_data, pl=TRUE,
                             rcov = ~ us(trait):units, random = ~us(trait):Species,nitt=number_iterations,thin=sample_interval,burnin=burnin_iterations, verbose=FALSE, DIC=TRUE, ginverse=list(Species=treeAinv))
          output<-single_class_model_formatter(std_mod1,cancer_data,input_file)
          write.csv(output,paste0(dirname(input_file),"/mcmcglmm_output/std_single_intercept_single_slope.csv"),row.names=FALSE)
          saveRDS(std_mod1,paste0(dirname(input_file),"/mcmcglmm_output/std_fitted_model.RDS"))
          
        }
        
        if (no_classes == 2){
          pr1<- list(R=list(V=diag(2), nu=2, fix=2),B=list(mu=(rep(0,12)),V=diag(12)*1e8),
                     G=list(G1=list(V=diag(2), nu=2, alpha.mu=cbind(0,0), alpha.V=diag(2)*25^2)))
          
          std_mod1<-MCMCglmm(cbind(ben,mal) ~ trait+trait:mammal+trait:std_log_rec + trait:std_log_bm + trait:std_log_dn_pl + trait:std_log_ds_pl, family=cbind("poisson","poisson"), pr=TRUE,prior = pr1, data=cancer_data, pl=TRUE,
                             rcov = ~ us(trait):units, random = ~us(trait):Species,nitt=number_iterations,thin=sample_interval,burnin=burnin_iterations, verbose=FALSE, DIC=TRUE, ginverse=list(Species=treeAinv))
          output<-two_class_model_formatter(std_mod1,cancer_data,input_file)
          write.csv(output,paste0(dirname(input_file),"/mcmcglmm_output/std_diff_intercept_single_slope.csv"),row.names=FALSE)
          saveRDS(std_mod1,paste0(dirname(input_file),"/mcmcglmm_output/std_fitted_model.RDS"))
          
        }
        
        find_episodic_species(input_file,input_data)
        eff_sample_covariance_formatter(std_mod1,input_file)
        eff_sample_covariate_formatter(std_mod1,input_file)
        
        phylo_prediction(std_mod1,cancer_data_full[cancer_data_full$Class %in% cancer_data$Class,],input_file)
      }
    }
  }
}

setwd("~/Documents/GitHub/genome_wide_cancer_dynamics/3_mcmc_modelling/")

main("./ARL14EP/hyphy_aBSREL_ARL14EP/output_ARL14EP_absrel.txt")
