#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(TRSC))
source("research/simulations/R/scenarios.R")
source("research/simulations/R/metrics.R")
source("research/simulations/R/benchmarks.R")
source("research/simulations/R/run-study.R")
source("research/simulations/R/post-treatment.R")
source("research/simulations/R/reporting.R")
source("research/simulations/R/frozen-study.R")
args<-commandArgs(trailingOnly=TRUE)
has<-function(x)x%in%args
value<-function(x,d=NULL){i<-match(x,args);if(is.na(i))return(d);if(i==length(args)||startsWith(args[i+1],"--"))stop("Missing value for ",x);args[i+1]}
splitv<-function(x)if(is.null(x))NULL else strsplit(x,",",fixed=TRUE)[[1]]
root<-"research/simulations";f<-trsc_sim_validate_frozen_registries(root)
if(has("--preflight")){z<-trsc_sim_stage_a_preflight(root,value("--output-dir",file.path(root,"output","preflight")),value("--ci-tests","not_recorded"),value("--ci-build","not_recorded"),value("--ci-check","not_recorded"));print(z$manifest);quit(status=0)}
out<-value("--output-dir",file.path(root,"results"))
allow<-normalizePath(file.path(root,"results"),mustWork=FALSE);req<-normalizePath(out,mustWork=FALSE)
inside<-identical(req,allow)||startsWith(req,paste0(allow,.Platform$file.sep))
if(!inside)stop("Scientific output must remain below research/simulations/results.")
cell<-splitv(value("--cell"));panel<-splitv(value("--panel"));si<-value("--shard-index");sn<-value("--shard-count")
if(xor(is.null(si),is.null(sn)))stop("Shard index/count must be supplied together.")
cells<-trsc_sim_select_cells(f,cell,panel,if(is.null(si))NULL else as.integer(si),if(is.null(sn))NULL else as.integer(sn))
expected<-trsc_sim_attempt_table(f);expected<-expected[expected$cell_id%in%cells$cell_id,,drop=FALSE]
if(has("--combine")){z<-trsc_sim_combine_frozen(out,f,expected);dir.create(out,recursive=TRUE,showWarnings=FALSE);write.csv(z$attempts,file.path(out,"combined_attempts.csv"),row.names=FALSE);saveRDS(z$records,file.path(out,"combined_records.rds"));cat("Combined",nrow(z$attempts),"declared attempts.\n");quit(status=0)}
if(!nrow(cells))stop("No frozen cells selected.")
cache<-file.path(out,"_fit_cache")
for(i in seq_len(nrow(cells)))for(r in seq_len(as.integer(cells$repetitions[i])))trsc_sim_run_frozen_attempt(cells[i,,drop=FALSE],r,f,out,cache)
cat("Completed selected frozen attempts.\n")
