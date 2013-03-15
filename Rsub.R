MakeShell <- function(pbs.file, job.name, cmd, mem, clean) {
  write("#$ -cwd", pbs.file, append=FALSE)
  write(paste("#$ -l vf=", mem, sep=""), pbs.file, append=TRUE)
  write(paste("#$ -N", job.name), pbs.file, append=TRUE)
  
  write("#$ -m eas", pbs.file, append=TRUE)
  write("#$ -M enda.byrne@uq.edu.au", pbs.file, append=TRUE)
  write(cmd, pbs.file, append=TRUE)  
}

gmdr <- 'java -Xmx5G -jar /clusterdata/gc5k/bin/gmdr.jar'
gcta <- '/clusterdata/apps/gcta/gcta64'
gcta_test <- '/clusterdata/gc5k/bin/gcta64_test'
gcta_beta <- '/clusterdata/jyang/bin/gcta64_test'
plink <- '/clusterdata/apps/plink/plink-1.07-x86_64/plink'
polygenic <- 'java -jar /clusterdata/gc5k/bin/polygenic.jar'
HE <- '/clusterdata/gc5k/bin/HE.jar'
###################################################################

args <- commandArgs(TRUE)
CMD <- c()
nm <- c()
mem <- "10G"
runit <- TRUE

for (i in 1:length(args)) {
  flag <- TRUE
  breakFlag <- FALSE  
  cmdbit <- args[i]
  if (substr(args[i], nchar(args[i]), nchar(args[i])) == ",") {
    breakFlag <- TRUE
    cmdbit <- substr(args[i], 1, nchar(args[i])-1)
  }

  if (substr(cmdbit, 1, 1) == "@") {
    nm <- substr(cmdbit, 2, nchar(cmdbit))
    flag <- FALSE
  } else if (substr(cmdbit, 1, 1) == "%") {
    mem <- substr(cmdbit,2, nchar(cmdbit))
    flag <- FALSE
  } else if (substr(cmdbit, 1, 2) == "TT" && nchar(cmdbit)==2) {
    runit <- FALSE
    flag <- FALSE
  } else if (cmdbit == 'gmdr') {
    cmd <- gmdr
  } else if (cmdbit == 'HE' ) {
    cmd <- HE
  } else if (cmdbit == 'gcta_beta' ) {
    cmd <- gcta_beta
  } else if (cmdbit == 'gcta_test' ) {
    cmd <- gcta_test
  } else if (cmdbit == 'gcta') {
    cmd <- gcta
  } else if (cmdbit == 'plink') {
    cmd <- plink
  } else if (cmdbit == 'polygenic') {
    cmd <- polygenic
  } else {
    cmd <- cmdbit
  }

  if(flag) {
    CMD <- paste(CMD, cmd)
    if(breakFlag) {
      CMD <- paste(CMD, "\n")
    }
  }

}

print(CMD)

if (length(nm) > 0) {
  pbs.file <- paste(nm, ".sh", sep="")
  job.name <- paste(nm, sep="")
} else {
  pbs.file <- paste(args[1], args[length(args)], ".sh", sep="-")
  job.name <- paste(args[1], args[length(args)], sep="-")
}

MakeShell(pbs.file, job.name, CMD, mem)
sys_qsub <- paste("qsub", pbs.file)
if(runit) {
  system(sys_qsub)
}
