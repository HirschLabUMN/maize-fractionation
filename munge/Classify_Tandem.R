#   Script to classify a tandem duplicate cluster as likely not real tandem
#   duplication based on pairwise sequence similarity. Will also flag some
#   clusters for manual inspection. These depend on tables of sequence summary
#   statistics calculated with 'compute' from the 'analysis' package, written
#   by K. Thornton.

#   Read the data tables
b73_ancestral_table <- read.table("/Users/tomkono/Dropbox/Projects/Fractionation/Data/Tandem_Divergence/B73_Outgroup_Summary.txt", header=TRUE, skip=3)
ph207_ancestral_table <- read.table("/Users/tomkono/Dropbox/Projects/Fractionation/Data/Tandem_Divergence/PH207_Outgroup_Summary.txt", header=TRUE, skip=3)

b73_ingroup_table <- read.table("/Users/tomkono/Dropbox/Projects/Fractionation/Data/Tandem_Divergence/B73_MzOnly_Summary.txt", header=TRUE, skip=3)
ph207_ingroup_table <- read.table("/Users/tomkono/Dropbox/Projects/Fractionation/Data/Tandem_Divergence/PH207_MzOnly_Summary.txt", header=TRUE, skip=3)

b73_clusters <- read.table("/Users/tomkono/Dropbox/Projects/Fractionation/Data/Tandem_Divergence/B73_Cluster_Numbers.txt", header=FALSE)
ph207_clusters <- read.table("/Users/tomkono/Dropbox/Projects/Fractionation/Data/Tandem_Divergence/PH207_Cluster_Numbers.txt", header=FALSE)

#   Then, apply the filters
pi_thresh <- 0.25
gapped <- 0.5

b73_fail_pi <- b73_ancestral_table$nsam == 2 &
    !is.na(b73_ancestral_table$ThetaPi) &
    b73_ancestral_table$ThetaPi > pi_thresh
b73_fail_gap <- (b73_ingroup_table$nsites_ug/b73_ingroup_table$nsites) < gapped
#   Get the names of those that fail for similarity and gapping
b73_fail_pi_names <- as.character(b73_ancestral_table$locus[b73_fail_pi])
b73_fail_gap_names <- as.character(b73_ingroup_table$locus[b73_fail_gap])
#   We need to strip the 'Anc_' from the front of the names of those that fail
#   on the similarity filter
b73_fail_pi_names <- gsub("Anc_", "", b73_fail_pi_names)

#   Get the names of those that are NA in similarity - these require some manual
#   inspection
b73_manual_names <- b73_ancestral_table$locus[is.na(b73_ancestral_table$ThetaPi)]
b73_manual_names <- gsub("Anc_", "", as.character(b73_manual_names))

#   What is left after these two are probably tandem duplicates
b73_real <- !(gsub("Anc_", "", as.character(b73_ancestral_table$locus)) %in% c(b73_manual_names, b73_fail_pi_names, b73_fail_gap_names))
b73_real <- gsub("Anc_", "", as.character(b73_ancestral_table$locus[b73_real]))

#   Make a table that holds the status of each filter
b73_pi_flt <- sapply(as.character(b73_ingroup_table$locus), function(x) {
    if(x %in% b73_fail_pi_names & b73_ingroup_table$nsam[b73_ingroup_table$locus ==x] == 2) {
        return("Fail")
    } else {
        return("Pass")
    }
    })
b73_cov_flt <- sapply(as.character(b73_ingroup_table$locus), function(x) {
    if(x %in% b73_fail_gap_names & b73_ingroup_table$nsam[b73_ingroup_table$locus == x] == 2) {
        return("Fail")
    } else {
        return("Pass")
    }
    })
b73_manual <- sapply(as.character(b73_ingroup_table$locus), function(x) {
    if(x %in% b73_manual_names | ((x %in% b73_fail_pi_names | x %in% b73_fail_gap_names) & b73_ingroup_table$nsam[b73_ingroup_table$locus == x] != 2)) {
        return("True")
    } else {
        return("False")
    }
    })
b73_syn <- sapply(as.character(b73_ingroup_table$locus), function(x) {
    anc_name <- paste("Anc_", x, sep="")
    if(anc_name %in% b73_ancestral_table$locus) {
        return("True")
    } else {
        return("False")
    }
    })
b73_looksreal <- sapply(as.character(b73_ingroup_table$locus), function(x) {
    if(x %in% b73_real) {
        return("True")
    } else {
        return("False")
    }
    })

b73_status <- data.frame(
    ClusterName=gsub(".fasta", "", as.character(b73_ingroup_table$locus)),
    Similarity=b73_pi_flt,
    Coverage=b73_cov_flt,
    Manual_Check=b73_manual,
    Has_Syn_Orthologue=b73_syn,
    Looks_Real=b73_looksreal,
    Genes=b73_clusters$V2
    )

#   And write it out
write.table(b73_status, file="B73_Cluster_Status.txt", sep="\t", quote=FALSE, row.names=FALSE)


#   Repeat for PH207
ph207_fail_pi <- ph207_ancestral_table$nsam == 2 &
    !is.na(ph207_ancestral_table$ThetaPi) &
    ph207_ancestral_table$ThetaPi > pi_thresh
ph207_fail_gap <- (ph207_ingroup_table$nsites_ug/ph207_ingroup_table$nsites) < gapped
#   Get the names of those that fail for similarity and gapping
ph207_fail_pi_names <- as.character(ph207_ancestral_table$locus[ph207_fail_pi])
ph207_fail_gap_names <- as.character(ph207_ingroup_table$locus[ph207_fail_gap])
#   We need to strip the 'Anc_' from the front of the names of those that fail
#   on the similarity filter
ph207_fail_pi_names <- gsub("Anc_", "", ph207_fail_pi_names)

#   Get the names of those that are NA in similarity - these require some manual
#   inspection
ph207_manual_names <- ph207_ancestral_table$locus[is.na(ph207_ancestral_table$ThetaPi)]
ph207_manual_names <- gsub("Anc_", "", as.character(ph207_manual_names))

#   What is left after these two are probably tandem duplicates
ph207_real <- !(gsub("Anc_", "", as.character(ph207_ancestral_table$locus)) %in% c(ph207_manual_names, ph207_fail_pi_names, ph207_fail_gap_names))
ph207_real <- gsub("Anc_", "", as.character(ph207_ancestral_table$locus[ph207_real]))

#   Make a table that holds the status of each filter
ph207_pi_flt <- sapply(as.character(ph207_ingroup_table$locus), function(x) {
    if(x %in% ph207_fail_pi_names & ph207_ingroup_table$nsam[ph207_ingroup_table$locus ==x] == 2) {
        return("Fail")
    } else {
        return("Pass")
    }
    })
ph207_cov_flt <- sapply(as.character(ph207_ingroup_table$locus), function(x) {
    if(x %in% ph207_fail_gap_names & ph207_ingroup_table$nsam[ph207_ingroup_table$locus == x] == 2) {
        return("Fail")
    } else {
        return("Pass")
    }
    })
ph207_manual <- sapply(as.character(ph207_ingroup_table$locus), function(x) {
    if(x %in% ph207_manual_names | ((x %in% ph207_fail_pi_names | x %in% ph207_fail_gap_names) & ph207_ingroup_table$nsam[ph207_ingroup_table$locus == x] != 2)) {
        return("True")
    } else {
        return("False")
    }
    })
ph207_syn <- sapply(as.character(ph207_ingroup_table$locus), function(x) {
    anc_name <- paste("Anc_", x, sep="")
    if(anc_name %in% ph207_ancestral_table$locus) {
        return("True")
    } else {
        return("False")
    }
    })
ph207_looksreal <- sapply(as.character(ph207_ingroup_table$locus), function(x) {
    if(x %in% ph207_real) {
        return("True")
    } else {
        return("False")
    }
    })

ph207_status <- data.frame(
    ClusterName=gsub(".fasta", "", as.character(ph207_ingroup_table$locus)),
    Similarity=ph207_pi_flt,
    Coverage=ph207_cov_flt,
    Manual_Check=ph207_manual,
    Has_Syn_Orthologue=ph207_syn,
    Looks_Real=ph207_looksreal,
    Genes=ph207_clusters$V2
    )

#   And write it out
write.table(ph207_status, file="PH207_Cluster_Status.txt", sep="\t", quote=FALSE, row.names=FALSE)
