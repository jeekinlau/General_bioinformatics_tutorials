####################################################################################################################################
#
#
#   Polyploid computational support meeting 12/05/2023
#   
#   Jeekin Lau
#
# Install the latest update 
# Most of the updates we are going to talk about are in the newest developmental build of mappoly
#
# 1) introduce new functions for filtering individuals filter_individuals()
# 2) show how to use new function for small inversions
#      a) use MDS to find inversions, then fix them for mapping
# 3) show how to use new pair of functions framework_map() and update_framework_map()
# 4) show a new function for showing the changes the hmm chain makes on the genotypic matrix plot_progeny_dosage_change()
#
#
#####################################################################################################################################



# if you dont have the latest version
devtools::install_github("mmollina/mappoly")
# Load the library
library(mappoly)


# Read and filter the data without conforming markers reason: it does not work well if we filter first
dat = read_geno_csv(file="https://raw.githubusercontent.com/jeekinlau/mapping_example/main/data/BExMG_subset_with_contaminants_2075mrks.csv",
                    ploidy = 4, filter.non.conforming = F)
print(dat, detailed = T) 


##### new QC step for finding individuals that dont fit the expected assumptions in the mapping populations
# show both ways to filter
?filter_individuals
dat <- filter_individuals(dat, type = "Gmat" ) 
dat <- filter_individuals(dat, type = "PCA")  # Select "X16009_N001" "X16009_N005" "X16009_N006" "X16009_N009" "X16009_N010" "SW"

print(dat, detailed = T) 

# Read and filter the data with conforming markers
dat = read_geno_csv(file="https://raw.githubusercontent.com/jeekinlau/mapping_example/main/data/BExMG_subset_with_contaminants_2075mrks.csv",
                    ploidy = 4, filter.non.conforming = TRUE)
dat <- filter_individuals(dat, ind.to.remove = c("X16009_N001", "X16009_N005", "X16009_N006", "X16009_N009", "X16009_N010", "SW"))


# lets take a look at the distribution of markers
seq.init = make_seq_mappoly(dat, arg="all")
go <- get_genomic_order(input.seq = seq.init) ## get genomic order of the sequence
plot(go)

# lets look at chr 3 because i know there is an inversion in the reference genome
seq_chr <- make_seq_mappoly(seq.init, arg = dat$mrk.names[which(dat$chrom=="3")])

tpt <- est_pairwise_rf(seq_chr,ncpus = 4)
#seq.filt <- rf_snp_filter(tpt, probs = c(0.05, 0.95))
mat <- rf_list_to_matrix(tpt)
mat


seq_test_mds <- mds_mappoly(mat)
seq_mds <- make_seq_mappoly(seq_test_mds)



# the mds based map
seq_mds.map <- est_rf_hmm_sequential(input.seq = seq_mds,
                                  start.set = 5,
                                  thres.twopt = 10, 
                                  thres.hmm = 10,
                                  extend.tail = 30,
                                  info.tail = TRUE, 
                                  twopt = tpt,
                                  sub.map.size.diff.limit = 8, 
                                  phase.number.limit = 20,
                                  reestimate.single.ph.configuration = TRUE,
                                  tol = 10e-3,
                                  tol.final = 10e-4)

plot(seq_mds.map)
plot_map_list(seq_mds.map)
plot_genome_vs_map(list(seq_mds.map))



edit_seq <- edit_order(input.seq = seq_mds)



# inverting the marker genomic positions in the dat object 
# inverted_markers=which(dat$mrk.names %in% edit_seq$inverted)
# inverted_markers_pos = dat$genome.pos[inverted_markers]
# dat2=dat
# dat2$genome.pos[inverted_markers]=as.numeric(rev(inverted_markers_pos))
# edit_seq$data.name="dat2"


edit_seq_mds = make_seq_mappoly(edit_seq)
edit_seq_mds



# map using flipping the inversion and running it using the physical location of the fixed inversion
edit_seq.map <- est_rf_hmm_sequential(input.seq = edit_seq_mds,
                                     start.set = 5,
                                     thres.twopt = 10, 
                                     thres.hmm = 10,
                                     extend.tail = 30,
                                     info.tail = TRUE, 
                                     twopt = tpt,
                                     sub.map.size.diff.limit = 8, 
                                     phase.number.limit = 20,
                                     reestimate.single.ph.configuration = TRUE,
                                     tol = 10e-3,
                                     tol.final = 10e-4)

plot_map_list(list(seq_mds.map,edit_seq.map))
plot_genome_vs_map(list(seq_mds.map,edit_seq.map))



lg3.seq = make_seq_mappoly(dat, arg="seq3")
physorder.map <- est_rf_hmm_sequential(input.seq = lg3.seq,
                                       start.set = 5,
                                       thres.twopt = 10, 
                                       thres.hmm = 10,
                                       extend.tail = 30,
                                       info.tail = TRUE, 
                                       twopt = tpt,
                                       sub.map.size.diff.limit = 8, 
                                       phase.number.limit = 20,
                                       reestimate.single.ph.configuration = TRUE,
                                       tol = 10e-3,
                                       tol.final = 10e-4)

plot_map_list(list(seq_mds.map, edit_seq.map, physorder.map))
plot_genome_vs_map(list(seq_mds.map, edit_seq.map, physorder.map))
summary_maps(list(seq_mds.map,edit_seq.map, physorder.map))









# New faster mapping 

?framework_map

edit_seq.map_framework=framework_map(edit_seq_mds,
                                     twopt = tpt,
                                     start.set = 5)

edit_seq.map_framework_update=update_framework_map(edit_seq.map_framework,input.seq = edit_seq_mds,twopt = tpt)
edit_seq.map_framework_update

plot_map_list(list(seq_mds.map,edit_seq.map,physorder.map,edit_seq.map_framework_update$both$map.list[[5]]))
plot_genome_vs_map(list(seq_mds.map,edit_seq.map,physorder.map,edit_seq.map_framework_update$both$map.list[[5]]))
summary_maps(list(seq_mds.map,edit_seq.map,physorder.map,edit_seq.map_framework_update$both$map.list[[3]]))

?plot_progeny_dosage_change

plot_progeny_dosage_change(list(physorder.map),error = 0.05, output_corrected = T) #caveat doesnt work well yet with MDS or edited sequences





