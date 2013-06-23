

library(BioGeoBEARS)

# Fix -- dmat=dispersal_multipliers_matrix, not matrix of d vals
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_univ_model_v1.R', chdir = TRUE)

#######################################################
# Get the default Psychotria dataset
#######################################################

# Set up data and get tree
# Get the example files directory
extdata_dir = system.file("extdata/", package="BioGeoBEARS")
# tmp hard code: extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"

# Set the filenames (Hawaiian Psychotria from Ree & Smith 2008)
trfn = paste(extdata_dir, "Psychotria_5.2.newick", sep="")
tr = read.tree(file=trfn)

geogfn = paste(extdata_dir, "Psychotria_geog.data", sep="")

# Look at the tree and ranges, for kicks
getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tr



#######################################################
# Set up the run object
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$num_cores_to_use = 2
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01

# Basic info
tips = 1:length(tr$tip.label)
nodes = (length(tr$tip.label)+1):(length(tr$tip.label)+tr$Nnode)

tipranges = getranges_from_LagrangePHYLIP(BioGeoBEARS_run_object$geogfn)
areas = getareas_from_tipranges_object(tipranges)
areas

numstates = numstates_from_numareas(numareas=length(areas))
states_list = areas_list_to_states_list_new(areas)
states_list
states_list_0based_index = rcpp_areas_list_to_states_list(areas)
states_list_0based_index

# Run the ML search
results_object = bears_optim_run(BioGeoBEARS_run_object)
#results_object = model_results
names(results_object)

results_object$ML_marginal_prob_each_state_at_branch_top_AT_node


# calculate the ML marginal probs of states at the base of each branch
# above each split (local, non-joint probs, global model)
results_object = get_MLsplitprobs_from_results(results_object)
names(results_object)

#######################################################
# Get the marginal probs of the splits (global ML probs, not local)
# (These are marginal, rather than joint probs; but not local optima)
#######################################################
leftright_nodes_matrix = get_leftright_nodes_matrix_from_results(tr, results_object)

# This gets you the prob. of each state at the left base above each node, and
# at the right base above each node
marprobs = results_object$ML_marginal_prob_each_split_at_branch_bottom_BELOW_node
left_ML_marginals_by_node = marprobs[leftright_nodes_matrix[, 2], ]
right_ML_marginals_by_node = marprobs[leftright_nodes_matrix[, 1], ]
right_ML_marginals_by_node




#######################################################
# Extract the outputs ancestral states at nodes, and plot!
#######################################################
relprobs_matrix = results_object$ML_marginal_prob_each_state_at_branch_top_AT_node
relprobs_matrix

statenames = areas_list_to_states_list_new(areas, split_ABC=FALSE)
statenames


MLprobs = get_ML_probs(relprobs_matrix)
MLprobs
MLstates = get_ML_states_from_relprobs(relprobs_matrix, statenames, returnwhat="indices")
foo = function(MLstate, statenames) {statenames[MLstate]}
MLstates = unlist(sapply(X=MLstates, FUN=foo, statenames))
MLstates

# Set up colors
colors_matrix = get_colors_for_numareas(length(areas))
colors_list_for_states = mix_colors_for_states(colors_matrix, states_list_0based_index, exclude_null=TRUE)
colors_list_for_states

possible_ranges_list_txt = areas_list_to_states_list_new(areas, split_ABC=FALSE, include_null_range=FALSE)
cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, colors_list_for_states, MLstates)

# Plot to screen
plot(tr, show.tip.label=TRUE, label.offset=0.2)
axisPhylo()
nodelabels(MLstates[nodes], bg=cols_byNode[nodes], cex=0.75)
tiplabels(MLstates[tips], bg=cols_byNode[tips])
title("Hawaiian Psychotria, DEC+J 3-parameter, marginal ML states")
mtext(text="Millions of years ago", side=1, line=2)



#######################################################
# Also add the splits to the plot
#######################################################
# First, get the corner coordinates
coords = corner_coords(tr)
coords

# LEFT SPLITS
relprobs_matrix = left_ML_marginals_by_node
MLprobs = get_ML_probs(relprobs_matrix)
MLprobs
MLstates = get_ML_states_from_relprobs(relprobs_matrix, statenames, returnwhat="indices")
foo = function(MLstate, statenames) {statenames[MLstate]}
MLstates = unlist(sapply(X=MLstates, FUN=foo, statenames))
MLstates

# Set up colors
colors_matrix = get_colors_for_numareas(length(areas))
colors_list_for_states = mix_colors_for_states(colors_matrix, states_list_0based_index, exclude_null=TRUE)
colors_list_for_states

possible_ranges_list_txt = areas_list_to_states_list_new(areas, split_ABC=FALSE, include_null_range=FALSE)
cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, colors_list_for_states, MLstates)

cornerlabels(text=MLstates, coords=coords$leftcorns, bg=cols_byNode, cex=0.75)


# RIGHT SPLITS
relprobs_matrix = right_ML_marginals_by_node
MLprobs = get_ML_probs(relprobs_matrix)
MLprobs
MLstates = get_ML_states_from_relprobs(relprobs_matrix, statenames, returnwhat="indices")
foo = function(MLstate, statenames) {statenames[MLstate]}
MLstates = unlist(sapply(X=MLstates, FUN=foo, statenames))
MLstates

# Set up colors
colors_matrix = get_colors_for_numareas(length(areas))
colors_list_for_states = mix_colors_for_states(colors_matrix, states_list_0based_index, exclude_null=TRUE)
colors_list_for_states

possible_ranges_list_txt = areas_list_to_states_list_new(areas, split_ABC=FALSE, include_null_range=FALSE)
cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, colors_list_for_states, MLstates)

cornerlabels(text=MLstates, coords=coords$rightcorns, bg=cols_byNode, cex=0.75)

