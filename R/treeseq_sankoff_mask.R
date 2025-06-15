#' Compute geographic ancestry coefficients through time (masked)
#'
#' @description
#' For each combination of sample subset and geographic region, calculates the 
#' proportion of genetic material inherited from ancestors in that region at 
#' different points in time, masking a subset of trees. This allows tracking how the geographic distribution 
#' of ancestry changes as we move backwards in time.
#'
#' @param ts A \code{treeseq} object
#' @param obj Result object from \code{treeseq_discrete_mpr}
#' @param cost_matrix Symmetric numeric matrix of migration costs between states
#' @param adjacency_matrix Binary matrix specifying allowed transitions between states
#' @param times Numeric vector of time points at which to calculate ancestry,
#'   must be ordered from present (0) to past
#' @param state_sets Integer vector grouping states into regions. Length must match
#'   number of states, values indicate region membership (1-based)
#' @param sample_sets Integer vector grouping samples into subsets. Length must match
#'   number of samples, values indicate subset membership (1-based). Use 0 to exclude
#'   samples.
#' @param tree_mask Integer (or logical) vector across all trees, where (0) means 
#'   the tree should be skipped, and one (1) means the tree should be included
#'
#' @return A 3-dimensional array with dimensions:
#'   [region, sample_subset, time_point]
#' Values represent the proportion of (non-masked) genetic material that sample_subset inherits
#' from ancestors in region at each time_point.
#'
#' @details
#' This function traces genetic material backwards in time to determine its
#' geographic location at different time points. For each sample subset, time point,
#' and geographic region, it calculates the fraction of the genome that was located
#' in that region at that time.
#'
#' States can be grouped into regions using state_sets. For example, multiple
#' states might represent different locations within the same continent. Similarly,
#' samples can be grouped into subsets using sample_sets, allowing calculation
#' of ancestry coefficients for different populations or sampling locations.
#'
#' The function uses the migration paths sampled by edge_history to determine
#' locations of genetic material through time. Different random samplings of
#' migration paths may give slightly different results when multiple equally
#' parsimonious paths exist.
#'
#' @seealso
#' \code{\link{treeseq_discrete_mpr_edge_history}} for the underlying migration paths
#' \code{\link{treeseq_discrete_mpr_ancestry_flux}} for tracking migration between regions
#'
#' @examples
#' # Load tree sequence
#' ts = treeseq_load(system.file("extdata", "test.trees", package="gaia"))
#'
#' # Set up states, costs, and adjacency
#' state = c(2L,1L,1L)
#' samples = cbind(node_id=0:2, state_id=state)
#' costs = matrix(c(0,1,1,0), 2, 2)
#' adjacency = matrix(1, 2, 2)
#' diag(adjacency) = 0
#'
#' # Compute base MPR
#' mpr_costs = treeseq_discrete_mpr(ts, samples, costs)
#'
#' # Define timepoints to examine ancestry
#' # Must cover node times from 0 to 1.0
#' times = seq(0, 1.0, by=0.2)
#'
#' # Compute ancestry coefficients
#' ancestry = treeseq_discrete_mpr_ancestry(
#'   ts, mpr_costs, costs, adjacency, times, 
#'   state_sets=1:2,     # Keep states separate
#'   sample_sets=1:3     # Keep samples separate
#' )
treeseq_discrete_mpr_ancestry_mask = function(ts, obj, cost_matrix,
    adjacency_matrix, times, state_sets, sample_sets, tree_mask)
{
    stopifnot(inherits(ts, "treeseq"))
    stopifnot(inherits(obj, "discrete") && inherits(obj, "mpr"))
    stopifnot(!is.unsorted(times))
    stopifnot(all(times >= 0))
    num_states = nrow(cost_matrix)
    num_samples = treeseq_num_samples(ts)
    if (missing(state_sets))
        state_sets = 1:num_states
    if (missing(sample_sets))
        sample_sets = 1:num_samples
    stopifnot(is.integer(state_sets))
    stopifnot(is.integer(sample_sets))
    stopifnot(all(state_sets > 0))
    stopifnot(all(sample_sets >= 0))
    stopifnot(all(tabulate(state_sets) > 0))
    stopifnot(all(tabulate(sample_sets) > 0))
    stopifnot(length(state_sets) == num_states)
    stopifnot(length(sample_sets) == num_samples)
    if (is.logical(tree_mask) || is.numeric(tree_mask))
        tree_mask = as.integer(tree_mask)
    stopifnot(is.integer(tree_mask))
    e = treeseq_discrete_mpr_edge_history(
        ts, obj, cost_matrix, adjacency_matrix, FALSE)
    .Call(
        C_treeseq_discrete_mpr_ancestry_mask
        , ts@tree
        , attr(e, "path.offset")
        , e$state_id
        , e$time
        , attr(e, "node.state")
        , max(state_sets)
        , state_sets - 1L
        , max(sample_sets)
        , sample_sets - 1L
        , times
        , tree_mask
    )
}


# #' Compute migration flux between geographic regions through time
# #'
# #' @description
# #' Tracks the movement of genetic material between geographic regions over time by
# #' computing the proportion of each sample subset's genome that migrated between
# #' each pair of regions during different time intervals. (masked)
# #'
# #' @param ts A \code{treeseq} object
# #' @param obj Result object from \code{treeseq_discrete_mpr}
# #' @param cost_matrix Symmetric numeric matrix of migration costs between states
# #' @param adjacency_matrix Binary matrix specifying allowed transitions between states
# #' @param times Numeric vector of time points defining intervals for flux calculation,
# #'   must be ordered from present (0) to past
# #' @param state_sets Integer vector grouping states into regions. Length must match
# #'   number of states, values indicate region membership (1-based)
# #' @param sample_sets Integer vector grouping samples into subsets. Length must match
# #'   number of samples, values indicate subset membership (1-based). Use 0 to exclude
# #'   samples.
# #' @param tree_mask Integer (or logical) vector across all trees, where (0) means 
# #'   the tree should be skipped, and one (1) means the tree should be included
# #'
# #' @return A 4-dimensional array with dimensions:
# #'   [source_region, dest_region, sample_subset, time_interval]
# #' Values represent the proportion of sample_subset's genome that moved from
# #' source_region to dest_region during each time_interval.
# #'
# #' @details
# #' This function extends ancestry coefficient calculations by explicitly tracking
# #' migrations between regions. For each time interval, it identifies genetic material
# #' that changed regions during that interval and records the source and destination
# #' regions.
# #'
# #' States can be grouped into regions using state_sets. For example, multiple
# #' states might represent different locations within the same continent. Similarly,
# #' samples can be grouped into subsets using sample_sets, allowing calculation
# #' of migration flux for different populations or sampling locations.
# #'
# #' The function uses the migration paths sampled by edge_history to determine
# #' when migrations occurred. Different random samplings of migration paths may
# #' give slightly different results when multiple equally parsimonious paths exist.
# #'
# #' Time intervals are defined by consecutive pairs of values in the times parameter.
# #' The flux value for an interval represents all migrations that occurred between
# #' the start and end of that interval.
# #'
# #' @seealso
# #' \code{\link{treeseq_discrete_mpr_edge_history}} for the underlying migration paths
# #' \code{\link{treeseq_discrete_mpr_ancestry}} for static ancestry proportions
# #'
# #' @examples
# #' # Load tree sequence
# #' ts = treeseq_load(system.file("extdata", "test.trees", package="gaia"))
# #'
# #' # Set up states, costs, and adjacency 
# #' state = c(2L,1L,1L)
# #' samples = cbind(node_id=0:2, state_id=state)
# #' costs = matrix(c(0,1,1,0), 2, 2)
# #' adjacency = matrix(1, 2, 2)
# #' diag(adjacency) = 0
# #'
# #' # Compute base MPR
# #' mpr_costs = treeseq_discrete_mpr(ts, samples, costs)
# #'
# #' # Define timepoints for flux intervals
# #' # Times should span node times (0-1.0)
# #' times = seq(0, 1.0, by=0.2)
# #'
# #' # Compute migration flux between states over time
# #' flux = treeseq_discrete_mpr_ancestry_flux(
# #'   ts, mpr_costs, costs, adjacency, times,
# #'   state_sets=1:2,     # Keep states separate
# #'   sample_sets=1:3     # Keep samples separate
# #' )
# treeseq_discrete_mpr_ancestry_flux_mask = function(ts, obj, cost_matrix,
#     adjacency_matrix, times, state_sets, sample_sets, tree_mask)
# {
#     stopifnot(inherits(ts, "treeseq"))
#     stopifnot(inherits(obj, "discrete") && inherits(obj, "mpr"))
#     stopifnot(!is.unsorted(times))
#     stopifnot(all(times >= 0))
#     num_states = nrow(cost_matrix)
#     num_samples = treeseq_num_samples(ts)
#     if (missing(state_sets))
#         state_sets = 1:num_states
#     if (missing(sample_sets))
#         sample_sets = rep(1L, num_samples)
#     stopifnot(is.integer(state_sets))
#     stopifnot(is.integer(sample_sets))
#     stopifnot(all(state_sets > 0))
#     stopifnot(all(sample_sets >= 0))
#     stopifnot(all(tabulate(state_sets) > 0))
#     stopifnot(all(tabulate(sample_sets) > 0))
#     stopifnot(length(state_sets) == num_states)
#     stopifnot(length(sample_sets) == num_samples)
#     num_state_sets = as.numeric(max(state_sets))
#     num_sample_sets = as.numeric(max(sample_sets))
#     num_time_bins = length(times) - 1
#     if (is.logical(tree_mask) || is.numeric(tree_mask))
#         tree_mask = as.integer(tree_mask)
#     stopifnot(is.integer(tree_mask))
#     storage = num_state_sets * num_state_sets * num_sample_sets * num_time_bins
#     if (storage > .Machine$integer.max)
#         stop("storage requirements too large")
#     e = treeseq_discrete_mpr_edge_history(
#         ts, obj, cost_matrix, adjacency_matrix, FALSE)
#     .Call(
#         C_treeseq_discrete_mpr_ancestry_flux
#         , ts@tree
#         , attr(e, "path.offset")
#         , e$state_id
#         , e$time
#         , as.integer(num_state_sets)
#         , state_sets - 1L
#         , as.integer(num_sample_sets)
#         , sample_sets - 1L
#         , times
#         , tree_mask
#     )
# }