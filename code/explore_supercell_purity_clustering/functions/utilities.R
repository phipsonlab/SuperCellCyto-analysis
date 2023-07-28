#' Get a list of clustering result directories
#'
#' @param root_dir
#' @param gamma_values
#' @param clust_algo
#'
#' @return
#' @export
#'
#' @examples
get_supercell_clustering_res_directories <- function(root_dir, gamma_values, clust_algo = c("louvain", "flowsom")) {
    clust_algo <- match.arg(clust_algo)

    res_dirnames <- lapply(gamma_values, function(gamma_value) {
        here(
            root_dir, paste0(clust_algo, "_supercell_runs"),
            paste0("gamma_", gamma_value)
        )
    })
    names(res_dirnames) <- sapply(gamma_values, function(g) {
        paste0("supercell_gamma", g, "_", clust_algo, "_savedir")
    })

    return(res_dirnames)
}


