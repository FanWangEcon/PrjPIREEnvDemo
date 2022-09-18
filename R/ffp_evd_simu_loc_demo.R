#' simulate demographics distributions by location
#'
#' @description there are m locations and n population groups, simulate the
#'   distribution location-specific population distribution as a m by n
#'   dataframe where the elements of the entire matrix sum jointly to 1.
#'
#' @param it_m_location integer the number of locations
#' @author fan wang, \url{http://fanwangecon.github.io}
#'
#' @return an array of tax-liabilities for particular kids count and martial
#'   status along an array of income levels
#' @references
#' \url{https://fanwangecon.github.io/prjenvdemo/articles/fv_rda_simu_loc_demo.html}
#' @export
#' @import dplyr tibble
#' @importFrom stats dbinom pnorm runif
#' @examples
#' ls_simu_loc_demo <- ffp_evd_simu_loc_demo_main()
#' mt_pop_data_frac <- ls_simu_loc_demo$mt_pop_data_frac
#' mt_loc_pop_bernoulli <- ls_simu_loc_demo$mt_loc_pop_bernoulli
#' df_demo_id_desc <- ls_simu_loc_demo$df_demo_id_desc
#' df_location_id_desc <- ls_simu_loc_demo$df_location_id_desc
#' print("share of population in location and demographic cell")
#' print(mt_pop_data_frac)
#' print(paste0(
#'   "location- and group-specific pbernoulli",
#'   " from pbinom(x|n, pbernoulli_loc_group)"
#' ))
#' print(mt_loc_pop_bernoulli)
#' print("population group id attributes")
#' print(df_demo_id_desc)
#' print("location id attributes")
#' print(df_location_id_desc)
#'
ffp_evd_simu_loc_demo_main <- function(it_M_location = 4,
                                       ar_it_N_pop_groups = c(3, 2),
                                       fl_sd_bound = 1,
                                       ar_st_outer_groups =
                                         c("lower_edu", "middle_edu", "higher_edu"),
                                       ar_st_dims =
                                         c("education", "age_group"),
                                       ar_outer_group_share =
                                         c(0.50, 0.25, 0.25),
                                       ar_fl_withingrp_runif_min =
                                         c(0.30, 0.20, 0.10),
                                       ar_fl_withingrp_runif_max =
                                         c(0.40, 0.80, 0.20),
                                       bl_inner_is_age = TRUE,
                                       it_age_group_year_span = 25,
                                       it_age_youngest_1st_group = 0) {


  # First, construct empty population share dataframe:
  # Assumed variable names for df_id
  it_N_pop_groups <- prod(ar_it_N_pop_groups)

  # Matrix of demographics by location
  mt_pop_data_frac <- matrix(
    data = NA,
    nrow = it_M_location,
    ncol = it_N_pop_groups
  )
  colnames(mt_pop_data_frac) <-
    paste0("popgrp", seq(1, it_N_pop_groups))
  rownames(mt_pop_data_frac) <-
    paste0("location", seq(1, it_M_location))

  mt_loc_pop_bernoulli <- array(NA,
    dim = c(it_M_location, ar_it_N_pop_groups[1])
  )
  colnames(mt_loc_pop_bernoulli) <- ar_st_outer_groups
  rownames(mt_loc_pop_bernoulli) <-
    paste0("location", seq(1, it_M_location))

  # Second, population distribution across locations.
  ar_fl_norm_z <-
    seq(-fl_sd_bound, fl_sd_bound, length.out = it_M_location + 1)
  ar_fl_norm_p_le_z <- pnorm(ar_fl_norm_z, 0, 1)
  ar_fl_norm_p_le_z_diff <- diff(ar_fl_norm_p_le_z)
  ar_fl_p_loc <-
    ar_fl_norm_p_le_z_diff / sum(ar_fl_norm_p_le_z_diff)
  ar_fl_p_loc

  # Third, population distribution within location among demographic (age) groups. Use a vector for binomial distribution. Do this separately for males and females.
  # Different bernoulli "win" probability for each location
  it_pop_inner <- ar_it_N_pop_groups[2]
  for (it_N_pop_outer_group in seq(1, ar_it_N_pop_groups[1])) {
    # Random draws
    set.seed(123 * (it_N_pop_outer_group))
    ar_fl_unif_draw <- sort(runif(it_M_location))

    # Scale with group-specific min and max
    fl_withingrp_runif_min <-
      ar_fl_withingrp_runif_min[it_N_pop_outer_group]
    fl_withingrp_runif_max <-
      ar_fl_withingrp_runif_max[it_N_pop_outer_group]
    fl_withingrp_runif_gap <-
      fl_withingrp_runif_max - fl_withingrp_runif_min
    fl_outer_group_share <-
      ar_outer_group_share[it_N_pop_outer_group]
    ar_fl_unif_prob <-
      sort(runif(it_M_location) * (fl_withingrp_runif_gap) + fl_withingrp_runif_min)
    mt_loc_pop_bernoulli[, it_N_pop_outer_group] <- ar_fl_unif_prob

    # Group-specific column start and end
    it_col_start <- it_pop_inner * (it_N_pop_outer_group - 1) + 1
    it_col_end <- it_pop_inner * (it_N_pop_outer_group)

    # Generate population proportion by locality
    for (it_loc in 1:it_M_location) {
      ar_p_pop_condi_loc <- dbinom(
        0:(it_pop_inner - 1),
        it_pop_inner - 1,
        ar_fl_unif_prob[it_loc]
      )
      mt_pop_data_frac[it_loc, it_col_start:it_col_end] <-
        ar_p_pop_condi_loc * ar_fl_p_loc[it_loc] * fl_outer_group_share
    }
  }

  # Sum of cells, should equal to 1
  # print(paste0("pop frac sum = ", sum(mt_pop_data_frac)))
  mt_pop_data_frac <- mt_pop_data_frac / sum(mt_pop_data_frac)



  # We define some location information, city name, region identifier, etc. Note that location is not the key attribute to show heterogeneities, we will just use a simple "city id" here.
  # Single column location id matrix
  mt_location_id_desc <- matrix(data = NA, nrow = it_M_location, ncol = 1)
  mt_location_id_desc[, 1] <- seq(1, it_M_location)

  # Dataframe with location ID
  set.seed(123)
  ar_st_varnames <- c("location_id")
  df_location_id_desc <- as_tibble(mt_location_id_desc) %>%
    rename_all(~ c(ar_st_varnames)) %>%
    mutate(location = paste0("location", location_id)) %>%
    rowwise() %>%
    mutate(
      countyname =
        paste0(
          "County ",
          paste0(sample(LETTERS, 5, replace = TRUE), collapse = "")
        )
    ) %>%
    select(location_id, location, countyname)


  # We define some population characteristics for the $N$ population groups, specifically, we focus on gender and age groups.
  # 2. Population Age Groups
  # Population Age Groups
  if (bl_inner_is_age) {
    it_age_groups <- ar_it_N_pop_groups[2]
    # Construct start and end ages
    ar_it_start_age <-
      seq(0, it_age_groups - 1) * it_age_group_year_span
    ar_it_end_age <-
      seq(1, it_age_groups) * it_age_group_year_span - 1
    ar_it_start_age <- ar_it_start_age + it_age_youngest_1st_group
    ar_it_end_age <- ar_it_end_age + it_age_youngest_1st_group
    # Age group strings
    ar_st_age_groups <-
      paste0("age_", ar_it_start_age, "_to_", ar_it_end_age)
  }

  # 3. use the base function expand.grid to mesh the two string arrays together
  mt_demo_id_desc <- expand.grid(
    sub_group = ar_st_age_groups,
    super_group = ar_st_outer_groups
  )

  # 4. Three columns, the first column is the group ID, convert to dateframe
  ar_st_varnames <- c("popgrp_id", ar_st_dims)
  df_demo_id_desc <- as_tibble(mt_demo_id_desc) %>%
    rowid_to_column(var = "id") %>%
    rename_all(~ c(ar_st_varnames)) %>%
    mutate(popgrp = paste0("popgrp", popgrp_id)) %>%
    select(popgrp_id, popgrp, ar_st_dims)


  # return
  return(
    list(
      mt_pop_data_frac = mt_pop_data_frac,
      mt_loc_pop_bernoulli = mt_loc_pop_bernoulli,
      df_demo_id_desc = df_demo_id_desc,
      df_location_id_desc = df_location_id_desc
    )
  )
}
