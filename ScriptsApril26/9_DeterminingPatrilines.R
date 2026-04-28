#######################################################################################################################
#************************************* Patriline determination *******************************
#######################################################################################################################
# Estimates the number of patrilines by comparing the paternally assigned haplotypes workers

#Again I'm taking two routes here: 

# - ROUTE 1: Using actual drone information using the true haplotypes to test the accuracy of the threshold method
#            - checks the paternity accuracy really closely but is unrealistic unless you have father-drone information. 
#            - uses sister thresholds and father_thresholds to visualise the most accurate thresholds required

# - ROUTE 2: Use the paternally assigned haplotypes from this script and ONLY sister thresholds to determine patriline numbers
#           - again using multiple sister thresholds to determine which of the Slov real data results are the most reliable. 

################################################################################
#******* Route 1 Functions *******************************

calc_nPaternity_Accuracy <- function(results, sister_threshold, pedigree, father_haplotypes, father_test_threshold) {
  
  # Extract paternal haplotypes
  results_paternal <- rownames(results)[grep(paste0("_paternal"), rownames(results))]
  results_paternal <- results[results_paternal, , drop = FALSE]
  
  # Create a loop grouping together sisters using pedigree 
  queen_ids <- unique(pedigree$mother)
  
  # Initialize a data frame to store the results
  nPaternity <- data.frame(queen_id = queen_ids, 
                           num_fathers_actual = rep(0, length(queen_ids)),
                           num_workers = rep(0, length(queen_ids)),
                           num_fathers_estimated = rep(0, length(queen_ids)),
                           num_fathers_correct = rep(0, length(queen_ids)),
                           father_accuracy_threshold = father_test_threshold,
                           sister_threshold = sister_threshold,
                           stringsAsFactors = FALSE)
  
  for (i in 1:length(queen_ids)) {
    print(i)
    queen_id <- queen_ids[i]
    sister_worker_ids <- t(pedigree[pedigree$mother == queen_id, "id"])
    
    # Pull out the results_paternal with the worker IDs
    pattern <- paste0("^(", paste(sister_worker_ids, collapse = "|"), ")_paternal$")
    matching_indices <- grep(pattern, rownames(results_paternal))
    sister_paternal <- results_paternal[matching_indices, , drop = FALSE]
    
    nPaternity$num_workers[i] <- nrow(sister_paternal)
    
    # Pull out the father IDs for this queen from the pedigree
    father_ids <- t(unique(pedigree[pedigree$mother == queen_id, "father"]))
    
    # Get the actual father haplotypes for this queen
    pattern_fathers <- paste0("^(", paste(father_ids, collapse="|"), ")_paternal$")
    matching_fathers_indices <- grep(pattern_fathers, rownames(father_haplotypes))
    actual_fathers_haplotypes <- father_haplotypes[matching_fathers_indices, , drop = FALSE]
    
    # Ensure column names (alleles) of haplotypes match before comparison
    common_columns <- intersect(colnames(sister_paternal), colnames(actual_fathers_haplotypes))
    
    sister_paternal <- sister_paternal[, common_columns, drop = FALSE]
    actual_fathers_haplotypes <- actual_fathers_haplotypes[, common_columns, drop = FALSE]
    
    # Store the actual number of fathers
    nPaternity$num_fathers_actual[i] <- length(father_ids)
    
    # Create similarity matrix 
    num_worker_haplotypes <- nrow(sister_paternal)
    num_father_haplotypes <- nrow(actual_fathers_haplotypes)
    similarity_matrix <- matrix(0, nrow = num_father_haplotypes, ncol = num_worker_haplotypes)
    rownames(similarity_matrix) <- rownames(actual_fathers_haplotypes)[1:num_father_haplotypes]
    colnames(similarity_matrix) <- rownames(sister_paternal)[1:num_worker_haplotypes]
    
    # Updated calculate_similarity function with NA handling
    calculate_similarity <- function(hap1, hap2) {
      # Remove positions where either haplotype has an NA value
      valid_indices <- !is.na(hap1) & !is.na(hap2)
      
      if (sum(valid_indices) == 0) {
        return(NA)  # Return NA if no valid comparisons can be made
      }
      
      # Calculate similarity only on valid indices
      sum(hap1[valid_indices] == hap2[valid_indices]) / length(hap1[valid_indices])
    }
    
    # Calculate the similarity matrix
    for (j in 1:num_father_haplotypes) {
      for (k in 1:num_worker_haplotypes) {
        haplotype_j <- actual_fathers_haplotypes[j, ]
        haplotype_k <- sister_paternal[k, ]
        match_percentage <- calculate_similarity(haplotype_j, haplotype_k)
        similarity_matrix[j, k] <- match_percentage
      }
    }
    
    # Step 1: Assign workers to father groups based on similarity
    father_groups <- list()
    unmatched_workers <- c()
    
    for (k in 1:num_worker_haplotypes) {
      matching_fathers <- c()
      
      for (j in 1:num_father_haplotypes) {
        if (similarity_matrix[j, k] >= father_test_threshold) {
          matching_fathers <- c(matching_fathers, rownames(actual_fathers_haplotypes)[j])
        }
      }
      
      if (length(matching_fathers) > 1) {
        stop("Threshold is too relaxed, multiple fathers match a single worker.")
      }
      
      if (length(matching_fathers) == 1) {
        father <- matching_fathers[1]
        if (!father %in% names(father_groups)) {
          father_groups[[father]] <- c()
        }
        father_groups[[father]] <- c(father_groups[[father]], rownames(sister_paternal)[k])
      } else {
        unmatched_workers <- c(unmatched_workers, rownames(sister_paternal)[k])
      }
    }
    
    # Step 2: Compare unmatched workers with matched workers in father groups
    still_unmatched <- c()
    
    for (worker in unmatched_workers) {
      matched_to_group <- FALSE
      
      for (father in names(father_groups)) {
        for (sister_worker in father_groups[[father]]) {
          similarity_score <- calculate_similarity(sister_paternal[worker, ], sister_paternal[sister_worker, ])
          if (!is.na(similarity_score) && similarity_score >= sister_threshold) {
            father_groups[[father]] <- c(father_groups[[father]], worker)
            matched_to_group <- TRUE
            break
          }
        }
        if (matched_to_group) {
          break
        }
      }
      
      if (!matched_to_group) {
        still_unmatched <- c(still_unmatched, worker)
      }
    }
    
    # Step 3: Compare still unmatched workers with each other
    num_still_unmatched <- length(still_unmatched)
    new_groups <- list()
    
    if (num_still_unmatched > 1) {
      for (a in 1:(num_still_unmatched - 1)) {
        for (b in (a + 1):num_still_unmatched) {
          worker_a <- still_unmatched[a]
          worker_b <- still_unmatched[b]
          
          similarity_score <- calculate_similarity(sister_paternal[worker_a, ], sister_paternal[worker_b, ])
          if (!is.na(similarity_score) && similarity_score >= sister_threshold) {
            new_groups[[paste0("group_", length(father_groups) + length(new_groups) + 1)]] <- c(worker_a, worker_b)
            still_unmatched <- still_unmatched[!(still_unmatched %in% c(worker_a, worker_b))]
          }
        }
      }
    }
    
    # Assign remaining unmatched workers to their own groups
    for (worker in still_unmatched) {
      new_groups[[paste0("group_", length(father_groups) + length(new_groups) + 1)]] <- c(worker)
    }
    
    # Combine father groups and new groups
    father_groups <- c(father_groups, new_groups)
    
    # Store the number of estimated fathers
    nPaternity$num_fathers_estimated[i] <- length(father_groups)
    
    # Count correct fathers (those that match the actual fathers)
    correct_fathers <- intersect(names(father_groups), rownames(actual_fathers_haplotypes))
    nPaternity$num_fathers_correct[i] <- length(correct_fathers)
  }
  
  return(nPaternity)
}
# Explanation of `calc_nPaternity_Accuracy()`:
{
# 1. Extract paternal haplotypes for all workers.
# 2. Loop through each queen (colony mother) to group her offspring.
# 3. For each queen:
#     - Compare worker haplotypes to true father haplotypes using a similarity threshold.
#     - Build a similarity matrix between all workers and fathers.
#     - Assign workers to father groups if similarity exceeds threshold.
#     - For unassigned workers, check if they resemble any assigned sisters.
#     - Still unmatched workers are compared with each other to form new groups.
#     - Remaining unmatched individuals are placed in their own group.
# 4. Count and store:
#     - Actual number of fathers
#     - Estimated number of fathers from grouping
#     - Number of correctly identified fathers
}

run_paternity_tests <- function(results_arg, sister_thresholds, father_test_thresholds, father_haplotypes) {
  
  # Initialize a list to store the results
  results_list <- list()
  
  # Loop over each combination of sister_threshold and father_test_threshold
  for (sister_threshold in sister_thresholds) {
    for (father_test_threshold in father_test_thresholds) {
      
      # Create a descriptive name for the result
      result_name <- paste0("sis", sister_threshold, "_fat", father_test_threshold)
      
      # Print to ensure thresholds are being passed correctly
      print(paste("Running for sister_threshold:", sister_threshold, "and father_test_threshold:", father_test_threshold))
      
      # Run the function separately for each combination of thresholds
      result <- calc_nPaternity_Accuracy(
        results = results_arg,  # Pass the dynamic results argument here
        sister_threshold = sister_threshold,  # sister threshold
        pedigree = pedigree,
        father_haplotypes = father_haplotypes,
        father_test_threshold = father_test_threshold  # father test threshold
      )
      
      # Store the result with the unique name
      results_list[[result_name]] <- result
      
      # Print to confirm each result has been stored
      print(paste("Finished:", result_name))
    }
  }
  
  # Return the list of results
  return(results_list)
}
# Explanantion of `run_paternity_tests()`:
{
# - Accepts a matrix of paternal haplotypes (`results_arg`), a list of `sister_thresholds`, 
#   and a list of `father_test_thresholds`.
# - For each combination of thresholds, runs `calc_nPaternity_Accuracy()`.
# - Stores results in a named list (`results_list`) for downstream analysis.
}

plot_paternity_number_grid_DRONES <- function(Results) {
  
  # Combine the list of results into one data frame
  df <- dplyr::bind_rows(Results, .id = "threshold_combination")
  
  # Extract thresholds from the list names if not already in the data
  df <- df %>%
    tidyr::separate(threshold_combination, into = c("sis", "fat"), sep = "_") %>%
    dplyr::mutate(
      sister_threshold = as.numeric(gsub("sis", "", sis)),
      father_accuracy_threshold = as.numeric(gsub("fat", "", fat))
    ) %>%
    dplyr::select(-sis, -fat)
  
  # Pivot to long format, excluding num_workers
  df_long <- df %>%
    tidyr::pivot_longer(
      cols = c(num_fathers_actual, num_fathers_estimated),
      names_to = "measure_type",
      values_to = "count"
    ) %>%
    dplyr::mutate(
      measure_type = dplyr::recode(measure_type,
                                   num_fathers_actual = "Actual nPatrilines",
                                   num_fathers_estimated = "Estimated nPatrilines"),
      measure_type = factor(measure_type, levels = c("Actual nPatrilines", "Estimated nPatrilines"))
    )
  
  # Custom labeller function for facets
  custom_labeller <- function(variable, value) {
    if (variable == "sister_threshold") {
      return(paste0("sis_threshold: ", value))
    } else if (variable == "father_accuracy_threshold") {
      return(paste0("drone_threshold: ", value))
    } else {
      return(as.character(value))
    }
  }
  
  # Plot
  plot <- ggplot(df_long, aes(x = as.factor(queen_id), y = count, color = measure_type, shape = measure_type)) +
    geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
    facet_grid(
      rows = vars(father_accuracy_threshold),
      cols = vars(sister_threshold),
      labeller = labeller(
        sister_threshold = function(x) paste0("sis_threshold: ", x),
        father_accuracy_threshold = function(x) paste0("drone_threshold: ", x)
      )
    ) +
    scale_color_manual(name = "Measure Type",
                       values = c("Actual nPatrilines" = "red", "Estimated nPatrilines" = "blue")) +
    scale_shape_manual(name = "Measure Type",
                       values = c("Actual nPatrilines" = 16, "Estimated nPatrilines" = 17)) +
    labs(x = "Colony ID", y = "Number of Patrilines", title = "Estimated vs Actual number of Patrilines") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      strip.text.x = element_text(size = 12),
      strip.text.y = element_text(size = 12),
      strip.placement = "outside",
      strip.background = element_rect(fill = NA),
      plot.title = element_text(size = 18, face = "bold")
    ) +
    labs(
      subtitle = NULL,
      caption = NULL
    ) 
  
  return(plot)
}

################################################################################
#******* Route 2 Functions *******************************

SisterONLY_calc_nPaternity_Accuracy <- function(results, sister_threshold, pedigree, simulated = NULL) {
  
  # Extract paternal haplotypes (still named '_paternal' in the results)
  results_paternal <- rownames(results)[grep(paste0("_paternal"), rownames(results))]
  results_paternal <- results[results_paternal, , drop = FALSE]
  
  # Create a loop grouping together sisters using pedigree 
  queen_ids <- unique(pedigree$mother)
  
  # Initialize a data frame to store the results
  if(simulated == TRUE){
    nPaternity <- data.frame(queen_id = queen_ids, 
                             num_workers = rep(0, length(queen_ids)),
                             num_sister_groups_estimated = rep(0, length(queen_ids)),
                             stringsAsFactors = FALSE,
                             sister_threshold = sister_threshold,
                             actual_number_fathers = rep(0, length(queen_ids)))
  }
  else if(simulated == FALSE){
    nPaternity <- data.frame(queen_id = queen_ids, 
                             num_workers = rep(0, length(queen_ids)),
                             num_sister_groups_estimated = rep(0, length(queen_ids)),
                             stringsAsFactors = FALSE,
                             sister_threshold = sister_threshold)
  }
  
  for (i in 1:length(queen_ids)) {
    print(i)
    queen_id <- queen_ids[i]
    sister_worker_ids <- t(pedigree[pedigree$mother == queen_id, "id"])
    
    if(simulated == TRUE){
      queen_group <- pedigree[pedigree$mother == queen_id,]
      father_ids <- t(unique(pedigree[pedigree$mother == queen_id, "father"]))
      nPaternity$actual_number_fathers[i] <- length(father_ids)
    }
    
    # Pull out the results_paternal with the worker IDs
    pattern <- paste0("^(", paste(sister_worker_ids, collapse = "|"), ")_paternal$")
    matching_indices <- grep(pattern, rownames(results_paternal))
    sister_paternal <- results_paternal[matching_indices, , drop = FALSE]
    
    # Store the number of workers for this queen
    nPaternity$num_workers[i] <- nrow(sister_paternal)
    
    # Ensure we have data to work with
    if (nrow(sister_paternal) == 0) {
      next
    }
    
    # Initialize a list to hold sister groups (clusters)
    sister_groups <- list()
    still_unmatched <- rownames(sister_paternal)
    
    # Define the similarity function to compare workers
    calculate_similarity <- function(hap1, hap2) {
      valid_indices <- !is.na(hap1) & !is.na(hap2)
      
      if (sum(valid_indices) == 0) {
        return(NA)  # Return NA if no valid comparisons can be made
      }
      
      sum(hap1[valid_indices] == hap2[valid_indices]) / length(hap1[valid_indices])
    }
    
    # Step 1: Group workers based on similarity (sister threshold)
    for (a in 1:(length(still_unmatched) - 1)) {
      for (b in (a + 1):length(still_unmatched)) {
        worker_a <- still_unmatched[a]
        worker_b <- still_unmatched[b]
        
        similarity_score <- calculate_similarity(sister_paternal[worker_a, ], sister_paternal[worker_b, ])
        
        # If they meet the sister threshold, group them together
        if (!is.na(similarity_score) && similarity_score >= sister_threshold) {
          group_found <- FALSE
          # Check if worker_a or worker_b already belongs to a group
          for (group_name in names(sister_groups)) {
            if (worker_a %in% sister_groups[[group_name]] || worker_b %in% sister_groups[[group_name]]) {
              sister_groups[[group_name]] <- unique(c(sister_groups[[group_name]], worker_a, worker_b))
              group_found <- TRUE
              break
            }
          }
          # If no group was found, create a new group for them
          if (!group_found) {
            new_group_name <- paste0("group_", length(sister_groups) + 1)
            sister_groups[[new_group_name]] <- c(worker_a, worker_b)
          }
        }
      }
    }
    
    # Step 2: Ensure all remaining workers (if any) are assigned to their own groups
    unmatched_workers <- setdiff(still_unmatched, unlist(sister_groups))
    for (worker in unmatched_workers) {
      new_group_name <- paste0("group_", length(sister_groups) + 1)
      sister_groups[[new_group_name]] <- c(worker)
    }
    
    # Store the number of estimated sister groups (clusters)
    nPaternity$num_sister_groups_estimated[i] <- length(sister_groups)
  }
  
  return(nPaternity)
}
# Explanation `sisterONLY_calc_nPaternity_Accuracy()`:
{
# 1. Extract paternal haplotypes for all workers.
# 2. Loop through each queen (colony mother) to group her offspring.
# 3. For each queen:
#     - (If simulated) record the true number of fathers from the pedigree.
#     - Subset worker paternal haplotypes belonging to that queen.
#     - Compute pairwise similarity between all worker haplotypes.
#     - Group workers into sister clusters if similarity exceeds the sister threshold.
#         * Workers are merged into existing groups if they match any group member.
#         * Otherwise, new groups are created.
#     - Assign any remaining unmatched workers to their own singleton groups.
# 4. Count and store:
#     - Number of workers
#     - Estimated number of sister groups (clusters)
#     - (If simulated) actual number of fathers
}

run_sister_clustering_tests <- function(results_arg, sister_thresholds, pedigree, simulated = NULL) {
  
  # Initialize a list to store the results
  results_list <- list()
  
  # Loop over each sister_threshold
  for (sister_threshold in sister_thresholds) {
    
    # Create a descriptive name for the result
    result_name <- paste0("sis", sister_threshold)
    
    # Print to ensure the threshold is being passed correctly
    print(paste("Running for sister_threshold:", sister_threshold))
    
    # Run the function for each sister_threshold
    result <- SisterONLY_calc_nPaternity_Accuracy(
      results = results_arg,      # Pass the dynamic results argument here
      sister_threshold = sister_threshold,  # sister threshold
      pedigree = pedigree,# Pass the pedigree information
      simulated = simulated
    )
    
    # Store the result with the unique name
    results_list[[result_name]] <- result
    
    # Print to confirm each result has been stored
    print(paste("Finished:", result_name))
  }
  
  # Return the list of results
  return(results_list)
}
# Explanation of `run_sister_clustering_tests()`:
{
# 1. Initialize an empty list to store results.
# 2. Loop through each provided sister similarity threshold.
# 3. For each threshold:
#     - Create a unique name for the result (e.g., "sis0.8").
#     - Print progress to track execution.
#     - Run `SisterONLY_calc_nPaternity_Accuracy()` using:
#         * The input results data
#         * Current sister threshold
#         * Pedigree information
#         * Optional simulation flag
#     - Store the output in the results list under the unique name.
#     - Print confirmation after completion.
# 4. Return the full list of results across all thresholds.
}

plot_paternity_number_grid_SISTERONLY <- function(...) {
  result_sets <- list(...)
  
  # Assign meaningful dataset labels
  names(result_sets) <- c("True_Route2", "nGE_Route2", "GE_Route2")
  
  # Combine all into one dataframe
  df <- purrr::imap_dfr(result_sets, function(results_list, dataset_name) {
    dplyr::bind_rows(results_list, .id = "threshold_combination") %>%
      mutate(dataset = dataset_name)
  })
  
  # Process threshold
  df <- df %>%
    tidyr::separate(threshold_combination, into = c("sis"), sep = "_", fill = "right", extra = "drop") %>%
    dplyr::mutate(
      sister_threshold = as.numeric(gsub("sis", "", sis))
    ) %>%
    dplyr::mutate(
      sister_threshold = factor(sister_threshold, levels = sort(unique(sister_threshold), decreasing = TRUE))
    ) %>%
    dplyr::select(-sis)
  
  # Long format
  df_long <- df %>%
    tidyr::pivot_longer(
      cols = c(actual_number_fathers, num_sister_groups_estimated),
      names_to = "measure_type",
      values_to = "count"
    ) %>%
    dplyr::mutate(
      measure_type = dplyr::recode(measure_type,
                                   actual_number_fathers = "Actual nPatrilines",
                                   num_sister_groups_estimated = "Estimated nPatrilines"),
      measure_type = factor(measure_type, levels = c("Actual nPatrilines", "Estimated nPatrilines"))
    )
  
  # Plot
  plot <- ggplot(df_long, aes(x = as.factor(queen_id), y = count, color = measure_type, shape = measure_type)) +
    geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
    facet_grid(
      rows = vars(dataset),
      cols = vars(sister_threshold),
      labeller = labeller(
        sister_threshold = function(x) paste0("sis_threshold: ", x),
        dataset = function(x) paste0(x)
      )
    ) +
    scale_color_manual(name = "Measure Type",
                       values = c("Actual nPatrilines" = "red", "Estimated nPatrilines" = "blue")) +
    scale_shape_manual(name = "Measure Type",
                       values = c("Actual nPatrilines" = 16, "Estimated nPatrilines" = 17)) +
    labs(x = "Colony ID", y = "Number of Patrilines", title = "Estimated vs Actual number of Patrilines") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      strip.text.x = element_text(size = 12),
      strip.text.y = element_text(size = 12),
      strip.placement = "outside",
      strip.background = element_rect(fill = NA),
      plot.title = element_text(size = 18, face = "bold")
    )
  
  return(plot)
}

################################################################################
#******* ROUTE 1 *******************************
################################################################################
workingDir <- "/home/jana/github/lstrachan_patrilines/"
workingDir = "~/Desktop/lstrachan_patrilines"
setwd(workingDir)

#••••• Simulated ••••• Can't do this with Real data since we don't have the father/drone info 

#Load your Pop object that contains the fathers and the pedigree with the fathers id in there 
Sim_pedigree_mat <- read.csv("Data/worker_pedigree.csv")
fathers_id_all <- Sim_pedigree_mat$father
father_id_1 <- paste(fathers_id_all, "1", sep = "_")

Simulated_pop <- load("Data/Pop_withFathers.Rdata")
true_haplotypes <- pullSnpHaplo(PopMerged)

#since they're coded as identical diploids we can just use one of them
father_haplotypes_1 <- true_haplotypes[rownames(true_haplotypes) %in% father_id_1,]

#replace _1 with _paternal so that it can all be compared 
father_all_unique <- unique(fathers_id_all)
father_id_paternal <- father_all_unique[order(father_all_unique)]
father_id_paternal <- paste(father_id_paternal, "paternal", sep = "_")

father_haplotype_pat <- father_haplotypes_1
rownames(father_haplotype_pat) <- father_id_paternal

#Load the outputs of 7_Haplotype_ParentageAssignments script
# --- Simulated assigned haplotypes from 7_Haplotype_ParentAssignments script 

# Load Rdata from the #7 script <---------------------------------------------------LOOK HERE
Haplo_R1_SimTrue <-  Route1_SimTrue
Haplo_R1_NoGE_SNP2k <-  Route1_NoGE_SNP2k
Haplo_R1_WithGE_SNP2k <-  Route1_WithGE_SNP2k
Haplo_R1_NoGE_SNP50k <-  Route1_NoGE_SNP50k
Haplo_R1_WithGE_SNP50k <-  Route1_WithGE_SNP50k


# Define the thresholds you want to test (takes a good while to run)
sister_thresholds <- c(1.0, 0.95, 0.90, 0.85, 0.80, 0.75)
father_test_thresholds <- c(0.9,0.95, 1)

#Run 2 things: 
#1. Run with haplotypes that have gone through route1 of haplotype parental assignment
#2. Run with haplotypes that have gone through route2 of haplotype parental assignment 

# Haplotypes route 1 : with recon pedigree and both parents 
PatR1_SimTrue_Haplo1 <- run_paternity_tests(results_arg = Haplo_R1_SimTrue$real_results_flipped,
                                     sister_thresholds = sister_thresholds,
                                     father_test_thresholds = father_test_thresholds,
                                     father_haplotypes = father_haplotype_pat)

# Haplotypes route 2 : with maternal pedigree and only maternal info 
PatR1_SimTrue_Haplo2 <- run_paternity_tests(results_arg = Haplo_R2_SimTrue$real_results_flipped,
                                     sister_thresholds = sister_thresholds,
                                     father_test_thresholds = father_test_thresholds,
                                     father_haplotypes = father_haplotype_pat)

#Do this for the others too if it works 


plot_paternity_number_grid_DRONES(PatR1_SimTrue_Haplo1)
plot_paternity_number_grid_DRONES(PatR1_SimTrue_Haplo2)
#Do this for the others too if it works 

################################################################################
#******* ROUTE 2 *******************************
################################################################################

#••• SIMULATED ••• 

#Load the outputs of 7_Haplotype_ParentageAssignments script <--------------------- LOOK HERE 
Haplo_R2_SimTrue <-  Route2_SimTrue
Haplo_R2_NoGE_SNP2k <-  Route2_NoGE_SNP2k
Haplo_R2_WithGE_SNP2k <-  Route2_WithGE_SNP2k
Haplo_R2_NoGE_SNP50k <-  Route2_NoGE_SNP50k
Haplo_R2_WithGE_SNP50k <-  Route2_WithGE_SNP50k

#Run 2 things: 
#1. Run with haplotypes that have gone through route1 of haplotype parental assignment
#2. Run with haplotypes that have gone through route2 of haplotype parental assignment 

# Haplotypes route 1 : with recon pedigree and both parents 

PatR2_SimTrue_Haplo1 <- run_sister_clustering_tests(results_arg = Haplo_R1_SimTrue$real_results_flipped,
                                                    sister_thresholds = sister_thresholds,
                                                    pedigree = simulated_pedigree,
                                                    simulated = TRUE)

# Haplotypes route 2 : with maternal pedigree and only maternal info 
PatR2_SimTrue_Haplo2 <- run_sister_clustering_tests(results_arg = Haplo_R2_SimTrue$real_results_flipped,
                                                    sister_thresholds = sister_thresholds,
                                                    pedigree = simulated_pedigree,
                                                    simulated = TRUE)

#Do this for the others too if it works 

plot_paternity_number_grid_SISTERONLY(
  Dataset1 = PatR2_SimTrue_Haplo1,
  Dataset2 = PatR2_NoGE_SNP2k_Haplo1,
  Dataset3= PatR2_WithGE_SNP2k_Haplo1)

plot_paternity_number_grid_SISTERONLY(
  Dataset1 = PatR2_SimTrue_Haplo2,
  Dataset2 = PatR2_NoGE_SNP2k_Haplo2,
  Dataset3= PatR2_WithGE_SNP2k_Haplo2)



#••• REAL ••• 

PatR2_Real_Haplo1 <- run_sister_clustering_tests(results_arg = Haplo_R1_Real$results$flipped,
                                                 sister_thresholds = sister_thresholds,
                                                 pedigree = Slov_pedigree_rec_filtered,
                                                 simulated = FALSE)


PatR2_Real_Haplo2 <- run_sister_clustering_tests(results_arg = Haplo_R2_Real$results$flipped,
                                                 sister_thresholds = sister_thresholds,
                                                 pedigree = Slov_pedigree_mat_filtered,
                                                 simulated = FALSE)






