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

calc_nPaternity_Accuracy <- function(results, sister_threshold, pedigree, father_haplotypes, father_test_threshold, recon_pedigree,data_type= NULL, haplo_assignment_type= NULL) {
  
  # Extract paternal haplotypes
  results_paternal <- rownames(results)[grep(paste0("_paternal"), rownames(results))]
  results_paternal <- results[results_paternal, , drop = FALSE]
  
  # Create a loop grouping together sisters using pedigree 
  queen_ids <- unique(pedigree$mother)
  
  # Initialize a data frame to store the results
  nPaternity <- data.frame(queen_id = queen_ids, 
                           num_dpqs_actual = rep(0, length(queen_ids)),
                           num_dpqs_recon_ped = rep(0, length(queen_ids)),
                           num_fathers_actual = rep(0, length(queen_ids)),
                           num_workers_actual = rep(0, length(queen_ids)),
                           num_workers_recon = rep(0, length(queen_ids)),
                           num_fathers_estimated = rep(0, length(queen_ids)),
                           num_fathers_correct = rep(0, length(queen_ids)),
                           father_accuracy_threshold = father_test_threshold,
                           sister_threshold = sister_threshold,
                           data_type = data_type,
                           haplo_assignment_type = haplo_assignment_type,
                           stringsAsFactors = FALSE)
  
  for (i in 1:length(queen_ids)) {
    print(i)
    queen_id <- queen_ids[i]
    sister_worker_ids <- t(pedigree[pedigree$mother == queen_id, "id"])
    nPaternity$num_workers_actual[i]  <- nrow(pedigree[pedigree$mother == queen_id, "id"])
    num_dpqs_actual <- pedigree[pedigree$mother == queen_id,]
    
    nPaternity$num_dpqs_actual[i] <- length(unique(num_dpqs_actual$dpc))
    
    # Check if recon_pedigree is missing, null, or the string "NA"
    if (is.null(recon_pedigree) || (length(recon_pedigree) == 1 && is.na(recon_pedigree))) {
      nPaternity$num_dpqs_recon_ped[i] <- NA
      nPaternity$num_workers_recon[i] <- NA
    } else {
      # It is a dataframe, proceed with filtering
      matches <- recon_pedigree[recon_pedigree$mother == queen_id, ]
      nPaternity$num_dpqs_recon_ped[i] <- length(unique(matches$dpc))
      nPaternity$num_workers_recon[i] <- length(unique(matches$id))
    }
    
    # Pull out the results_paternal with the worker IDs
    pattern <- paste0("^(", paste(sister_worker_ids, collapse = "|"), ")_paternal$")
    matching_indices <- grep(pattern, rownames(results_paternal))
    sister_paternal <- results_paternal[matching_indices, , drop = FALSE]
  
    
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

run_paternity_tests <- function(results_arg, sister_thresholds, father_test_thresholds, father_haplotypes, pedigree, recon_pedigree, data_type= NULL, haplo_assignment_type= NULL) {
  
  results_list <- list()
  
  for (sister_threshold in sister_thresholds) {
    for (father_test_threshold in father_test_thresholds) {
      
      result_name <- paste0("sis", sister_threshold, "_fat", father_test_threshold)
      print(paste("Running:", result_name))
      
      # Wrap the core logic in tryCatch to prevent crashing
      result <- tryCatch({
        calc_nPaternity_Accuracy(
          results = results_arg,
          sister_threshold = sister_threshold,
          pedigree = pedigree,
          recon_pedigree = recon_pedigree,
          father_haplotypes = father_haplotypes,
          father_test_threshold = father_test_threshold,
          data_type = data_type,
          haplo_assignment_type = haplo_assignment_type
          )
      }, error = function(e) {
        # This code runs only if an error (like 'stop') occurs
        message(paste("!!! Error in", result_name, ":", e$message))
        return(paste("FAILED:", e$message)) # Return a string instead of a dataframe
      })
      
      results_list[[result_name]] <- result
      print(paste("Completed:", result_name))
    }
  }
  
  return(results_list)
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

SisterONLY_calc_nPaternity_Accuracy <- function(results, sister_threshold, pedigree, recon_pedigree, simulated = NULL, data_type= NULL, haplo_assignment_type= NULL) {
  
  # Extract paternal haplotypes (still named '_paternal' in the results)
  results_paternal <- rownames(results)[grep(paste0("_paternal"), rownames(results))]
  results_paternal <- results[results_paternal, , drop = FALSE]
  
  # Create a loop grouping together sisters using pedigree 
  queen_ids <- unique(pedigree$mother)
  
  # Initialize a data frame to store the results
  if(simulated == TRUE){
    nPaternity <- data.frame(queen_id = queen_ids, 
                             num_dpqs_actual = rep(0, length(queen_ids)),
                             num_dpqs_recon_ped = rep(0, length(queen_ids)),
                             num_workers_actual = rep(0, length(queen_ids)),
                             num_workers_recon = rep(0, length(queen_ids)),
                             num_sister_groups_estimated = rep(0, length(queen_ids)),
                             stringsAsFactors = FALSE,
                             sister_threshold = sister_threshold,
                             actual_number_fathers = rep(0, length(queen_ids)),
                             data_type = data_type,
                             haplo_assignment_type = haplo_assignment_type)
  }
  if(simulated == FALSE){
    nPaternity <- data.frame(queen_id = queen_ids, 
                             num_dpqs_recon_ped = rep(0, length(queen_ids)),
                             num_workers_actual = rep(0, length(queen_ids)),
                             num_workers_recon = rep(0, length(queen_ids)),
                             num_sister_groups_estimated = rep(0, length(queen_ids)),
                             stringsAsFactors = FALSE,
                             sister_threshold = sister_threshold,
                             data_type = data_type,
                             haplo_assignment_type = haplo_assignment_type)
  }
  
  for (i in 1:length(queen_ids)) {
    print(i)
    queen_id <- queen_ids[i]
    sister_worker_ids <- t(pedigree[pedigree$mother == queen_id, "id"])
    if (simulated == TRUE){
      nPaternity$num_workers_actual[i]  <- nrow(pedigree[pedigree$mother == queen_id, "id"])    
      }
    if (simulated == FALSE){
    nPaternity$num_workers_actual[i]  <- length(sister_worker_ids)
     }
    
    if(simulated == TRUE){
      queen_group <- pedigree[pedigree$mother == queen_id,]
      father_ids <- t(unique(pedigree[pedigree$mother == queen_id, "father"]))
      nPaternity$actual_number_fathers[i] <- length(father_ids)
      nPaternity$num_dpqs_actual[i] <- length(unique(queen_group$dpc))
    }
    
    # Check if recon_pedigree is missing, null, or the string "NA"
    if (is.null(recon_pedigree) || (length(recon_pedigree) == 1 && is.na(recon_pedigree))) {
      nPaternity$num_dpqs_recon_ped[i] <- NA
      nPaternity$num_workers_recon[i] <- NA
    } else {
      # It is a dataframe, proceed with filtering
      matches <- recon_pedigree[recon_pedigree$mother == queen_id, ]
      nPaternity$num_dpqs_recon_ped[i] <- length(unique(matches$dpc))
      nPaternity$num_workers_recon[i] <- length(unique(matches$id))
    }
    
    # Pull out the results_paternal with the worker IDs
    pattern <- paste0("^(", paste(sister_worker_ids, collapse = "|"), ")_paternal$")
    matching_indices <- grep(pattern, rownames(results_paternal))
    sister_paternal <- results_paternal[matching_indices, , drop = FALSE]
    
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

run_sister_clustering_tests <- function(results_arg, sister_thresholds, pedigree, recon_pedigree, simulated = NULL, data_type= NULL, haplo_assignment_type= NULL) {
  
  # Initialize a list to store the results
  results_list <- list()
  
  # Loop over each sister_threshold
  for (sister_threshold in sister_thresholds) {
    
    # Create a descriptive name for the result
    result_name <- paste0("sis", sister_threshold)
    
    # Print to ensure the threshold is being passed correctly
    print(paste("Running for sister_threshold:", sister_threshold))
    
    # Run the function for each sister_threshold
    result <- tryCatch({
      SisterONLY_calc_nPaternity_Accuracy(
      results = results_arg,      # Pass the dynamic results argument here
      sister_threshold = sister_thresholds,  # sister threshold
      pedigree = pedigree,# Pass the pedigree information
      recon_pedigree = recon_pedigree, #Pass on reconstructed pedigree info
      simulated = simulated,
      data_type= data_type, 
      haplo_assignment_type= haplo_assignment_type
    )
    }, error = function(e) {
      # This code runs only if an error (like 'stop') occurs
      message(paste("!!! Error in", result_name, ":", e$message))
      return(paste("FAILED:", e$message)) # Return a string instead of a dataframe
    })
    # Store the result with the unique name
    results_list[[result_name]] <- result
    
    # Print to confirm each result has been stored
    print(paste("Finished:", result_name))
  }
  
  # Return the list of results
  return(results_list)
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

#Load the outputs of 7_Haplotype_ParentageAssignments script
load("Data/Pipeline/7_Haplotype_ParentAssignments.Rdata")

#••••• Simulated ••••• Can't do this with Real data since we don't have the father/drone info 

#Load your Pop object that contains the fathers and the pedigree with the fathers id in there 
Sim_pedigree_mat <- read.csv("Data/worker_pedigree.csv")
fathers_id_all <- Sim_pedigree_mat$father
father_id_1 <- paste(fathers_id_all, "1", sep = "_")

true_haplotypes <- pullSnpHaplo(PopMerged)

#since they're coded as identical diploids we can just use one of them
father_haplotypes_1 <- true_haplotypes[rownames(true_haplotypes) %in% father_id_1,]

#replace _1 with _paternal so that it can all be compared 
father_all_unique <- unique(fathers_id_all)
father_id_paternal <- father_all_unique[order(father_all_unique)]
father_id_paternal <- paste(father_id_paternal, "paternal", sep = "_")

father_haplotype_pat <- father_haplotypes_1
rownames(father_haplotype_pat) <- father_id_paternal

# Define the thresholds you want to test (takes a good while to run)
sister_thresholds <- c(1.0, 0.95, 0.90, 0.85, 0.80, 0.75)
father_test_thresholds <- c(1,0.95, 0.9)

#Run 2 things: 
#1. Run with haplotypes that have gone through route1 of haplotype parental assignment
#2. Run with haplotypes that have gone through route2 of haplotype parental assignment 

# Haplotypes route 1 : with recon pedigree and both parents 
PatR1_SimTrue_Haplo1 <- run_paternity_tests(results_arg = Route1_SimTrue$real_results_flipped,
                                               sister_thresholds = sister_thresholds,
                                               father_test_thresholds = father_test_thresholds,
                                               pedigree = Sim_pedigree_mat,
                                               recon_pedigree = NULL,
                                               father_haplotypes = father_haplotype_pat,
                                            data_type = "True",
                                            haplo_assignment_type = "Used recon pedigree")


PatR1_NoGE_SNP2k_Haplo1 <- run_paternity_tests(results_arg = Route1_NoGE_SNP2k$results_flipped,
                                     sister_thresholds = sister_thresholds,
                                     father_test_thresholds = father_test_thresholds,
                                     pedigree = Sim_pedigree_mat,
                                     recon_pedigree = Rec_pedigree_2k_NoGE_filtered,
                                     father_haplotypes = father_haplotype_pat,
                                     data_type = "NoGE_2k",
                                     haplo_assignment_type = "Used recon pedigree")


PatR1_NoGE_SNP50k_Haplo1 <- run_paternity_tests(results_arg = Route1_NoGE_SNP50k$results_flipped,
                                               sister_thresholds = sister_thresholds,
                                               father_test_thresholds = father_test_thresholds,
                                               pedigree = Sim_pedigree_mat,
                                               recon_pedigree = Rec_pedigree_50k_NoGE_filtered,
                                               father_haplotypes = father_haplotype_pat,
                                               data_type = "NoGE_50k",
                                               haplo_assignment_type = "Used recon pedigree")

PatR1_WithGE_SNP2k_Haplo1 <- run_paternity_tests(results_arg = Route1_WithGE_SNP2k$results_flipped,
                                               sister_thresholds = sister_thresholds,
                                               father_test_thresholds = father_test_thresholds,
                                               pedigree = Sim_pedigree_mat,
                                               recon_pedigree = Rec_pedigree_2k_WithGE_filtered,
                                               father_haplotypes = father_haplotype_pat,
                                               data_type = "WithGE_SNP2k",
                                               haplo_assignment_type = "Used recon pedigree")


PatR1_WithGE_SNP50k_Haplo1 <- run_paternity_tests(results_arg = Route1_WithGE_SNP50k$results_flipped,
                                                sister_thresholds = sister_thresholds,
                                                father_test_thresholds = father_test_thresholds,
                                                pedigree = Sim_pedigree_mat,
                                                recon_pedigree = Rec_pedigree_50k_WithGE_filtered,
                                                father_haplotypes = father_haplotype_pat,
                                                data_type = "WithGE_50k",
                                                haplo_assignment_type = "Used recon pedigree")

# Haplotypes route 2 : with maternal pedigree and only maternal info 
PatR1_SimTrue_Haplo2 <- run_paternity_tests(results_arg = Route2_SimTrue$real_results_flipped,
                                            sister_thresholds = sister_thresholds,
                                            father_test_thresholds = father_test_thresholds,
                                            pedigree = Sim_pedigree_mat,
                                            recon_pedigree = NULL,
                                            father_haplotypes = father_haplotype_pat,
                                            data_type = "True",
                                            haplo_assignment_type = "Used mat pedigree")


PatR1_NoGE_SNP2k_Haplo2 <- run_paternity_tests(results_arg = Route2_NoGE_SNP2k$phased_results_flipped,
                                               sister_thresholds = sister_thresholds,
                                               father_test_thresholds = father_test_thresholds,
                                               pedigree = Sim_pedigree_mat,
                                               recon_pedigree = NULL,
                                               father_haplotypes = father_haplotype_pat,
                                               data_type = "NoGE_2k",
                                               haplo_assignment_type = "Used mat pedigree")


PatR1_NoGE_SNP50k_Haplo2 <- run_paternity_tests(results_arg = Route2_NoGE_SNP50k$phased_results_flipped,
                                                sister_thresholds = sister_thresholds,
                                                father_test_thresholds = father_test_thresholds,
                                                pedigree = Sim_pedigree_mat,
                                                recon_pedigree = NULL,
                                                father_haplotypes = father_haplotype_pat,
                                                data_type = "NoGE_50k",
                                                haplo_assignment_type = "Used mat pedigree")

PatR1_WithGE_SNP2k_Haplo2 <- run_paternity_tests(results_arg = Route2_WithGE_SNP2k$phased_results_flipped,
                                                 sister_thresholds = sister_thresholds,
                                                 father_test_thresholds = father_test_thresholds,
                                                 pedigree = Sim_pedigree_mat,
                                                 recon_pedigree = NULL,
                                                 father_haplotypes = father_haplotype_pat,
                                                 data_type = "WithGE_SNP2k",
                                                 haplo_assignment_type = "Used mat pedigree")


PatR1_WithGE_SNP50k_Haplo2 <- run_paternity_tests(results_arg = Route2_WithGE_SNP50k$phased_results_flipped,
                                                  sister_thresholds = sister_thresholds,
                                                  father_test_thresholds = father_test_thresholds,
                                                  pedigree = Sim_pedigree_mat,
                                                  recon_pedigree = NULL,
                                                  father_haplotypes = father_haplotype_pat,
                                                  data_type = "WithGE_SNP50k",
                                                  haplo_assignment_type = "Used mat pedigree")




# plot_paternity_number_grid_DRONES(PatR1_SimTrue_Haplo1)
# plot_paternity_number_grid_DRONES(PatR1_SimTrue_Haplo2)
#Do this for the others too if it works 

################################################################################
#******* ROUTE 2 *******************************
################################################################################

#••• SIMULATED ••• 

#Run 2 things: 
#1. Run with haplotypes that have gone through route1 of haplotype parental assignment
#2. Run with haplotypes that have gone through route2 of haplotype parental assignment 

# Haplotypes route 1 : with recon pedigree and both parents 

sister_thresholds <- c(1.0, 0.95, 0.90, 0.85, 0.80, 0.75, 0.7)

PatR2_SimTrue_Haplo1 <- run_sister_clustering_tests(results_arg = Route1_SimTrue$real_results_flipped,
                                                       sister_thresholds = sister_thresholds,
                                                       pedigree = Sim_pedigree_mat,
                                                       recon_pedigree = NULL,
                                                       simulated = TRUE,
                                                    data_type = "True",
                                                    haplo_assignment_type = "Used recon pedigree")

PatR2_NoGE_SNP2k_Haplo1 <- run_sister_clustering_tests(results_arg = Route1_NoGE_SNP2k$results_flipped,
                                                    sister_thresholds = sister_thresholds,
                                                    pedigree = Sim_pedigree_mat,
                                                    recon_pedigree = Rec_pedigree_2k_NoGE_filtered,
                                                    simulated = TRUE,
                                                    data_type = "NoGE_SNP2k",
                                                    haplo_assignment_type = "Used recon pedigree")

PatR2_NoGE_SNP50k_Haplo1 <- run_sister_clustering_tests(results_arg = Route1_NoGE_SNP50k$results_flipped,
                                                       sister_thresholds = sister_thresholds,
                                                       pedigree = Sim_pedigree_mat,
                                                       recon_pedigree = Rec_pedigree_50k_NoGE_filtered,
                                                       simulated = TRUE,
                                                       data_type = "NoGE_SNP50k",
                                                       haplo_assignment_type = "Used recon pedigree")

PatR2_WithGE_SNP2k_Haplo1 <- run_sister_clustering_tests(results_arg = Route1_WithGE_SNP2k$results_flipped,
                                                       sister_thresholds = sister_thresholds,
                                                       pedigree = Sim_pedigree_mat,
                                                       recon_pedigree = Rec_pedigree_2k_WithGE_filtered,
                                                       simulated = TRUE,
                                                       data_type = "WithGE_SNP2k",
                                                       haplo_assignment_type = "Used recon pedigree")

PatR2_WithGE_SNP50k_Haplo1 <- run_sister_clustering_tests(results_arg = Route1_WithGE_SNP50k$results_flipped,
                                                        sister_thresholds = sister_thresholds,
                                                        pedigree = Sim_pedigree_mat,
                                                        recon_pedigree = Rec_pedigree_50k_WithGE_filtered,
                                                        simulated = TRUE,
                                                        data_type = "WithGE_SNP50k",
                                                        haplo_assignment_type = "Used recon pedigree")

# Haplotypes route 2 : with maternal pedigree and only maternal info 
PatR2_SimTrue_Haplo2 <- run_sister_clustering_tests(results_arg = Route2_SimTrue$real_results_flipped,
                                                    sister_thresholds = sister_thresholds,
                                                    pedigree = Sim_pedigree_mat,
                                                    recon_pedigree = NULL,
                                                    simulated = TRUE,
                                                    data_type = "True",
                                                    haplo_assignment_type = "Used mat pedigree")

PatR2_NoGE_SNP2k_Haplo2 <- run_sister_clustering_tests(results_arg = Route2_NoGE_SNP2k$phased_results_flipped,
                                                       sister_thresholds = sister_thresholds,
                                                       pedigree = Sim_pedigree_mat,
                                                       recon_pedigree = NULL,
                                                       simulated = TRUE,
                                                       data_type = "NoGE_2k",
                                                       haplo_assignment_type = "Used mat pedigree")

PatR2_NoGE_SNP50k_Haplo2 <- run_sister_clustering_tests(results_arg = Route2_NoGE_SNP50k$phased_results_flipped,
                                                        sister_thresholds = sister_thresholds,
                                                        pedigree = Sim_pedigree_mat,
                                                        recon_pedigree = NULL,
                                                        simulated = TRUE,
                                                        data_type = "NoGE_50k",
                                                        haplo_assignment_type = "Used mat pedigree")

PatR2_WithGE_SNP2k_Haplo2 <- run_sister_clustering_tests(results_arg = Route2_WithGE_SNP2k$phased_results_flipped,
                                                         sister_thresholds = sister_thresholds,
                                                         pedigree = Sim_pedigree_mat,
                                                         recon_pedigree = NULL,
                                                         simulated = TRUE,
                                                         data_type = "WithGE_2k",
                                                         haplo_assignment_type = "Used mat pedigree")


PatR2_WithGE_SNP50k_Haplo2 <- run_sister_clustering_tests(results_arg = Route2_WithGE_SNP50k$phased_results_flipped,
                                                          sister_thresholds = sister_thresholds,
                                                          pedigree = Sim_pedigree_mat,
                                                          recon_pedigree = NULL,
                                                          simulated = TRUE,
                                                          data_type = "WithGE_50k",
                                                          haplo_assignment_type = "Used mat pedigree")

#Do this for the others too if it works 

# plot_paternity_number_grid_SISTERONLY(
#   Dataset1 = PatR2_SimTrue_Haplo1,
#   Dataset2 = PatR2_NoGE_SNP2k_Haplo1,
#   Dataset3= PatR2_WithGE_SNP2k_Haplo1)
# 
# plot_paternity_number_grid_SISTERONLY(
#   Dataset1 = PatR2_SimTrue_Haplo2,
#   Dataset2 = PatR2_NoGE_SNP2k_Haplo2,
#   Dataset3= PatR2_WithGE_SNP2k_Haplo2)



#** ••• REAL ••• **

colnames(Slov_pedigree_mat_filtered) <- c("id", "dpc", "mother")
colnames(Slov_pedigree_rec_filter) <- c("id", "dpc", "mother")

PatR2_Real_Haplo1 <- run_sister_clustering_tests(results_arg = Route1_Real$results_flipped,
                                                 sister_thresholds = sister_thresholds,
                                                 pedigree = Slov_pedigree_mat_filtered,
                                                 recon_pedigree = Slov_pedigree_rec_filter,
                                                 simulated = FALSE,
                                                 data_type = "Real",
                                                 haplo_assignment_type = "Used recon pedigree")


PatR2_Real_Haplo2 <- run_sister_clustering_tests(results_arg = Route2_Real$phased_results_flipped,
                                                 sister_thresholds = sister_thresholds,
                                                 pedigree = Slov_pedigree_mat_filtered,
                                                 recon_pedigree = NULL,
                                                 simulated = FALSE,
                                                 data_type = "Real",
                                                 haplo_assignment_type = "Used mat pedigree")



save.image(file = paste0(workingDir, "/Data/Pipeline/9_DeterminingPatrilines.Rdata"))



