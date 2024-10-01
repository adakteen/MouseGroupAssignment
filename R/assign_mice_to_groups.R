#' Assign Mice to Groups with Balanced Tumor Sizes
#'
#' This function assigns mice to groups while ensuring balanced mean tumor sizes, medians, standard errors, and IQRs.
#' It also generates a ggplot showing the distribution of tumor sizes in each group.
#'
#' @param miceData A data frame containing the mice data.
#' @param num_groups The number of groups to assign the mice to.
#' @param tumor_col The name of the tumor size column (as a string).
#' @param cage_col The name of the cage ID column (as a string).
#' @param mouse_col The name of the mouse ID column (as a string).
#' @return A list containing:
#'   \describe{
#'     \item{assignment}{A data frame with the group assignments for each mouse.}
#'     \item{group_stats}{A data frame with statistics for each group.}
#'     \item{outliers}{A data frame with identified outliers.}
#'     \item{plot}{A ggplot object showing the tumor size distribution by group.}
#'   }
#' @details
#' The function first identifies outliers in the tumor size data using the 1.5×IQR rule and removes them.
#' It then assigns mice to groups ensuring that each group has an equivalent number of mice and balanced tumor size metrics.
#' The function attempts to maintain cage integrity by keeping mice from the same cage together when possible.
#'
#' @examples
#' # Load your data
#' miceData <- data.frame(
#'   ID = 1:30,
#'   Cage_Number = rep(1:6, each = 5),
#'   Size = c(2.5, 2.6, 2.7, 2.8, 15.0, 3.0, 3.1, 3.2, 3.1, 3.0,
#'            2.9, 2.8, 2.7, 2.6, 0.5, 3.2, 3.3, 3.4, 3.3, 3.2,
#'            2.7, 2.8, 2.9, 2.8, 2.7, 3.0, 3.1, 3.0, 3.1, 3.0)
#' )
#' tumor_col <- "Size"
#' cage_col <- "Cage_Number"
#' mouse_col <- "ID"
#' num_groups <- 3
#' result <- assign_mice_to_groups(miceData, num_groups, tumor_col, cage_col, mouse_col)
#' print(result$assignment)
#' print(result$group_stats)
#' print(result$plot)
#'
#' @importFrom dplyr %>% filter arrange group_by mutate summarize count full_join left_join ungroup select first n
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter stat_summary labs theme_minimal theme mean_se
#' @importFrom rlang sym enquo .data as_string
#' @export
assign_mice_to_groups <- function(miceData, num_groups, tumor_col, cage_col, mouse_col) {
  # Convert column names to symbols for non-standard evaluation
  tumor_col <- sym(tumor_col)
  cage_col <- sym(cage_col)
  mouse_col <- sym(mouse_col)

  # Identify outliers
  outlier_result <- identify_outliers(miceData, tumor_col)
  data_no_outliers <- outlier_result$data_no_outliers
  outliers <- outlier_result$outliers

  # Total number of mice after removing outliers
  total_mice <- nrow(data_no_outliers)

  # Calculate group sizes
  base_group_size <- floor(total_mice / num_groups)
  extra_mice <- total_mice %% num_groups  # Remaining mice to distribute
  group_sizes <- rep(base_group_size, num_groups)
  if (extra_mice > 0) {
    group_sizes[1:extra_mice] <- group_sizes[1:extra_mice] + 1
  }

  # Sort mice by tumor size
  data_no_outliers <- data_no_outliers %>% arrange(!!tumor_col)

  # Determine the number of strata (equal to group size)
  num_strata <- base_group_size
  if (extra_mice > 0) {
    num_strata <- num_strata + 1
  }

  # Assign strata labels to mice
  data_no_outliers$Stratum <- rep(1:num_strata, each = num_groups, length.out = total_mice)

  # Initialize group assignments
  data_no_outliers$Group <- NA
  data_no_outliers$MovedCage <- FALSE

  # Assign mice to groups by rotating group labels within each stratum
  for (stratum in unique(data_no_outliers$Stratum)) {
    stratum_mice <- data_no_outliers %>% filter(Stratum == stratum)
    groups_in_stratum <- sample(1:num_groups, nrow(stratum_mice), replace = FALSE)
    stratum_mice$Group <- groups_in_stratum
    data_no_outliers$Group[data_no_outliers[[as_string(mouse_col)]] %in% stratum_mice[[as_string(mouse_col)]]] <- stratum_mice$Group
  }

  # Ensure group sizes are correct
  group_counts <- data_no_outliers %>% count(Group)
  while (any(group_counts$n != group_sizes)) {
    # Find groups that are over or under the desired size
    over_groups <- group_counts$Group[group_counts$n > group_sizes]
    under_groups <- group_counts$Group[group_counts$n < group_sizes]

    # Move mice from over_groups to under_groups
    for (group_over in over_groups) {
      mice_in_over_group <- data_no_outliers %>% filter(Group == group_over)
      mice_to_move <- mice_in_over_group[1, ]  # Select one mouse to move
      group_under <- under_groups[1]
      data_no_outliers$Group[data_no_outliers[[as_string(mouse_col)]] == mice_to_move[[as_string(mouse_col)]]] <- group_under
      data_no_outliers$MovedCage[data_no_outliers[[as_string(mouse_col)]] == mice_to_move[[as_string(mouse_col)]]] <- TRUE

      # Update group counts
      group_counts <- data_no_outliers %>% count(Group)
      # Update over and under groups
      over_groups <- group_counts$Group[group_counts$n > group_sizes]
      under_groups <- group_counts$Group[group_counts$n < group_sizes]
      if (length(over_groups) == 0 || length(under_groups) == 0) {
        break
      }
    }
  }

  # Attempt to keep cage integrity
  # Identify mice that need to be moved to keep cages together
  data_no_outliers <- data_no_outliers %>%
    group_by(!!cage_col) %>%
    mutate(CageGroup = first(Group)) %>%
    ungroup()

  # If mice from the same cage are in different groups, try to adjust
  cages <- unique(data_no_outliers[[as_string(cage_col)]])
  for (cage in cages) {
    cage_data <- data_no_outliers %>% filter((!!cage_col) == cage)
    assigned_groups <- unique(cage_data$Group)
    if (length(assigned_groups) > 1) {
      # Attempt to move mice to the majority group for that cage, provided group sizes allow
      target_group <- assigned_groups[which.max(table(cage_data$Group))]
      other_groups <- setdiff(assigned_groups, target_group)
      for (grp in other_groups) {
        mice_to_move <- cage_data %>% filter(Group == grp)
        for (i in 1:nrow(mice_to_move)) {
          # Check if moving the mouse will not exceed group size
          if (sum(data_no_outliers$Group == target_group) < group_sizes[target_group]) {
            data_no_outliers$Group[data_no_outliers[[as_string(mouse_col)]] == mice_to_move[[as_string(mouse_col)]][i]] <- target_group
            data_no_outliers$MovedCage[data_no_outliers[[as_string(mouse_col)]] == mice_to_move[[as_string(mouse_col)]][i]] <- FALSE
          }
        }
      }
    }
  }

  # Recalculate group statistics
  group_stats <- data_no_outliers %>%
    group_by(Group) %>%
    summarize(
      GroupSize = n(),
      GroupMean = mean(!!tumor_col),
      GroupMedian = median(!!tumor_col),
      GroupSE = sd(!!tumor_col) / sqrt(n()),
      GroupIQR = IQR(!!tumor_col)
    )

  # Generate ggplot
  plot <- ggplot(data_no_outliers, aes(x = factor(Group), y = !!tumor_col)) +
    geom_boxplot(outlier.shape = NA, fill = "lightblue", alpha = 0.5) +
    geom_jitter(width = 0.2, size = 2, aes(color = factor(Group))) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, fill = "red") +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
    labs(title = "Tumor Sizes by Group",
         x = "Group",
         y = "Tumor Size") +
    theme_minimal() +
    theme(legend.position = "none")

  # Output results
  assignment <- data_no_outliers %>% select(!!mouse_col, !!cage_col, !!tumor_col, Group, MovedCage)
  # Rename columns to original names for clarity
  colnames(assignment)[1:3] <- c("MouseID", "CageID", "TumorSize")
  group_stats <- group_stats
  outliers <- outliers
  return(list(assignment = assignment, group_stats = group_stats, outliers = outliers, plot = plot))
}

#' Identify Outliers Using the 1.5*IQR Rule
#'
#' This function identifies outliers in the tumor size data.
#'
#' @param data A data frame containing the data.
#' @param tumor_col The name of the tumor size column (as a symbol).
#' @return A list containing:
#'   \item{data_no_outliers}{Data frame without outliers.}
#'   \item{outliers}{Data frame with outliers.}
#' @importFrom dplyr %>% filter
#' @importFrom rlang sym
#' @export
identify_outliers <- function(data, tumor_col) {
  # If tumor_col is a string, convert it to a symbol
  if (is.character(tumor_col)) {
    tumor_col <- sym(tumor_col)
  }
  Q1 <- quantile(data[[as_string(tumor_col)]], 0.25)
  Q3 <- quantile(data[[as_string(tumor_col)]], 0.75)
  IQR_value <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR_value
  upper_bound <- Q3 + 1.5 * IQR_value
  outliers <- data %>% filter((!!tumor_col < lower_bound) | (!!tumor_col > upper_bound))
  data_no_outliers <- data %>% filter((!!tumor_col >= lower_bound) & (!!tumor_col <= upper_bound))
  return(list(data_no_outliers = data_no_outliers, outliers = outliers))
}
