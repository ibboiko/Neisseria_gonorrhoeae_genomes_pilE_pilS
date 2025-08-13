#version 08: to handle special case when the distance from the pilE to the nearest SmaCla is longer than in the reference genomes, i.e.than 800 bp
# Load needed libraries
library(Biostrings)
library(openxlsx)
library(stringr)

# -------------------------------
# Parameters and target sequences
# -------------------------------
target_sequences <- list(
  dnaA = "ATGACATTAGCAGAGTTTTGGCCGCTGTGCCTCCGCCGTCTTCACGATATGTTGCCTCACGGGCAGTTTGCGCAATGGATTGCGCCCCTTACGGTTGGTGAGGAGGGTGGCGTATGGGTGGTGTACGGCAAGAACCAGTTTGCCTGCAATATGCTCAAGAGCCAGTTTGCCGGAAAAATAGAAGCGGTAAGGGAAGAGTTGGCTGCCGGCCGTCCCGCCTTCGTATTCAAACCGGGAGAAGGCGTGCGTTATGAGATGGCGGCGGTTGAAGGTGCTGTCGAACCTGCCGAGCCGTCCTTGCACGCGGGGTCGGAGGAGATGCCCGTGCAGGAGGTTCTGTTGGACGAGCTGCCGTCTGAAAAGCCTGTCAAACCCGCTGCGTCGAAAACGGCGGCGGATATTTTGGCGGAACGTATGAAAAACCTGCCGCACGAGCCGCGTCAGGCTGCCGGGCCTGCTTCCCGGCCGGAATCGGCGGCAGTTGCCAAAGCGCGGACGGATGCGCAGCGTGATGCGGAAGAAGCGCGTTACGAACAAACCAACCTGTCTCCGGATTACACGTTTGATACGTTGGTAGAAGGTAAGGGCAACCGCCTTGCGGCGGCTGCGGCGCAGGCGATTGCGGAAAACCCGGGGCAGAGTTACAACCCGTTCTTCCTGTACGGCAGCACGGGTTTGGGCAAAACCCACCTTGTGCAGGCGGTCGGCAACGAGCTGTTGAAAAACCGTCCCGATGCCAAAGTGCGCTATATGCATTCGGACGACTACATCCGCAGCTTTATGAAGGCGGTTCGCAACAATACCTACGACGTGTTCAAGCAGCAATACAAGCAATACGACCTGCTGATTATCGACGATATTCAGTTCATCAAAGGAAAAGACCGTACGATGGAAGAATTTTTCTATCTGTACAACCATTTTCACAATGAGAAAAAACAGCTCATCCTCACTTGCGATGTTTTACCCGCCAAAATCGAAGGTATGGACGACCGCCTCAAATCCCGCTTTTCGTGGGGGCTGACTTTGGAACTCGAGCCGCCCGAATTGGAAATGCGTATCGCCATTTTGCAGAAAAAGGCGGAAGCGGCGGGCATCAGTATCGAAGACGAAGCCGCGCTGTTCATTGCCAATCTGATCCGTTCCAACGTGCGCGAACTGGAAGGCGCGTTCAACCGTGTCGGAGCGAGCAGCCGCTTTATGAACCGTCCCGTCATCGACATCGATTTGGCGCGTACCGCTTTGCAGGACATTATTGCCGAGAAGCACAAAGTCATCACCGCCGACATCATCATCGATGCGGTGGCGAAATATTACCGCATCAAAATCAGCGACGTACTCGGCAAAAAACGCACGCGCAACATTGCCCGTCCGCGTCAGGTTGCCATGAGCCTGACCAAAGAATTGACCACTTTGAGCCTGCCGTCTATCGGCGATTCGTTCGGCGGACGCGACCATACGACCGTCATGCACGGCATCAGGGCGGTGGCGAAACTGCGCGAGGAAGACCCCGAGTTGGCGCAGGATTACGAGAAACTGCTGATTCTGATTCAAAACTGA",
  SmaCla = "ATCGATATATTATTTCCACCGGAACGGACGACCCCGCCCGCCTTGCAAACCCTTAAAAGACAAGCCGCCCGGG",
  garP = "TTGAGATTTTTGAATTTACGCGTTAGAAT",
  Leader_sequence_Class_I = "GCCTTTTTGAAGGGTATTCAT"
)

# Default similarity threshold for all targets
similarity_threshold <- 0.8  # 80% similarity threshold

# Specific similarity threshold for SmaCla (updated to 95%)
similarity_threshold_SmaCla <- 0.95  # 95% similarity threshold for SmaCla

# Maximum distance to search for pilE-associated SmaCla repeat from Leader sequence
max_distance_to_leader <- 1000  # 1000 bp flanking region; I increased from 800 to 1000 as we have 2 genomes where this distance is 973 (ID 77942) and 976 (ID 46276)

# Calculate mismatch allowance based on sequence type
get_max_mismatch <- function(target_name, target_seq) {
  target_length <- nchar(target_seq)
  
  if (target_name == "Leader_sequence_Class_I") {
    # Leader sequence: only 100% match
    return(0)
  } else if (target_name %in% c("dnaA", "SmaCla")) {
    # Allow 1 mismatch/gap/insertion per 10 bases
    return(floor(target_length / 10))
  } else {
    # Default: 80% similarity
    return(floor(target_length * (1 - similarity_threshold)))
  }
}

# -----------------------------------
# Similarity function for sequences
# -----------------------------------
calculate_similarity <- function(seq1, seq2) {
  len <- min(nchar(seq1), nchar(seq2))
  if(len == 0) return(0)
  seq1_sub <- substr(seq1, 1, len)
  seq2_sub <- substr(seq2, 1, len)
  matches <- sum(strsplit(seq1_sub, "")[[1]] == strsplit(seq2_sub, "")[[1]])
  return(matches / len)
}

# ---------------------------------------------------------------
# Function to search for target sequence matches in a genomic sequence.
# Uses Biostrings' matchPattern for approximate matching.
# Returns a data frame with match positions, match type, and similarity.
# For SmaCla: checks both forward and reverse complement simultaneously
# For other targets: checks reverse only if forward not found
# ---------------------------------------------------------------
find_target_positions <- function(genomic_seq, target_name, target, max_mismatch) {
  target_length <- nchar(target)
  genomic_dna <- DNAString(genomic_seq)
  target_dna <- DNAString(target)
  
  results <- data.frame(match_type = character(), 
                        Position = numeric(), 
                        similarity = numeric(),
                        trimmed_bp = numeric(),
                        stringsAsFactors = FALSE)
  
  # Special handling for SmaCla: search both orientations simultaneously
  # and use a specific similarity threshold for SmaCla
  if (target_name == "SmaCla") {
    # Use SmaCla-specific similarity threshold
    current_threshold <- similarity_threshold_SmaCla
    
    # Forward search
    f_matches <- matchPattern(target_dna, genomic_dna, max.mismatch = max_mismatch)
    if (length(f_matches) > 0) {
      for (pos in start(f_matches)) {
        if (!is.na(pos) && pos > 0 && (pos + target_length - 1) <= length(genomic_dna)) {
          candidate <- as.character(subseq(genomic_dna, pos, pos + target_length - 1))
          sim <- calculate_similarity(candidate, target)
          
          if (sim >= current_threshold) {
            results <- rbind(results, data.frame(match_type = "forward", 
                                                 Position = pos, 
                                                 similarity = sim,
                                                 trimmed_bp = 0,
                                                 stringsAsFactors = FALSE))
          }
        }
      }
    }
    
    # Reverse complement search (always for SmaCla)
    rc_target_dna <- reverseComplement(target_dna)
    rc_target <- as.character(rc_target_dna)
    
    r_matches <- matchPattern(rc_target_dna, genomic_dna, max.mismatch = max_mismatch)
    if (length(r_matches) > 0) {
      for (pos in start(r_matches)) {
        if (!is.na(pos) && pos > 0 && (pos + target_length - 1) <= length(genomic_dna)) {
          candidate <- as.character(subseq(genomic_dna, pos, pos + target_length - 1))
          # Compare to reverse complement of target for correct similarity calculation
          sim <- calculate_similarity(candidate, rc_target)
          
          if (sim >= current_threshold) {
            # Check if this position already exists in forward results
            # Only add if it's a unique position
            if (!any(results$Position == pos)) {
              results <- rbind(results, data.frame(match_type = "reverse", 
                                                   Position = pos, 
                                                   similarity = sim,
                                                   trimmed_bp = 0,
                                                   stringsAsFactors = FALSE))
            }
          }
        }
      }
    }
  } else {
    # Use default similarity threshold for other targets
    current_threshold <- similarity_threshold
    
    # Standard search for other sequences (unchanged from original)
    # Forward search
    f_matches <- matchPattern(target_dna, genomic_dna, max.mismatch = max_mismatch)
    if (length(f_matches) > 0) {
      for (pos in start(f_matches)) {
        if (!is.na(pos) && pos > 0 && (pos + target_length - 1) <= length(genomic_dna)) {
          candidate <- as.character(subseq(genomic_dna, pos, pos + target_length - 1))
          sim <- calculate_similarity(candidate, target)
          
          if ((target_name == "Leader_sequence_Class_I" && sim == 1) || 
              (target_name != "Leader_sequence_Class_I" && sim >= current_threshold)) {
            results <- rbind(results, data.frame(match_type = "forward", 
                                                 Position = pos, 
                                                 similarity = sim,
                                                 trimmed_bp = 0,
                                                 stringsAsFactors = FALSE))
          }
        }
      }
    }
    
    # Only search reverse if no forward matches with sufficient similarity were found
    if (nrow(results) == 0) {
      # Create reverse complement of target for proper matching
      rc_target_dna <- reverseComplement(target_dna)
      rc_target <- as.character(rc_target_dna)
      
      r_matches <- matchPattern(rc_target_dna, genomic_dna, max.mismatch = max_mismatch)
      if (length(r_matches) > 0) {
        for (pos in start(r_matches)) {
          if (!is.na(pos) && pos > 0 && (pos + target_length - 1) <= length(genomic_dna)) {
            candidate <- as.character(subseq(genomic_dna, pos, pos + target_length - 1))
            # Compare to reverse complement of target for correct similarity calculation
            sim <- calculate_similarity(candidate, rc_target)
            
            if ((target_name == "Leader_sequence_Class_I" && sim == 1) || 
                (target_name != "Leader_sequence_Class_I" && sim >= current_threshold)) {
              results <- rbind(results, data.frame(match_type = "reverse", 
                                                   Position = pos, 
                                                   similarity = sim,
                                                   trimmed_bp = 0,
                                                   stringsAsFactors = FALSE))
            }
          }
        }
      }
    }
  }
  
  return(results)
}

# ---------------------------------------------------------------
# Function to find the pilE-associated SmaCla repeat
# Finds the closest SmaCla repeat within a defined distance of a Leader sequence
# ---------------------------------------------------------------
find_pile_associated_SmaCla <- function(genomic_seq, leader_positions, max_distance) {
  # If no leader sequences found, return empty results
  if (nrow(leader_positions) == 0) {
    return(data.frame(
      match_type = character(),
      Position = numeric(),
      similarity = numeric(),
      distance_to_leader = numeric(),
      leader_position = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  
  # Find all SmaCla repeats
  smacla_target <- target_sequences[["SmaCla"]]
  max_mismatch <- get_max_mismatch("SmaCla", smacla_target)
  smacla_positions <- find_target_positions(genomic_seq, "SmaCla", smacla_target, max_mismatch)
  
  # If no SmaCla repeats found, return empty results
  if (nrow(smacla_positions) == 0) {
    return(data.frame(
      match_type = character(),
      Position = numeric(),
      similarity = numeric(),
      distance_to_leader = numeric(),
      leader_position = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  
  results <- data.frame(
    match_type = character(),
    Position = numeric(),
    similarity = numeric(),
    distance_to_leader = numeric(),
    leader_position = numeric(),
    stringsAsFactors = FALSE
  )
  
  # For each leader sequence, find the closest SmaCla repeat
  for (i in 1:nrow(leader_positions)) {
    leader_pos <- as.numeric(leader_positions$Position[i])
    
    # Calculate distance from each SmaCla repeat to this leader
    distances <- abs(as.numeric(smacla_positions$Position) - leader_pos)
    
    # Find SmaCla repeats within the max distance
    within_distance <- which(distances <= max_distance)
    
    if (length(within_distance) > 0) {
      # If we found SmaCla repeats within the distance, find the closest one
      closest_idx <- within_distance[which.min(distances[within_distance])]
      
      # Add to results
      results <- rbind(results, data.frame(
        match_type = smacla_positions$match_type[closest_idx],
        Position = smacla_positions$Position[closest_idx],
        similarity = smacla_positions$similarity[closest_idx],
        distance_to_leader = distances[closest_idx],
        leader_position = leader_pos,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(results)
}

# Function to create trimmed versions of dnaA sequence (10 to 200 bp trimmed)
create_trimmed_dnaA_sequences <- function(dnaA_seq) {
  trimmed_seqs <- list()
  max_trim <- min(200, nchar(dnaA_seq) - 100)  # Ensure we don't trim too much
  
  for (trim_size in seq(10, max_trim, by = 10)) {
    trimmed_seq <- substr(dnaA_seq, trim_size + 1, nchar(dnaA_seq))
    trimmed_seqs[[paste0("trimmed_", trim_size)]] <- trimmed_seq
  }
  
  return(trimmed_seqs)
}

# ---------------------------------------------------
# Function to extract the ORIGIN sequence from a GenBank file.
# Removes digits, whitespace, and any characters not in the valid nucleotide set.
# ---------------------------------------------------
get_origin_sequence <- function(file) {
  lines <- readLines(file)
  origin_start <- grep("^ORIGIN", lines, ignore.case = TRUE)
  if (length(origin_start) == 0) return("")
  seq_lines <- lines[(origin_start[1] + 1):length(lines)]
  seq_str <- toupper(gsub("[[:digit:][:space:]]", "", paste(seq_lines, collapse = "")))
  seq_str <- gsub("[^ACGTMRWSYKVHDBN]", "", seq_str)
  return(seq_str)
}

# ---------------------------------------------------------------
# Function to calculate relative position to dnaA
# Takes a position, dnaA position, and genome length
# Returns relative position (0-1) and recalculated position with dnaA at position 1
# ---------------------------------------------------------------
calculate_relative_to_dnaA <- function(position, dnaA_position, genome_length) {
  # Skip calculation for "not found" positions
  if (!is.numeric(position)) {
    return(list(
      relative_position = NA, 
      recalculated_position = NA,
      recalculated_relative_position = NA
    ))
  }
  
  # Calculate relative position (0-1)
  relative_position <- position / genome_length
  
  # If dnaA starts at position 1, no recalculation needed
  if (dnaA_position == 1) {
    recalculated_position <- position
  } else {
    # Calculate offset to make dnaA start at position 1
    offset <- dnaA_position - 1
    
    # Feature is after dnaA in the original sequence
    if (position >= dnaA_position) {
      recalculated_position <- position - offset
    } else {
      # Feature is before dnaA in the original sequence
      # Position it after the last recalculated position (they "wrap around")
      recalculated_position <- (genome_length - offset) + position
    }
  }
  
  # Calculate recalculated relative position
  recalculated_relative_position <- recalculated_position / genome_length
  
  return(list(
    relative_position = relative_position, 
    recalculated_position = recalculated_position,
    recalculated_relative_position = recalculated_relative_position
  ))
}

# Function to format genome ID based on special cases
format_genome_id <- function(file_name) {
  # Check for special case: FA1090
  if (grepl("FA1090_N160_CP115654_1_AllFeatures_04MAR2024\\.gb", file_name, ignore.case = TRUE)) {
    return("FA1090")
  }
  # Check for special case: MS11
  else if (grepl("MS11_HL122_CP115904_1_AllFeatures_04MAR2024\\.gb", file_name, ignore.case = TRUE)) {
    return("MS11")
  }
  # Default case: extract the numeric ID
  else {
    return(sub("sequence(\\d+)\\.gb", "\\1", file_name))
  }
}

# ------------------------------------------------------------------------------
# Main analysis function: Process each GenBank (.gb) file, check for each target,
# and output a table with correct positions.
# For each target not found, add a row with "not found".
# ------------------------------------------------------------------------------
analyze_genbank_files <- function(input_folder, output_file, target_sequences, threshold) {
  files <- list.files(input_folder, pattern = "\\.gb$", full.names = TRUE)
  
  results <- data.frame(
    Finished_genome_ID = character(),
    Genome_Length = numeric(),
    Feature = character(),
    Position = character(),  # using character to allow "not found"
    Strand = character(),
    Similarity = character(),
    Targeted_Sequence_Length = numeric(),
    Relative_position_to_dnaA = character(),  # new column
    Recalculated_position = character(),      # new column for debugging/testing
    Recalculated_Relative_position_to_dnaA = character(), # new column
    Notes = character(),
    stringsAsFactors = FALSE
  )
  
  for (file in files) {
    cat("Processing file:", file, "\n")
    genomic_seq <- get_origin_sequence(file)
    if (nchar(genomic_seq) == 0) {
      cat("No ORIGIN sequence found in", file, "\n")
      next
    }
    seq_length <- nchar(genomic_seq)
    
    # Use the new function to format genome ID
    finished_id <- format_genome_id(basename(file))
    
    # Storage for Leader sequences positions in this genome
    leader_matches <- data.frame()
    
    # Storage for file results before calculating relative positions
    file_results <- data.frame(
      Finished_genome_ID = character(),
      Genome_Length = numeric(),
      Feature = character(),
      Position = character(),
      Strand = character(),
      Similarity = character(),
      Targeted_Sequence_Length = numeric(),
      Relative_position_to_dnaA = character(),
      Recalculated_position = character(),
      Recalculated_Relative_position_to_dnaA = character(),
      Notes = character(),
      stringsAsFactors = FALSE
    )
    
    # Process each target
    for (target_name in names(target_sequences)) {
      cat(sprintf("  Searching for %s...\n", target_name))
      cat(sprintf("  Using similarity threshold: %s\n", 
                  ifelse(target_name == "SmaCla", 
                         similarity_threshold_SmaCla, 
                         similarity_threshold)))
      
      cur_target <- target_sequences[[target_name]]
      target_length <- nchar(cur_target)
      max_mismatch <- get_max_mismatch(target_name, cur_target)
      
      matches <- find_target_positions(genomic_seq, target_name, cur_target, max_mismatch)
      
      # Store Leader sequence positions for later use with pilE-associated SmaCla
      if (target_name == "Leader_sequence_Class_I" && nrow(matches) > 0) {
        leader_matches <- matches
      }
      
      # For dnaA, if not found, try with trimmed sequences
      trimmed_notes <- ""
      if (target_name == "dnaA" && nrow(matches) == 0) {
        cat("    Full dnaA not found; searching for trimmed dnaA sequences...\n")
        trimmed_sequences <- create_trimmed_dnaA_sequences(cur_target)
        
        # Try each trimmed sequence
        for (trim_name in names(trimmed_sequences)) {
          trim_size <- as.numeric(gsub("trimmed_", "", trim_name))
          trimmed_seq <- trimmed_sequences[[trim_name]]
          trimmed_length <- nchar(trimmed_seq)
          
          # Recalculate mismatch allowance for trimmed sequence
          trimmed_max_mismatch <- get_max_mismatch("dnaA", trimmed_seq)
          
          trimmed_matches <- find_target_positions(genomic_seq, "dnaA", trimmed_seq, trimmed_max_mismatch)
          
          if (nrow(trimmed_matches) > 0) {
            # Found match with trimmed sequence
            trimmed_notes <- paste("dnaA found with", trim_size, "bp trimmed from start")
            cat("    ", trimmed_notes, "\n")
            
            # Add trimmed_bp information
            trimmed_matches$trimmed_bp <- trim_size
            
            matches <- trimmed_matches
            target_length <- trimmed_length  # Update target length for the output
            break  # Stop once we find a match
          }
        }
      }
      
      if (nrow(matches) == 0) {
        file_results <- rbind(file_results, data.frame(
          Finished_genome_ID = finished_id,
          Genome_Length = seq_length,
          Feature = target_name,
          Position = "not found",
          Strand = "",
          Similarity = "",
          Targeted_Sequence_Length = target_length,
          Relative_position_to_dnaA = "NA",
          Recalculated_position = "NA",
          Recalculated_Relative_position_to_dnaA = "NA",
          Notes = "",
          stringsAsFactors = FALSE
        ))
      } else {
        # Sort by similarity (descending) to ensure we report the best matches first
        matches <- matches[order(matches$similarity, decreasing = TRUE), ]
        
        for (i in 1:nrow(matches)) {
          notes <- ""
          if (matches$trimmed_bp[i] > 0) {
            notes <- paste("dnaA found with", matches$trimmed_bp[i], "bp trimmed from start")
          }
          
          file_results <- rbind(file_results, data.frame(
            Finished_genome_ID = finished_id,
            Genome_Length = seq_length,
            Feature = target_name,
            Position = matches$Position[i],
            Strand = matches$match_type[i],
            Similarity = sprintf("%.3f", matches$similarity[i]),
            Targeted_Sequence_Length = target_length - matches$trimmed_bp[i], # Adjust length for trimmed sequences
            Relative_position_to_dnaA = "NA", # Placeholder, will be calculated later
            Recalculated_position = "NA", # Placeholder, will be calculated later
            Recalculated_Relative_position_to_dnaA = "NA", # Placeholder, will be calculated later
            Notes = notes,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    # After processing all targets, find pilE-associated SmaCla
    cat("  Searching for pilE-associated SmaCla repeats near Leader sequences...\n")
    pile_associated_matches <- find_pile_associated_SmaCla(genomic_seq, leader_matches, max_distance_to_leader)
    
    if (nrow(pile_associated_matches) == 0) {
      file_results <- rbind(file_results, data.frame(
        Finished_genome_ID = finished_id,
        Genome_Length = seq_length,
        Feature = "pilE-associated-SmaCla", # Updated feature name
        Position = "not found",
        Strand = "",
        Similarity = "",
        Targeted_Sequence_Length = nchar(target_sequences[["SmaCla"]]),
        Relative_position_to_dnaA = "NA",
        Recalculated_position = "NA",
        Recalculated_Relative_position_to_dnaA = "NA",
        Notes = "",
        stringsAsFactors = FALSE
      ))
    } else {
      # Sort by distance to leader (ascending)
      pile_associated_matches <- pile_associated_matches[order(pile_associated_matches$distance_to_leader), ]
      
      for (i in 1:nrow(pile_associated_matches)) {
        notes <- sprintf("Distance to Leader: %d bp; Leader position: %d", 
                         pile_associated_matches$distance_to_leader[i],
                         pile_associated_matches$leader_position[i])
        
        file_results <- rbind(file_results, data.frame(
          Finished_genome_ID = finished_id,
          Genome_Length = seq_length,
          Feature = "pilE-associated-SmaCla", # Updated feature name
          Position = pile_associated_matches$Position[i],
          Strand = pile_associated_matches$match_type[i],
          Similarity = sprintf("%.3f", pile_associated_matches$similarity[i]),
          Targeted_Sequence_Length = nchar(target_sequences[["SmaCla"]]),
          Relative_position_to_dnaA = "NA", # Placeholder, will be calculated later
          Recalculated_position = "NA", # Placeholder, will be calculated later
          Recalculated_Relative_position_to_dnaA = "NA", # Placeholder, will be calculated later
          Notes = notes,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Calculate relative positions based on dnaA position
    # First, find the dnaA position from our results
    dnaA_rows <- which(file_results$Feature == "dnaA")
    if (length(dnaA_rows) > 0) {
      # Use the first (best) dnaA match if multiple are found
      dnaA_position <- as.numeric(file_results$Position[dnaA_rows[1]])
      
      # If dnaA was found, calculate relative positions for all features
      for (i in 1:nrow(file_results)) {
        position <- file_results$Position[i]
        
        # Handle numeric positions and "not found" differently
        if (position != "not found") {
          rel_pos <- calculate_relative_to_dnaA(as.numeric(position), dnaA_position, seq_length)
          
          # Format with 3 decimal places
          file_results$Relative_position_to_dnaA[i] <- sprintf("%.3f", rel_pos$relative_position)
          file_results$Recalculated_position[i] <- as.character(rel_pos$recalculated_position)
          file_results$Recalculated_Relative_position_to_dnaA[i] <- sprintf("%.3f", rel_pos$recalculated_relative_position)
        } else {
          file_results$Relative_position_to_dnaA[i] <- "NA"
          file_results$Recalculated_position[i] <- "NA"
          file_results$Recalculated_Relative_position_to_dnaA[i] <- "NA"
        }
      }
    } else {
      # If dnaA was not found, set all relative positions to NA
      file_results$Relative_position_to_dnaA <- "NA"
      file_results$Recalculated_position <- "NA"
      file_results$Recalculated_Relative_position_to_dnaA <- "NA"
    }
    
    # Remove duplicate SmaCla entries that match with pilE-associated-SmaCla
    # Get positions of pilE-associated-SmaCla
    pile_associated_rows <- which(file_results$Feature == "pilE-associated-SmaCla")
    if (length(pile_associated_rows) > 0) {
      # For each pilE-associated-SmaCla, find and remove matching SmaCla entries
      for (row_idx in pile_associated_rows) {
        if (file_results$Position[row_idx] != "not found") {
          pile_pos <- as.numeric(file_results$Position[row_idx])
          pile_strand <- file_results$Strand[row_idx]
          
          # Find SmaCla entries with the same position and strand
          smacla_rows <- which(file_results$Feature == "SmaCla" & 
                                 file_results$Position == pile_pos & 
                                 file_results$Strand == pile_strand)
          
          if (length(smacla_rows) > 0) {
            # Remove the matching SmaCla entries
            file_results <- file_results[-smacla_rows, ]
          }
        }
      }
    }
    
    # Append this file's results to the overall results
    results <- rbind(results, file_results)
  }
  
  write.xlsx(results, output_file, rowNames = FALSE)
  cat("Output written to", output_file, "\n")
  return(results)
}

# -------------------------
# Set input and output paths
# -------------------------
input_dir <- "./Finished_genome_N65"
output_file <- "./v08_SmaCla_Position_Output_N65.xlsx"

# Run the analysis
results <- analyze_genbank_files(input_dir, output_file, target_sequences, similarity_threshold)