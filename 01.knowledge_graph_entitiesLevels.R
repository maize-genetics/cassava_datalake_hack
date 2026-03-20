# Knowledge Graph with entity parent-child relationships and biological metadata
# Bethany Econopouly, econopouly@cornell.edu
# With help from Claude-4
# In progress script (test)

#Load libraries
library(yaml)
library(visNetwork)
library(htmlwidgets)
library(dplyr)

# Function to calculate hierarchy levels
calculate_hierarchy_levels <- function(processes, all_process_ids) {
  
  # Build parent-child relationships
  children_map <- list()
  parents_map <- list()
  
  for(process in processes) {
    process_id <- process$process_id
    if(!is.null(process$child_processes) && process$child_processes != "") {
      children <- strsplit(as.character(process$child_processes), ",")[[1]]
      children <- trimws(children)
      children <- children[children %in% all_process_ids]
      
      children_map[[process_id]] <- children
      
      for(child in children) {
        parents_map[[child]] <- c(parents_map[[child]], process_id)
      }
    }
  }
  
  # Find root nodes
  root_nodes <- setdiff(all_process_ids, names(parents_map))
  
  # Calculate levels using BFS
  levels <- list()
  queue <- data.frame(id = root_nodes, level = 1, stringsAsFactors = FALSE)
  
  while(nrow(queue) > 0) {
    current <- queue[1, ]
    queue <- queue[-1, , drop = FALSE]
    
    if(is.null(levels[[current$id]])) {
      levels[[current$id]] <- current$level
      
      if(!is.null(children_map[[current$id]])) {
        for(child in children_map[[current$id]]) {
          if(is.null(levels[[child]])) {
            queue <- rbind(queue, data.frame(id = child, level = current$level + 1, stringsAsFactors = FALSE))
          }
        }
      }
    }
  }
  
  return(levels)
}

create_knowledge_graph <- function(yaml_file) {
  
  yaml_data <- read_yaml(yaml_file)
  
  # Create ALL nodes and edges with biological metadata columns
  nodes <- data.frame(
    id = character(),
    label = character(),
    title = character(),
    group = character(),
    size = numeric(),
    color = character(),
    node_type = character(),
    process_type = character(),
    supplied_accession_name = character(),
    sample_label = character(),
    ncbi_biosample = character(),
    origin = character(),
    provider = character(),
    sequencing_method = character(),
    fastq_location = character(),
    accession_type = character(),
    entity_description = character(),
    stringsAsFactors = FALSE
  )
  
  edges <- data.frame(
    from = character(),
    to = character(),
    arrows = character(),
    color = character(),
    stringsAsFactors = FALSE
  )
  
  all_process_ids <- c()
  all_process_types <- c()
  all_entity_ids <- c()
  used_ids <- c()
  
  # Store pending parent-child entity relationships
  pending_parent_edges <- list()
  
  # First pass: collect process types and IDs
  for(i in 1:length(yaml_data$processes)) {
    process <- yaml_data$processes[[i]]
    process_type <- process$type %||% "unknown"
    all_process_types <- c(all_process_types, process_type)
    
    original_id <- process$process_id %||% paste0("proc_", i)
    process_id <- original_id
    
    counter <- 1
    while(process_id %in% used_ids) {
      process_id <- paste0(original_id, "_dup_", counter)
      counter <- counter + 1
    }
    
    used_ids <- c(used_ids, process_id)
    all_process_ids <- c(all_process_ids, process_id)
  }
  
  # Add "entity" as a type
  all_process_types <- c(all_process_types, "entity")
  
  # Create color mapping
  unique_types <- unique(all_process_types)
  color_palette <- c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#E17055", 
                     "#FFEAA7", "#DDA0DD", "#F8BBD9", "#B8860B", "#F7DC6F",
                     "#FF8C94", "#98D8C8", "#F06292", "#AED581", "#FFB74D",
                     "#A1C4FD", "#C2E9FB", "#FCCB90", "#D1A3FF", "#FFE0B2")
  
  type_colors <- setNames(color_palette[1:length(unique_types)], unique_types)
  
  # Calculate hierarchy levels
  hierarchy_levels <- calculate_hierarchy_levels(yaml_data$processes, all_process_ids)
  
  # Second pass: create ALL nodes (processes AND entities)
  used_ids <- c()
  
  for(i in 1:length(yaml_data$processes)) {
    process <- yaml_data$processes[[i]]
    
    original_id <- process$process_id %||% paste0("proc_", i)
    process_id <- original_id
    
    counter <- 1
    while(process_id %in% used_ids) {
      process_id <- paste0(original_id, "_dup_", counter)
      counter <- counter + 1
    }
    used_ids <- c(used_ids, process_id)
    
    process_name <- process$name %||% "Unnamed Process"
    process_type <- process$type %||% "unknown"
    
    hierarchy_level <- hierarchy_levels[[process_id]] %||% 1
    node_size <- max(15, 50 - (hierarchy_level * 8))
    
    entity_count <- length(process$entities %||% list())
    
    # Add process node with biological metadata
    nodes <- rbind(nodes, data.frame(
      id = process_id,
      label = paste0(process_id, "\n", substr(process_name, 1, 25)),
      title = paste0("<b>", process_name, "</b><br>",
                     "Type: ", process_type, "<br>",
                     "ID: ", process_id, "<br>",
                     "Hierarchy Level: ", hierarchy_level, "<br>",
                     "Entities: ", entity_count, "<br>",
                     "Description: ", substr(process$description %||% "No description", 1, 150), "..."),
      group = process_type,
      size = node_size,
      color = type_colors[process_type],
      # Add process metadata columns
      node_type = "process",
      process_type = process_type,
      supplied_accession_name = NA,
      sample_label = NA,
      ncbi_biosample = NA,
      origin = NA,
      provider = NA,
      sequencing_method = NA,
      fastq_location = NA,
      accession_type = NA,
      entity_description = NA,
      stringsAsFactors = FALSE
    ))
    
    # Create edges to child processes
    if(!is.null(process$child_processes) && process$child_processes != "") {
      children <- strsplit(as.character(process$child_processes), ",")[[1]]
      children <- trimws(children)
      
      for(child in children) {
        if(child != "" && !is.na(child)) {
          edges <- rbind(edges, data.frame(
            from = process_id,
            to = child,
            arrows = "to",
            color = "#666666",
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    # Create individual entity nodes with biological metadata
    if(!is.null(process$entities) && length(process$entities) > 0) {
      for(j in 1:length(process$entities)) {
        entity <- process$entities[[j]]
        entity_id <- entity$entity_iid %||% paste0(process_id, ".entity.", j)
        
        counter <- 1
        original_entity_id <- entity_id
        while(entity_id %in% all_entity_ids) {
          entity_id <- paste0(original_entity_id, "_dup_", counter)
          counter <- counter + 1
        }
        all_entity_ids <- c(all_entity_ids, entity_id)
        
        # Extract biological metadata
        supplied_accession_name <- entity$metadata$organism$supplied_accession_name %||% 
          entity$metadata$sample$supplied_accession_name %||% NA
        sample_label <- entity$metadata$sample$label %||% 
          entity$metadata$assay$label %||% NA
        ncbi_biosample <- entity$metadata$sample$ncbi_biosample %||% 
          entity$metadata$biosample$ncbi_biosample %||% NA
        origin <- entity$metadata$sample$origin %||% 
          entity$metadata$organism$origin_notes %||% NA
        provider <- entity$metadata$sample$provider %||% NA
        sequencing_method <- entity$metadata$assay$sequencing_method %||% NA
        fastq_location <- entity$metadata$assay$fastq_location %||% NA
        accession_type <- entity$metadata$organism$accession_type %||% NA
        description <- entity$metadata$organism$description %||% NA
        
        entity_label <- supplied_accession_name %||% sample_label %||% entity_id
        
        # Create enhanced entity title
        entity_title <- paste0("<b>Entity: ", entity_id, "</b><br>")
        if(!is.na(supplied_accession_name)) {
          entity_title <- paste0(entity_title, "Accession: ", supplied_accession_name, "<br>")
        }
        if(!is.na(sample_label)) {
          entity_title <- paste0(entity_title, "Label: ", sample_label, "<br>")
        }
        if(!is.na(sequencing_method)) {
          entity_title <- paste0(entity_title, "Method: ", sequencing_method, "<br>")
        }
        if(!is.na(origin)) {
          entity_title <- paste0(entity_title, "Origin: ", origin, "<br>")
        }
        
        # Create enhanced nodes dataframe with biological metadata
        nodes <- rbind(nodes, data.frame(
          id = entity_id,
          label = substr(entity_label, 1, 15),
          title = entity_title,
          group = paste0("entity_", process_id),
          size = 8,
          color = type_colors["entity"],
          # Add biological metadata columns
          node_type = "entity",
          process_type = NA,
          supplied_accession_name = supplied_accession_name %||% NA,
          sample_label = sample_label %||% NA,
          ncbi_biosample = ncbi_biosample %||% NA,
          origin = origin %||% NA,
          provider = provider %||% NA,
          sequencing_method = sequencing_method %||% NA,
          fastq_location = fastq_location %||% NA,
          accession_type = accession_type %||% NA,
          entity_description = description %||% NA,
          stringsAsFactors = FALSE
        ))
        
        # Edge from process to entity
        edges <- rbind(edges, data.frame(
          from = process_id,
          to = entity_id,
          arrows = "to",
          color = "#CCCCCC",
          stringsAsFactors = FALSE
        ))
        
        # Store parent relationships for later processing
        if(!is.null(entity$parents)) {
          for(parent in entity$parents) {
            parent_id <- parent$entity_iid
            if(!is.null(parent_id) && parent_id != "") {
              pending_parent_edges <- c(pending_parent_edges, 
                                        list(list(from = parent_id, to = entity_id)))
            }
          }
        }
      }
    }
  }
  
  # Now add parent-child entity edges
  cat("Processing", length(pending_parent_edges), "potential entity parent-child relationships...\n")
  entity_parent_edges_added <- 0
  
  for(parent_edge in pending_parent_edges) {
    # Only add edge if both parent and child entities exist
    if(parent_edge$from %in% all_entity_ids && parent_edge$to %in% all_entity_ids) {
      edges <- rbind(edges, data.frame(
        from = parent_edge$from,
        to = parent_edge$to,
        arrows = "to",
        color = "#FFB6C1",  # Light pink for entity-entity relationships
        stringsAsFactors = FALSE
      ))
      entity_parent_edges_added <- entity_parent_edges_added + 1
    }
  }
  
  cat("Added", entity_parent_edges_added, "entity parent-child relationships\n")
  
  # Check for duplicates
  if(any(duplicated(nodes$id))) {
    duplicated_ids <- nodes$id[duplicated(nodes$id)]
    cat("ERROR: Duplicate IDs:", paste(duplicated_ids, collapse = ", "), "\n")
    return(NULL)
  }
  
  # Remove edges to non-existent nodes
  missing_targets <- setdiff(edges$to, nodes$id)
  if(length(missing_targets) > 0) {
    cat("Warning: Removing edges to missing nodes:", paste(missing_targets, collapse = ", "), "\n")
    edges <- edges[edges$to %in% nodes$id, ]
  }
  
  cat("Created", nrow(nodes), "nodes and", nrow(edges), "edges\n")
  
  # Create network
  network <- visNetwork(nodes, edges, width = "100%", height = "100vh") %>%
    visNodes(font = list(size = 10)) %>%
    visEdges(arrows = list(to = list(enabled = TRUE, scaleFactor = 0.8)),
             smooth = list(enabled = TRUE, type = "cubicBezier")) %>%
    visOptions(highlightNearest = list(enabled = TRUE, hover = TRUE, degree = 1),
               nodesIdSelection = TRUE) %>%
    visLayout(randomSeed = 123) %>%
    visInteraction(navigationButtons = TRUE, 
                   dragNodes = TRUE, 
                   dragView = TRUE, 
                   zoomView = TRUE) %>%
    visPhysics(enabled = TRUE, 
               stabilization = list(iterations = 100),
               adaptiveTimestep = TRUE,
               timestep = 0.5,
               maxVelocity = 30) %>%
    # Cluster entities by their parent process
    visClusteringByGroup(groups = grep("^entity_", unique(nodes$group), value = TRUE)) %>%
    visOptions(collapse = list(enabled = TRUE,
                               fit = TRUE,
                               resetHighlight = TRUE,
                               clusterOptions = list(borderWidth = 3,
                                                     childBorderWidth = 2,
                                                     font = list(size = 14))))
  
  # Store nodes and edges data for lookup and saving
  attr(network, "nodes_data") <- nodes
  attr(network, "edges_data") <- edges
  attr(network, "type_colors") <- type_colors
  return(network)
}

# Simple legend function
add_colored_legend <- function(widget_file, type_colors) {
  
  html_content <- readLines(widget_file)
  head_end <- grep("</head>", html_content)
  
  legend_css <- c(
    "<style>",
    "body { margin: 0; padding: 0; font-family: Arial, sans-serif; }",
    ".html-widget { width: 100vw !important; height: 100vh !important; }",
    ".legend-container {",
    "  position: absolute;",
    "  top: 60px;",
    "  right: 20px;",
    "  background: rgba(255, 255, 255, 0.95);",
    "  border: 2px solid #ddd;",
    "  border-radius: 8px;",
    "  padding: 15px;",
    "  font-size: 12px;",
    "  z-index: 1000;",
    "  box-shadow: 0 4px 8px rgba(0,0,0,0.1);",
    "  max-width: 200px;",
    "  max-height: 70vh;",
    "  overflow-y: auto;",
    "}",
    ".legend-title {",
    "  font-weight: bold;",
    "  margin-bottom: 10px;",
    "  color: #333;",
    "  border-bottom: 2px solid #eee;",
    "  padding-bottom: 5px;",
    "  font-size: 14px;",
    "}",
    ".legend-item {",
    "  display: flex;",
    "  align-items: center;",
    "  margin: 6px 0;",
    "  color: #555;",
    "  padding: 2px 0;",
    "}",
    ".legend-color {",
    "  width: 12px;",
    "  height: 12px;",
    "  border-radius: 50%;",
    "  margin-right: 8px;",
    "  border: 1px solid #999;",
    "  flex-shrink: 0;",
    "}",
    ".legend-text {",
    "  font-size: 11px;",
    "  line-height: 1.2;",
    "  word-wrap: break-word;",
    "}",
    "</style>"
  )
  
  legend_items <- c()
  for(type_name in names(type_colors)) {
    color <- type_colors[type_name]
    legend_items <- c(legend_items, 
                      paste0('  <div class="legend-item">',
                             '<div class="legend-color" style="background-color: ', color, ';"></div>',
                             '<div class="legend-text">', type_name, '</div>',
                             '</div>'))
  }
  
  legend_html <- c(
    '<div class="legend-container">',
    '  <div class="legend-title">Process Types</div>',
    legend_items,
    '  <div style="margin-top: 10px; font-size: 10px; color: #666;">',
    '    🔗 Light pink arrows = Entity parent-child relationships<br>',
    '    🔗 Gray arrows = Process-entity relationships<br>',
    '    🔗 Dark gray arrows = Process parent-child relationships',
    '  </div>',
    '</div>'
  )
  
  html_content <- c(html_content[1:(head_end-1)], 
                    legend_css,
                    html_content[head_end:length(html_content)])
  
  body_start <- grep("<body", html_content)
  html_content <- c(html_content[1:body_start], 
                    legend_html,
                    html_content[(body_start+1):length(html_content)])
  
  writeLines(html_content, widget_file)
}

# Create and save the graph
knowledge_graph <- create_knowledge_graph("input/catalogue15.yaml")

if(!is.null(knowledge_graph)) {
  type_colors <- attr(knowledge_graph, "type_colors")
  
  saveWidget(
    knowledge_graph,
    "cassava_knowledge_graph_with_entity_parents.html",
    selfcontained = TRUE,
    title = "Cassava Data Lake Knowledge Graph - With Entity Parents"
  )
  
  add_colored_legend("cassava_knowledge_graph_with_entity_parents.html", type_colors)
  
  cat("Graph saved with entity parent-child relationships!\n")
  knowledge_graph
} else {
  cat("Graph creation failed\n")
}

## Save data 
nodes_df <- attr(knowledge_graph, "nodes_data") 
edges_df <- attr(knowledge_graph, "edges_data")
write.csv(nodes_df, "cassava_nodes.csv", row.names = FALSE)
write.csv(edges_df, "cassava_edges.csv", row.names = FALSE)

## Look at data to see if makes sense and has relevent fields
# Extract the data
nodes_df <- attr(knowledge_graph, "nodes_data") 
edges_df <- attr(knowledge_graph, "edges_data")

# Test 1: Check if biological metadata columns exist
cat("=== COLUMN CHECK ===\n")
cat("Nodes columns:", paste(names(nodes_df), collapse = ", "), "\n\n")

# Look at accession names
cat("=== ENTITIES WITH ACCESSION NAMES ===\n")
entities_with_accessions <- nodes_df[nodes_df$node_type == "entity" & !is.na(nodes_df$supplied_accession_name), ]
if(nrow(entities_with_accessions) > 0) {
  for(i in 1:min(5, nrow(entities_with_accessions))) {
    cat("Entity:", entities_with_accessions$id[i], 
        "-> Accession:", entities_with_accessions$supplied_accession_name[i], "\n")
  }
  cat("Found", nrow(entities_with_accessions), "entities with accession names\n\n")
} else {
  cat("No entities with accession names found\n\n")
}

# Check for the specific examples 
cat("=== SPECIFIC TEST CASES ===\n")
test_entities <- c("pr00178.15", "pr00134.51", "pr00173.67")
for(entity_id in test_entities) {
  entity_row <- nodes_df[nodes_df$id == entity_id, ]
  if(nrow(entity_row) > 0) {
    cat("Found", entity_id, ":\n")
    cat("  - Accession:", entity_row$supplied_accession_name, "\n")
    cat("  - Sample Label:", entity_row$sample_label, "\n")
    cat("  - Origin:", entity_row$origin, "\n")
    cat("  - Provider:", entity_row$provider, "\n\n")
  } else {
    cat("Entity", entity_id, "not found\n\n")
  }
}

# Check processes have types
cat("=== PROCESS TYPES ===\n")
process_types <- table(nodes_df$process_type[nodes_df$node_type == "process"])
cat("Process type counts:\n")
print(process_types)
cat("\n")

# Check for entity parent-child relationships
cat("=== ENTITY PARENT-CHILD RELATIONSHIPS ===\n")
pink_edges <- edges_df[edges_df$color == "#FFB6C1", ]
if(nrow(pink_edges) > 0) {
  cat("Found", nrow(pink_edges), "entity parent-child relationships\n")
  for(i in 1:min(3, nrow(pink_edges))) {
    parent_name <- nodes_df$supplied_accession_name[nodes_df$id == pink_edges$from[i]]
    child_name <- nodes_df$supplied_accession_name[nodes_df$id == pink_edges$to[i]]
    cat("  ", pink_edges$from[i], "(", parent_name, ") -> ", 
        pink_edges$to[i], "(", child_name, ")\n")
  }
} else {
  cat("No entity parent-child relationships found\n")
}

# Example, for a child entity get accession from parent entity
get_accession_name <- function(entity_id, nodes_df, edges_df) {
  # First check if entity has its own accession name
  entity_row <- nodes_df[nodes_df$id == entity_id, ]
  if(nrow(entity_row) == 0) {
    return("Entity not found")
  }
  
  if(!is.na(entity_row$supplied_accession_name)) {
    return(entity_row$supplied_accession_name)
  }
  
  # If no accession name, look for parent with accession name
  parent_edges <- edges_df[edges_df$to == entity_id & edges_df$color == "#FFB6C1", ]
  
  if(nrow(parent_edges) > 0) {
    for(i in 1:nrow(parent_edges)) {
      parent_id <- parent_edges$from[i]
      parent_row <- nodes_df[nodes_df$id == parent_id, ]
      if(!is.na(parent_row$supplied_accession_name)) {
        return(paste0(parent_row$supplied_accession_name, " (from parent ", parent_id, ")"))
      }
    }
  }
  
  return("No accession name found")
}

## Test accession name retrieval 
get_accession_name("pr00199.6", nodes_df, edges_df)


# Function to find entities belonging to processes of a specific type
find_entities_by_process_type <- function(nodes_df, edges_df, target_process_type) {
  
  # Find all processes of the target type
  target_processes <- nodes_df[nodes_df$node_type == "process" & 
                                 grepl(target_process_type, nodes_df$process_type, ignore.case = TRUE), ]
  
  cat("Found", nrow(target_processes), "processes of type containing '", target_process_type, "'\n")
  
  if(nrow(target_processes) == 0) {
    return(data.frame())
  }
  
  # Print the process types found
  cat("Process types found:\n")
  print(unique(target_processes$process_type))
  
  # Find all entities that belong to these processes
  target_entities <- data.frame()
  
  for(i in 1:nrow(target_processes)) {
    process_id <- target_processes$id[i]
    
    # Find entities that have edges FROM this process
    process_entities <- edges_df[edges_df$from == process_id & 
                                   edges_df$color == "#CCCCCC", ]  # Process-to-entity edges
    
    if(nrow(process_entities) > 0) {
      entity_ids <- process_entities$to
      entities <- nodes_df[nodes_df$id %in% entity_ids & nodes_df$node_type == "entity", ]
      
      if(nrow(entities) > 0) {
        entities$parent_process_id <- process_id
        entities$parent_process_type <- target_processes$process_type[i]
        target_entities <- rbind(target_entities, entities)
      }
    }
  }
  
  return(target_entities)
}

# Find WGS assay entities and info
find_parent_study_and_investigation <- function(entity_id, nodes_df, edges_df, max_depth = 10) {
  
  study_id <- NA
  investigation_id <- NA
  
  # Start from the entity and work up the hierarchy
  current_nodes <- c(entity_id)
  visited <- c()
  
  for(depth in 1:max_depth) {
    if(length(current_nodes) == 0) break
    
    next_nodes <- c()
    
    for(node_id in current_nodes) {
      if(node_id %in% visited) next
      visited <- c(visited, node_id)
      
      # Check if this node is a study or investigation
      node_info <- nodes_df[nodes_df$id == node_id, ]
      if(nrow(node_info) > 0) {
        node_type <- node_info$process_type[1]
        if(!is.na(node_type)) {
          if(grepl("study", node_type, ignore.case = TRUE) && is.na(study_id)) {
            study_id <- node_id
          }
          if(grepl("investigation", node_type, ignore.case = TRUE) && is.na(investigation_id)) {
            investigation_id <- node_id
          }
        }
      }
      
      # Find parent nodes (nodes that point TO this node)
      parent_edges <- edges_df[edges_df$to == node_id, ]
      if(nrow(parent_edges) > 0) {
        parents <- parent_edges$from
        next_nodes <- c(next_nodes, parents)
      }
    }
    
    current_nodes <- unique(next_nodes)
    
    # Stop if we found both study and investigation
    if(!is.na(study_id) && !is.na(investigation_id)) {
      break
    }
  }
  
  return(list(
    study_id = if(is.na(study_id)) "Not found" else study_id,
    investigation_id = if(is.na(investigation_id)) "Not found" else investigation_id
  ))
}

# Enhanced function to get study and investigation names
get_process_name <- function(process_id, nodes_df) {
  if(is.na(process_id) || process_id == "Not found") {
    return("Not found")
  }
  
  process_info <- nodes_df[nodes_df$id == process_id, ]
  if(nrow(process_info) > 0) {
    # Extract name from title (it's in the format "process_id\nname")
    title <- process_info$title[1]
    if(!is.na(title)) {
      # Parse the title to extract the process name
      title_parts <- strsplit(title, "<br>")[[1]]
      name_line <- title_parts[1]
      name_line <- gsub("<b>|</b>", "", name_line)  # Remove bold tags
      return(name_line)
    }
  }
  return("Unknown")
}

# Re-run the WGS entity extraction with enhanced details
cat("=== WGS ENTITY EXTRACTION WITH STUDY/INVESTIGATION INFO ===\n")

# Find WGS assay entities (reusing previous function)
wgs_entities <- find_entities_by_process_type(nodes_df, edges_df, "assay.*wgs")

cat("Found", nrow(wgs_entities), "entities from WGS assay processes\n")

if(nrow(wgs_entities) > 0) {
  # Enhanced details with study and investigation info
  enhanced_details <- c()
  accession_names <- c()
  
  cat("Processing entities and finding parent studies/investigations...\n")
  
  for(i in 1:nrow(wgs_entities)) {
    if(i %% 100 == 0) cat("Processed", i, "of", nrow(wgs_entities), "entities...\n")
    
    entity_id <- wgs_entities$id[i]
    
    # Get accession name
    accession <- get_accession_name(entity_id, nodes_df, edges_df)
    clean_accession <- gsub(" \\(from parent.*\\)", "", as.character(accession))
    
    # Only process entities with valid accession names
    if(length(clean_accession) > 0 && 
       clean_accession != "No accession name found" && 
       clean_accession != "Entity not found" &&
       clean_accession != "NA" &&
       !is.na(clean_accession) &&
       nchar(clean_accession) > 0) {
      
      accession_names <- c(accession_names, clean_accession)
      
      # Find parent study and investigation
      parents <- find_parent_study_and_investigation(entity_id, nodes_df, edges_df)
      
      # Get readable names for study and investigation
      study_name <- get_process_name(parents$study_id, nodes_df)
      investigation_name <- get_process_name(parents$investigation_id, nodes_df)
      
      # Create enhanced detail row
      enhanced_details <- c(enhanced_details, 
                            paste0(entity_id, "\t", 
                                   clean_accession, "\t",
                                   wgs_entities$parent_process_id[i], "\t",
                                   wgs_entities$parent_process_type[i], "\t",
                                   parents$study_id, "\t",
                                   study_name, "\t",
                                   parents$investigation_id, "\t",
                                   investigation_name))
    }
  }
  
  # Create enhanced output files
  unique_accessions <- sort(unique(as.character(accession_names)))
  unique_accessions <- unique_accessions[!is.na(unique_accessions) & nchar(unique_accessions) > 0]
  
  cat("Found", length(unique_accessions), "unique accession names from WGS assays\n")
  
  if(length(enhanced_details) > 0) {
    # Write enhanced detailed file
    enhanced_header <- "Entity_ID\tAccession_Name\tWGS_Process_ID\tWGS_Process_Type\tStudy_ID\tStudy_Name\tInvestigation_ID\tInvestigation_Name"
    writeLines(c(enhanced_header, enhanced_details), "wgs_assay_enhanced_detailed.txt")
    
    # Also update the simple accession names file
    writeLines(unique_accessions, "wgs_assay_accession_names.txt")
    
    cat("Created enhanced files:\n")
    cat("- wgs_assay_accession_names.txt (just accession names)\n") 
    cat("- wgs_assay_enhanced_detailed.txt (with study and investigation info)\n")
    
    # Show a preview of the enhanced data
    cat("\nPreview of enhanced data (first 3 entries):\n")
    cat("Entity_ID\tAccession_Name\tWGS_Process_ID\tWGS_Process_Type\tStudy_ID\tStudy_Name\tInvestigation_ID\tInvestigation_Name\n")
    for(i in 1:min(3, length(enhanced_details))) {
      cat(enhanced_details[i], "\n")
    }
    
    # Summary statistics
    enhanced_df <- read.table(text = paste(c(enhanced_header, enhanced_details), collapse = "\n"), 
                              header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
    
    cat("\n=== SUMMARY STATISTICS ===\n")
    cat("Total entities:", nrow(enhanced_df), "\n")
    cat("Unique studies found:", length(unique(enhanced_df$Study_ID[enhanced_df$Study_ID != "Not found"])), "\n")
    cat("Unique investigations found:", length(unique(enhanced_df$Investigation_ID[enhanced_df$Investigation_ID != "Not found"])), "\n")
    
    cat("\nStudies with WGS data:\n")
    study_counts <- table(enhanced_df$Study_Name[enhanced_df$Study_Name != "Not found"])
    print(head(sort(study_counts, decreasing = TRUE), 10))
    
  } else {
    cat("No enhanced details to write\n")
  }
} else {
  cat("No WGS entities found\n")
}

# Re-run extraction but for all DNA & RNA sequence data
# Define the process patterns for DNA sequencing and analysis
dna_patterns <- c(
  "assay.*wgs",           # Whole genome sequencing assays
  "assay.*gbs",           # Genotyping-by-sequencing assays  
  "analysis.*genome_assembly"  # Genome assembly analyses
)

# Define the process patterns for RNA sequencing and analysis
rna_patterns <- c(
  "assay.*rnaseq",                    # RNA sequencing assays
  "analysis.*rnaseq_quantification"   # RNA-seq quantification analyses
)

# Enhanced function to extract entities from multiple process types
extract_multi_type_entities <- function(nodes_df, edges_df, process_patterns, output_prefix) {
  
  cat("=== EXTRACTING ENTITIES FOR:", paste(process_patterns, collapse = ", "), "===\n")
  
  all_entities <- data.frame()
  
  # Find entities for each process pattern
  for(pattern in process_patterns) {
    cat("Searching for pattern:", pattern, "\n")
    entities <- find_entities_by_process_type(nodes_df, edges_df, pattern)
    
    if(nrow(entities) > 0) {
      entities$pattern_matched <- pattern
      all_entities <- rbind(all_entities, entities)
      cat("  Found", nrow(entities), "entities for pattern:", pattern, "\n")
    } else {
      cat("  No entities found for pattern:", pattern, "\n")
    }
  }
  
  if(nrow(all_entities) == 0) {
    cat("No entities found for any patterns\n")
    return()
  }
  
  # Remove duplicates (same entity from multiple patterns)
  all_entities <- all_entities[!duplicated(all_entities$id), ]
  cat("Total unique entities after deduplication:", nrow(all_entities), "\n")
  
  # Enhanced details with study and investigation info
  enhanced_details <- c()
  accession_names <- c()
  
  cat("Processing entities and finding parent studies/investigations...\n")
  
  for(i in 1:nrow(all_entities)) {
    if(i %% 50 == 0) cat("Processed", i, "of", nrow(all_entities), "entities...\n")
    
    entity_id <- all_entities$id[i]
    
    # Get accession name
    accession <- get_accession_name(entity_id, nodes_df, edges_df)
    clean_accession <- gsub(" \\(from parent.*\\)", "", as.character(accession))
    
    # Only process entities with valid accession names
    if(length(clean_accession) > 0 && 
       clean_accession != "No accession name found" && 
       clean_accession != "Entity not found" &&
       clean_accession != "NA" &&
       !is.na(clean_accession) &&
       nchar(clean_accession) > 0) {
      
      accession_names <- c(accession_names, clean_accession)
      
      # Find parent study and investigation
      parents <- find_parent_study_and_investigation(entity_id, nodes_df, edges_df)
      
      # Get readable names for study and investigation
      study_name <- get_process_name(parents$study_id, nodes_df)
      investigation_name <- get_process_name(parents$investigation_id, nodes_df)
      
      # Get additional entity metadata
      entity_row <- all_entities[all_entities$id == entity_id, ]
      sequencing_method <- if(!is.na(entity_row$sequencing_method)) entity_row$sequencing_method else "Not specified"
      origin <- if(!is.na(entity_row$origin)) entity_row$origin else "Not specified"
      provider <- if(!is.na(entity_row$provider)) entity_row$provider else "Not specified"
      
      # Create enhanced detail row
      enhanced_details <- c(enhanced_details, 
                            paste0(entity_id, "\t", 
                                   clean_accession, "\t",
                                   all_entities$parent_process_id[i], "\t",
                                   all_entities$parent_process_type[i], "\t",
                                   all_entities$pattern_matched[i], "\t",
                                   sequencing_method, "\t",
                                   origin, "\t",
                                   provider, "\t",
                                   parents$study_id, "\t",
                                   study_name, "\t",
                                   parents$investigation_id, "\t",
                                   investigation_name))
    }
  }
  
  # Create enhanced output files
  unique_accessions <- sort(unique(as.character(accession_names)))
  unique_accessions <- unique_accessions[!is.na(unique_accessions) & nchar(unique_accessions) > 0]
  
  cat("Found", length(unique_accessions), "unique accession names\n")
  
  if(length(enhanced_details) > 0) {
    # Write enhanced detailed file
    enhanced_header <- "Entity_ID\tAccession_Name\tProcess_ID\tProcess_Type\tPattern_Matched\tSequencing_Method\tOrigin\tProvider\tStudy_ID\tStudy_Name\tInvestigation_ID\tInvestigation_Name"
    writeLines(c(enhanced_header, enhanced_details), paste0(output_prefix, "_enhanced_detailed.txt"))
    
    # Also update the simple accession names file
    writeLines(unique_accessions, paste0(output_prefix, "_accession_names.txt"))
    
    cat("Created files:\n")
    cat("-", paste0(output_prefix, "_accession_names.txt"), "(just accession names)\n") 
    cat("-", paste0(output_prefix, "_enhanced_detailed.txt"), "(with study and investigation info)\n")
    
    # Show a preview of the enhanced data
    cat("\nPreview of enhanced data (first 3 entries):\n")
    cat("Entity_ID\tAccession_Name\tProcess_ID\tProcess_Type\tPattern_Matched\tSequencing_Method\tOrigin\tProvider\tStudy_ID\tStudy_Name\tInvestigation_ID\tInvestigation_Name\n")
    for(i in 1:min(3, length(enhanced_details))) {
      cat(enhanced_details[i], "\n")
    }
    
    # Summary statistics
    enhanced_df <- read.table(text = paste(c(enhanced_header, enhanced_details), collapse = "\n"), 
                              header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
    
    cat("\n=== SUMMARY STATISTICS ===\n")
    cat("Total entities:", nrow(enhanced_df), "\n")
    cat("Unique studies found:", length(unique(enhanced_df$Study_ID[enhanced_df$Study_ID != "Not found"])), "\n")
    cat("Unique investigations found:", length(unique(enhanced_df$Investigation_ID[enhanced_df$Investigation_ID != "Not found"])), "\n")
    
    cat("\nProcess types breakdown:\n")
    print(table(enhanced_df$Process_Type))
    
    cat("\nPattern matched breakdown:\n")
    print(table(enhanced_df$Pattern_Matched))
    
    cat("\nStudies with this data type:\n")
    study_counts <- table(enhanced_df$Study_Name[enhanced_df$Study_Name != "Not found"])
    print(head(sort(study_counts, decreasing = TRUE), 10))
    
    cat("\nSequencing methods found:\n")
    method_counts <- table(enhanced_df$Sequencing_Method[enhanced_df$Sequencing_Method != "Not specified"])
    if(length(method_counts) > 0) {
      print(method_counts)
    } else {
      cat("No sequencing methods specified\n")
    }
    
    return(enhanced_df)
    
  } else {
    cat("No enhanced details to write\n")
    return(NULL)
  }
} 

# Extract DNA sequencing and analysis entities
cat("=== EXTRACTING DNA SEQUENCING AND ANALYSIS ENTITIES ===\n")

dna_results <- extract_multi_type_entities(nodes_df, edges_df, dna_patterns, "dna_sequencing_analysis")

# Extract RNA sequencing and analysis entities  
cat("=== EXTRACTING RNA SEQUENCING AND ANALYSIS ENTITIES ===\n")

rna_results <- extract_multi_type_entities(nodes_df, edges_df, rna_patterns, "rna_sequencing_analysis")

# Create a combined summary
cat("=== COMBINED SUMMARY ===\n")

if(!is.null(dna_results)) {
  cat("DNA Sequencing/Analysis entities:", nrow(dna_results), "\n")
  cat("DNA Unique accessions:", length(unique(dna_results$Accession_Name)), "\n")
}

if(!is.null(rna_results)) {
  cat("RNA Sequencing/Analysis entities:", nrow(rna_results), "\n")
  cat("RNA Unique accessions:", length(unique(rna_results$Accession_Name)), "\n")
}

# Check for overlap between DNA and RNA datasets
if(!is.null(dna_results) && !is.null(rna_results)) {
  dna_accessions <- unique(dna_results$Accession_Name)
  rna_accessions <- unique(rna_results$Accession_Name)
  
  overlap <- intersect(dna_accessions, rna_accessions)
  cat("Accessions with both DNA and RNA data:", length(overlap), "\n")
  
  if(length(overlap) > 0) {
    cat("Examples of accessions with both data types:\n")
    for(i in 1:min(5, length(overlap))) {
      cat(" -", overlap[i], "\n")
    }
    
    # Save the overlap list
    writeLines(sort(overlap), "accessions_with_both_dna_and_rna.txt")
    cat("Saved overlap list to: accessions_with_both_dna_and_rna.txt\n")
  }
}

# Create a master list of all sequencing accessions
if(!is.null(dna_results) && !is.null(rna_results)) {
  all_accessions <- sort(unique(c(dna_results$Accession_Name, rna_results$Accession_Name)))
  writeLines(all_accessions, "all_sequencing_accession_names.txt")
  cat("Created master list: all_sequencing_accession_names.txt with", length(all_accessions), "unique accessions\n")
} else if(!is.null(dna_results)) {
  all_accessions <- sort(unique(dna_results$Accession_Name))
  writeLines(all_accessions, "all_sequencing_accession_names.txt")
  cat("Created master list: all_sequencing_accession_names.txt with", length(all_accessions), "unique DNA accessions\n")
} else if(!is.null(rna_results)) {
  all_accessions <- sort(unique(rna_results$Accession_Name))
  writeLines(all_accessions, "all_sequencing_accession_names.txt")
  cat("Created master list: all_sequencing_accession_names.txt with", length(all_accessions), "unique RNA accessions\n")
}

# Create summary tables for DNA and RNA
cat("=== CREATING SUMMARY TABLES ===\n")

# DNA Summary Table
if(!is.null(dna_results)) {
  cat("DNA Summary Table:\n")
  
  # Create summary by pattern
  dna_summary <- data.frame(
    Type = character(),
    Num_Entities = numeric(),
    Unique_Accessions = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(pattern in dna_patterns) {
    pattern_data <- dna_results[dna_results$Pattern_Matched == pattern, ]
    
    # Extract type name (remove "assay.*" or "analysis.*" prefix)
    type_name <- gsub("assay\\.*|analysis\\.*", "", pattern)
    
    num_entities <- nrow(pattern_data)
    unique_accessions <- length(unique(pattern_data$Accession_Name))
    
    dna_summary <- rbind(dna_summary, data.frame(
      Type = type_name,
      Num_Entities = num_entities,
      Unique_Accessions = unique_accessions,
      stringsAsFactors = FALSE
    ))
  }
  
  print(dna_summary)
  write.csv(dna_summary, "dna_summary_table.csv", row.names = FALSE)
  cat("DNA summary saved to: dna_summary_table.csv\n\n")
} else {
  cat("No DNA results to summarize\n\n")
}

# RNA Summary Table
if(!is.null(rna_results)) {
  cat("RNA Summary Table:\n")
  
  # Create summary by pattern
  rna_summary <- data.frame(
    Type = character(),
    Num_Entities = numeric(),
    Unique_Accessions = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(pattern in rna_patterns) {
    pattern_data <- rna_results[rna_results$Pattern_Matched == pattern, ]
    
    # Extract type name (remove "assay.*" or "analysis.*" prefix)
    type_name <- gsub("assay\\.*|analysis\\.*", "", pattern)
    
    num_entities <- nrow(pattern_data)
    unique_accessions <- length(unique(pattern_data$Accession_Name))
    
    rna_summary <- rbind(rna_summary, data.frame(
      Type = type_name,
      Num_Entities = num_entities,
      Unique_Accessions = unique_accessions,
      stringsAsFactors = FALSE
    ))
  }
  
  print(rna_summary)
  write.csv(rna_summary, "rna_summary_table.csv", row.names = FALSE)
  cat("RNA summary saved to: rna_summary_table.csv\n")
} else {
  cat("No RNA results to summarize\n")
}