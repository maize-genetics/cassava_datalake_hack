# Conversational Cassava Data Assistant with conversation memory
# Bethany Econopouly, econopouly@cornell.edu
# With help from Claude-4
# In progress script (test)
# Use with knowledge_graph_entitiesLevels.R to get cassava_nodes.csv and cassava_edges.csv
# Load libraries
library(httr2)
library(jsonlite)
library(dplyr)

# Load data from knowledge graph
nodes_df <- read.csv("~/path/to/cassava_nodes.csv", stringsAsFactors = FALSE)
edges_df <- read.csv("~/path/to/cassava_edges.csv", stringsAsFactors = FALSE)

# Set your API key
Sys.setenv(AI_API_KEY = "*****")

# CONVERSATION MEMORY STORAGE
conversation_history <- list()

# Function to call Claude with full conversation history
call_claude_with_history <- function(messages, system_prompt = NULL, model = "anthropic.claude-4-sonnet") {
  
  # Build request body with ALL messages (conversation history)
  body <- list(
    model = model,
    max_tokens = 2000,
    messages = messages  # This now contains full conversation
  )
  
  # Add system prompt
  if(!is.null(system_prompt)) {
    body$system <- system_prompt
  }
  
  # Make API request
  tryCatch({
    resp <- request("https://api.ai.it.cornell.edu/v1/messages") %>%
      req_headers(
        "Content-Type" = "application/json",
        "Authorization" = paste("Bearer", Sys.getenv("CORNELL_AI_API_KEY"))
      ) %>%
      req_body_json(body) %>%
      req_perform()
    
    # Parse response
    result <- resp %>% resp_body_json()
    
    # Extract content
    if("content" %in% names(result) && length(result$content) > 0) {
      return(result$content[[1]]$text)
    } else if("choices" %in% names(result)) {
      return(result$choices[[1]]$message$content)
    } else {
      return(as.character(result))
    }
    
  }, error = function(e) {
    if(inherits(e, "httr2_http")) {
      error_content <- resp_body_string(e$resp)
      return(paste("API Error:", resp_status(e$resp), "-", error_content))
    } else {
      return(paste("Request Error:", e$message))
    }
  })
}

# Function to create system context
create_data_context <- function() {
  paste0(
    "You are helping analyze a cassava genomic knowledge graph dataset.\n\n",
    
    "DATA STRUCTURE:\n",
    "- nodes_df columns: ", paste(colnames(nodes_df), collapse = ", "), "\n",
    "- edges_df columns: ", paste(colnames(edges_df), collapse = ", "), "\n",
    "- Total nodes: ", nrow(nodes_df), " (", sum(nodes_df$node_type == "entity"), " entities, ", sum(nodes_df$node_type == "process"), " processes)\n\n",
    
    "GRAPH STRUCTURE:\n",
    "- Processes connect TO entities via edges (from=process, to=entity)\n",
    "- Use pattern: filter processes -> get IDs -> filter edges -> get connected entities\n\n",
    
    "EXAMPLE R CODE:\n",
    "wgs_processes <- nodes_df %>% filter(process_type == 'assay|wgs')\n",
    "wgs_process_ids <- wgs_processes$id\n",
    "wgs_entity_ids <- edges_df %>% filter(from %in% wgs_process_ids) %>% pull(to)\n",
    "wgs_count <- length(unique(wgs_entity_ids))\n\n",
    
    "You can have conversations and generate R code when requested."
  )
}

# Function to execute R code
safe_execute <- function(code_string) {
  tryCatch({
    result <- eval(parse(text = code_string))
    return(list(success = TRUE, result = result, error = NULL))
  }, error = function(e) {
    return(list(success = FALSE, result = NULL, error = as.character(e)))
  })
}

# System context
system_context <- create_data_context()

# MAIN FUNCTION WITH CONVERSATION MEMORY
ask_cassava <- function(user_question, execute_code = TRUE) {
  
  tryCatch({
    cat("\n💭 Thinking...\n")
    
    # Add current question to conversation history
    conversation_history <<- append(conversation_history, 
                                    list(list(role = "user", content = user_question)))
    
    # Call Claude with FULL conversation history
    response_text <- call_claude_with_history(conversation_history, system_context)
    
    # Check for errors
    if(grepl("API Error:|Request Error:", response_text)) {
      cat("\n❌", response_text, "\n")
      return(list(error = response_text, success = FALSE))
    }
    
    # Add Claude's response to conversation history
    conversation_history <<- append(conversation_history, 
                                    list(list(role = "assistant", content = response_text)))
    
    cat("\n Claude:\n")
    cat(response_text, "\n")
    
    # Execute R code if present
    if(execute_code && grepl("```r|```R", response_text)) {
      
      r_code <- gsub(".*```r\\s*\n", "", response_text)
      r_code <- gsub("\n```.*", "", r_code)
      
      if(nchar(trimws(r_code)) > 0) {
        cat("\n⚡ Executing R code...\n")
        cat("---\n")
        cat(r_code)
        cat("\n---\n\n")
        
        result <- safe_execute(r_code)
        
        if(result$success) {
          cat(" Results:\n")
          print(result$result)
          
          # Add code execution result to conversation memory
          exec_summary <- paste("Code executed successfully. Result:", capture.output(print(result$result))[1])
          conversation_history <<- append(conversation_history, 
                                          list(list(role = "user", content = exec_summary)))
          
          return(list(response = response_text, code = r_code, result = result$result, success = TRUE))
        } else {
          cat("❌ Code execution error:\n")
          cat(result$error, "\n")
          return(list(response = response_text, code = r_code, error = result$error, success = FALSE))
        }
      }
    }
    
    return(list(response = response_text, success = TRUE))
    
  }, error = function(e) {
    cat("❌ Function Error:", e$message, "\n")
    return(list(error = paste("Function Error:", e$message), success = FALSE))
  })
}

# CONVERSATION MANAGEMENT FUNCTIONS
reset_conversation <- function() {
  conversation_history <<- list()
  cat("Conversation history cleared.\n")
}

show_conversation <- function() {
  if(length(conversation_history) == 0) {
    cat("No conversation history.\n")
  } else {
    cat("Conversation History:\n")
    cat("========================\n")
    for(i in 1:length(conversation_history)) {
      msg <- conversation_history[[i]]
      role_emoji <- if(msg$role == "user") "👤" else "🤖"
      preview <- substr(msg$content, 1, 80)
      if(nchar(msg$content) > 80) preview <- paste0(preview, "...")
      cat(sprintf("%s %s: %s\n", role_emoji, toupper(msg$role), preview))
    }
  }
}

# Initialize
cat("Cassava Genomic Data Assistant (with Memory!)\n")
cat("=================================================\n")
cat("Use reset_conversation() to start fresh\n")
cat("Use show_conversation() to see history\n\n")

# Test it out - now with memory!
ask_cassava("Hello! What can you tell me about this cassava dataset?")
ask_cassava("Can you tell me how many entities have WGS")
ask_cassava("Can you tell me how many of those entities with WGS data have accession names?")
ask_cassava("294 is not correct. Are you looking at parent entities of the entities with wgs to get the accession name from the biosample data?")
ask_cassava("could you make a csv file with all 2790 entities with wgs data and include metadata such as accession name that you can get from that entity, parent entity, and parent processes from the parent entity?")
ask_cassava("5282 is more than you originally found (2790)- why is that?")
ask_cassava("can you make a new csv file with the correction?")
ask_cassava("thank you, great work.")

head(wgs_metadata_corrected)
