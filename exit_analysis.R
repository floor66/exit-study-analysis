rm(list=ls())

# Load database (CSV file)
raw_database_csv <- read.csv("20220426 EXIT database.csv", sep=";", header=TRUE)

# Calculate pooled proportions for outcomes?
pooled_proportions <- FALSE
pooled_proportion_outcomes <- c("n_extraction_site_hernia", "n_ssi", "n_sso", "n_wound_dehiscense")

# Save figures as .svg files automagically
save <- FALSE







































library(shiny)
library(dplyr)
library(ggplot2)
library(meta)
library(robvis)
library(grid)
library(svglite)

capitalize <- function(x) {
  return(paste(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)), sep=""));
}

# Exclude excluded studies
excludes <- raw_database_csv[which(raw_database_csv$exclude_include==1),]
n_excludes <- length(unique(excludes$index))
print(paste("Excluded ", n_excludes ," studies", sep=""))
database <- raw_database_csv[which(raw_database_csv$exclude_include==0),]

# Fix holes in data
for (row in 1:nrow(database)) {
  idx <- database[row, "index"]
  parents <- database[which(database$index == idx),]
  
  if(database[row, "first_author"] == "") {
    database[row, "first_author"] <- unique(parents[which(parents$first_author != ""), "first_author"])
  }
  
  if(is.na(database[row, "year"])) {
    database[row, "year"] <- unique(parents[which(!is.na(parents$year)), "year"])
  }
}

# Construct study labels out of first author + year
database$label <- paste(database$first_author, ", ", database$year, sep="")

# Create a readable, factored study type variable
database$study_type_real <- ifelse(
  database$study_type == 1, "Randomized controlled trial",
  ifelse(
    database$study_type == 2, "Prospective cohort",
    ifelse(
      database$study_type == 3, "Retrospective cohort", "NA"
    )
  )
)

# Create study type titles for subgroup analysis
database$study_type_real <- factor(
  database$study_type_real,
  levels = c("Randomized controlled trial", "Prospective cohort", "Retrospective cohort")
)

# Risk of bias data
robs <- database[which(nchar(database$title) > 0), c(grep("label|rob_*", names(database)))]

for(i in 1:8) {
  rob <- robs[[paste("rob_", ifelse(i == 8, "overall", paste("d", i, sep="")), sep="")]]
  rob <- trimws(rob[which(nchar(rob) > 0)])
  
  print(paste(sep="", "Domain ", i, ":"))
  print(paste(sep="", "Low: ", (sum(rob == "Low") / length(rob) * 100), "%"))
  print(paste(sep="", "Moderate: ", (sum(rob == "Moderate") / length(rob) * 100), "%"))
  print(paste(sep="", "Serious: ", (sum(rob == "Serious") / length(rob) * 100), "%"))
  print(paste(sep="", "Critical: ", (sum(rob == "Critical") / length(rob) * 100), "%"))
  print("")
}

# write.csv(robs, file = "robs.csv")

# Remove studies if they have <10 patients
database <- database[which(database$n_subgroup >= 5),]

# Remove 'open' study groups
database <- database[which(regexpr("pen", database$type_of_surgery) == -1),]

# Pooled proportions
pooled_forest <- function(database,
                          discriminator, discriminator_label,
                          experiment_val, experiment_label,
                          outcome, outcome_label,
                          save = FALSE) {

  experiment <- database[which(database[[discriminator]] == experiment_val),]
  experiment <- experiment[which(!is.na(experiment[[outcome]])),]
  
  if(nrow(experiment) > 1) {
    pooled_props <- metaprop(
      studlab = experiment$label,
      event = as.numeric(experiment[[outcome]]),
      n = experiment$n_subgroup,
      method = "Inverse",
      method.tau = "DL",
      comb.fixed = FALSE,
      comb.random = TRUE,
      prediction = TRUE
    )
    
    print(paste(length(unique(experiment$label)), "studies,", sum(experiment$n_subgroup, na.rm=TRUE), "patients"))
    
    if(save) {
      filename <- paste("Figures/by_category/", discriminator_label, "=", experiment_label, "_", outcome_label, ".svg", sep="")
      svglite(filename, width = 15, height = 15, pointsize = 12)
    }
    
    # Forest plot
    forest(
      pooled_props,
      digits.pval = 3,
      digits = 3,
      xlim = c(0, 0.5)
    )
    
    if(save) {
      grid.text(
        paste(discriminator_label, ": ", experiment_label, "\n", outcome_label, sep=""),
        0.5, 0.95, gp=gpar(cex=2)
      )

      dev.off()
      
      writeLines(
        gsub(pattern = "textLength='(.*)' ", replace = "", readLines(filename)),
        con=filename
      )
    }
    
    return(pooled_props)
  } else {
    print(paste("WARNING: ", "No data available for ", outcome_label, " (", discriminator, "=", experiment_val ,")", sep=""));
  }
}

# By extraction site
if(pooled_proportions) {
  categories <- unique(database$extraction_site_category)

  for(outcome in pooled_proportion_outcomes) {
    for(cat in categories) {
      pooled_forest(
        database,
        "extraction_site_category",
        "Extraction site",
        cat,
        capitalize(cat),
        outcome,
        capitalize(gsub("^n_", "", outcome)),
        save = save
      )
    }
  }
}

# Comparative studies
comparative_forest = function(database, experiment, experiment_label, control, control_label,
                              outcome, outcome_label,
                              save = FALSE, aggregate = FALSE, aggregate_by = NULL) {

  # Only comparative studies
  intersection <- intersect(experiment$label, control$label)
  experiment <- experiment[experiment$label %in% intersection,]
  control <- control[control$label %in% intersection,]
  
  experiment[[outcome]] <- as.numeric(experiment[[outcome]])
  control[[outcome]] <- as.numeric(control[[outcome]])
  
  experiment <- experiment[which(!is.na(experiment[[outcome]])),]
  control <- control[which(!is.na(control[[outcome]])),]
  
  if(nrow(experiment) > 1 && nrow(control) > 1) {
    if(aggregate) {
      # Aggregate midline and non-midline based on 'extraction_in_midline'
      agg_cols <- c(
        "index",
        "label",
        aggregate_by,
        "study_type_real",
        outcome,
        "n_subgroup"
      )
      
      experiment <- aggregate(
        cbind(experiment[[outcome]], experiment$n_subgroup),
        by = list(
          experiment$index,
          experiment$label,
          experiment[[aggregate_by]],
          experiment$study_type_real
        ),
        FUN = sum
      )
      
      control <- aggregate(
        cbind(control[[outcome]], control$n_subgroup),
        by = list(
          control$index,
          control$label,
          control[[aggregate_by]],
          control$study_type_real
        ),
        FUN = sum
      )
      
      colnames(experiment) <- agg_cols
      colnames(control) <- agg_cols
    }
    
    # Fix follow-up periods
    for(row in 1:nrow(experiment)) {
      fus <- as.numeric(gsub(",", ".", database[which(database$index == experiment[row, "index"]), "follow_up_simple"]))
      experiment[row, "follow_up_simple"] <- round(min(fus))
    }
    
    for(row in 1:nrow(control)) {
      fus <- as.numeric(gsub(",", ".", database[which(database$index == control[row, "index"]), "follow_up_simple"]))
      control[row, "follow_up_simple"] <- round(min(fus))
    }
    
    experiment[which(is.na(experiment$follow_up_simple)), "follow_up_simple"] <- "*"
    control[which(is.na(control$follow_up_simple)), "follow_up_simple"] <- "*"
    
    # Comparative analysis
    
    comparative <- metabin(
      data = experiment,
      event.e = experiment[[outcome]],
      n.e = experiment$n_subgroup,
      event.c = control[[outcome]],
      n.c = control$n_subgroup,
      method = "Inverse",
      method.tau = "DL",
      comb.fixed = FALSE,
      comb.random = TRUE,
      studlab = experiment$label,
      byvar = experiment$study_type_real,
      sm = "OR",
      label.left = paste("Favours", tolower(experiment_label)),
      label.right = paste("Favours", tolower(control_label)),
      label.e = experiment_label,
      label.c = control_label,
      overall = TRUE,
      allstudies = TRUE,
      print.byvar = FALSE
    )
    
    tmp <- sum(experiment$n_subgroup) + sum(control$n_subgroup)
    print(paste(length(unique(experiment$label)), "studies,", tmp, "patients"))

    if(save) {
      filename <- paste("Figures/funnel_plots/", experiment_label, "_vs_", control_label, "_", outcome, ".svg", sep="")
      svglite(filename, width = 10, height = 10, pointsize = 12)
    }
    
    funnel(
      comparative,
      studlab = TRUE,
      contour = c(0.9, 0.95, 0.99),
      col.contour = c("gray75", "gray85", "gray95")
    )
    
    legend(
      "topleft",
      legend = c("p < 0.1", "p < 0.05", "p < 0.01"),
      fill = c("gray75", "gray85", "gray95"),
      cex = 0.75
    )
    
    if(save) {
      grid.text(
        paste(experiment_label, " versus ", control_label, "\n", outcome_label, sep=""),
        0.5, 0.9625, gp=gpar(cex=2)
      )
      
      dev.off()
    }
    
    if(save) {
      filename <- paste("Figures/comparative/", experiment_label, "_vs_", control_label, "_", outcome, ".svg", sep="")
      svglite(filename, width = 15, height = 15, pointsize = 12)
    }
    
    forest(
      x = comparative,
      digits.pval = 3,
      digits = 3,
      digits.addcols = c(0),
      test.effect.subgroup = TRUE,
      test.overall.random = TRUE,
      sortvar = experiment$follow_up_simple,
      leftcols = c("studlab", "follow_up_simple", "event.e", "n.e", "event.c", "n.c"),
      rightcols = c("w.random", "effect.ci"),
      leftlabs = c("Study", "Follow-up", "Event", "N", "Event", "N", "Weight", "OR [95%-CI]")
    )

    if(save) {
      grid.text(
        paste(experiment_label, " versus ", control_label, "\n", outcome_label, sep=""),
        0.5, 0.95, gp=gpar(cex=2)
      )
      
      dev.off()
      
      writeLines(
        gsub(pattern = "textLength='(.*)' ", replace = "", readLines(filename)),
        con=filename
      )
    }
    
    return(comparative)
  } else {
    print(paste("WARNING: ", "No data available for ", outcome_label, " (", experiment_label, " vs. ", control_label, ")", sep=""));
  }
}

# tmp <- database
# database <- database[which(database$extraction_site_category != "pfannenstiel"),]
# database <- tmp

# comparative_forest(
#   database,
#   database[which(database$extraction_in_midline == 1),],
#   "Midline",
#   database[which(database$extraction_in_midline == 0),],
#   "Non-midline",
#   "n_wound_dehiscense",
#   "WD",
#   aggregate = TRUE,
#   aggregate_by = "extraction_in_midline",
#   save = TRUE
# )

# v <- database$extraction_length_simple
# 
# v <- gsub(",", "\\.", v)
# v <- as.numeric(v)
# median(v, na.rm=TRUE)
# min(v, na.rm=TRUE)
# max(v, na.rm=TRUE)

# sum(database[which(database$oncological_surgery == 0),]$n_subgroup)
# sum(database[which(database$oncological_surgery == 1),]$n_subgroup)
# sum(database[which(database$oncological_surgery == 2),]$n_subgroup)


database$crosscomp_group <- ifelse(database$extraction_site_category == "midline", "midline",
                             ifelse(database$extraction_site_category == "umbilical", "midline",
                              ifelse(database$extraction_site_category == "transverse", "non-midline",
                               ifelse(database$extraction_site_category == "paramedian", "non-midline",
                                ifelse(database$extraction_site_category == "pfannenstiel", "pfannenstiel", NA)))))

database$suture_group <- ifelse(regexpr("PDS", database$fascial_closure_technique) != -1, "pds",
                          ifelse(regexpr("icryl", database$fascial_closure_technique) != -1, "vicryl",
                           ifelse(regexpr("xon", database$fascial_closure_technique) != -1, "maxon", "misc")))

# Shiny for proportions
if(TRUE) {
  shinyApp(
    ui = fluidPage(
      titlePanel("Forest plot creator"),
      sidebarLayout(
        sidebarPanel(
          h2("Pooled proportions"),
          selectInput("discriminator", "Variable to determine group:", colnames(database)),
          textInput("discriminator_label", "Label"),
          uiOutput("select_1"),
          textInput("experiment_label", "Label"),
          selectInput("outcome", "Outcome variable:", colnames(database)),
          textInput("outcome_label", "Label"),
          sliderInput("plotHeight", label = "Plot height:", min = 0, max = 2000, step = 10, value = 500),
          actionButton("save", "Save plot (.svg)")
        ),
        mainPanel(
          plotOutput("title", height = "100px"),
          plotOutput("forest")
        )
      )
    ),
    server = function(input, output) {
      output$select_1 <- renderUI({
        selectInput("experiment", "Which group:", unique(database[[input$discriminator]]))
      })

      plotHeight <- reactive(input$plotHeight)
      
      observeEvent(input$save, {
        pooled_forest(
          database,
          input$discriminator,
          input$discriminator_label,
          input$experiment,
          input$experiment_label,
          input$outcome,
          input$outcome_label,
          save = TRUE
        )
        
        showModal(modalDialog(
          title = "Save",
          "Plot saved successfully.",
          easyClose = TRUE,
          footer = NULL
        ))
      })
      
      output$title <- renderPlot({
        grid.text(
          paste(input$discriminator_label, ": ", input$experiment_label, "\n", input$outcome_label, sep=""),
          0.5, 0.5, gp=gpar(cex=2)
        )
      })
      
      output$forest <- renderPlot({
        pooled_forest(
          database,
          input$discriminator,
          input$discriminator_label,
          input$experiment,
          input$experiment_label,
          input$outcome,
          input$outcome_label
        )
      }, height = function() plotHeight())
    }
  )
}

# Shiny for comparative
if(TRUE) {
  shinyApp(
    ui = fluidPage(
      titlePanel("Forest plot creator"),
      sidebarLayout(
        sidebarPanel(
          h2("Comparative"),
          selectInput("discriminator", "Variable to determine groups:", colnames(database)),
          uiOutput("select_1"),
          textInput("experiment_label", "Label"),
          uiOutput("select_2"),
          textInput("control_label", "Label"),
          selectInput("outcome", "Outcome variable:", colnames(database)),
          textInput("outcome_label", "Label"),
          checkboxInput("aggregate", "Aggregate groups?", TRUE),
          uiOutput("select_3"),
          sliderInput("plotHeight", label = "Plot height:", min = 0, max = 2000, step = 10, value = 500),
          actionButton("save", "Save plot (.svg)")
        ),
        mainPanel(
          plotOutput("title", height = "100px"),
          plotOutput("forest")
        )
      )
    ),
    server = function(input, output) {
      output$select_1 <- renderUI({
        selectInput("experiment", "Experiment:", unique(database[[input$discriminator]]))
      })
      
      output$select_2 <- renderUI({
        selectInput("control", "Control:", unique(database[[input$discriminator]]))
      })

      output$select_3 <- renderUI({
        selectInput("aggregate_by", "Aggregate by:", colnames(database), selected = input$discriminator)
      })

      plotHeight <- reactive(input$plotHeight)
      
      observeEvent(input$save, {
        comparative_forest(
          database,
          database[which(database[[input$discriminator]] == input$experiment),],
          input$experiment_label,
          database[which(database[[input$discriminator]] == input$control),],
          input$control_label,
          input$outcome,
          input$outcome_label,
          aggregate = input$aggregate,
          aggregate_by = input$aggregate_by,
          save = TRUE
        )

        showModal(modalDialog(
          title = "Save",
          "Plot saved successfully.",
          easyClose = TRUE,
          footer = NULL
        ))
      })
      
      output$title <- renderPlot({
        grid.text(
          paste(input$experiment_label, " versus ", input$control_label, "\n", input$outcome_label, sep=""),
          0.5, 0.5, gp=gpar(cex=2)
        )
      })
      
      output$forest <- renderPlot({
        result <- comparative_forest(
          database,
          database[which(database[[input$discriminator]] == input$experiment),],
          input$experiment_label,
          database[which(database[[input$discriminator]] == input$control),],
          input$control_label,
          input$outcome,
          input$outcome_label,
          aggregate = input$aggregate,
          aggregate_by = input$aggregate_by
        )
      }, height = function() plotHeight())
    }
  )
}
