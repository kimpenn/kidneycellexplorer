
selectizeIt <- shiny:::selectizeIt
checkAsIs <- shiny:::checkAsIs
htmlDependency <- htmltools:::htmlDependency
toJSON <- shiny:::toJSON
attachDependencies <- htmltools:::attachDependencies
choicesWithNames <- shiny:::choicesWithNames
selectOptions <- shiny:::selectOptions
#controlLabel <- shiny:::controlLabel
controlLabel <- shiny:::shinyInputLabel
htmlEscape <- htmltools:::htmlEscape

selectizeInput2 <- function (inputId, ..., options = NULL, width = NULL) 
{
    selectizeIt2(inputId, selectInput2(inputId, ..., selectize = FALSE, 
                                     width = width), options)
}

selectInput2 <- function (inputId, label, choices, selected = NULL, disabled = NULL,  multiple = FALSE, 
          selectize = TRUE, width = NULL, size = NULL) 
{
    selected <- restoreInput(id = inputId, default = selected)
    choices <- choicesWithNames(choices)
    if (is.null(selected)) {
        if (!multiple) 
            selected <- firstChoice(choices)
    }
    else selected <- as.character(selected)
    if (!is.null(size) && selectize) {
        stop("'size' argument is incompatible with 'selectize=TRUE'.")
    }
    selectTag <- tags$select(id = inputId, class = if (!selectize) 
        "form-control", size = size, selectOptions2(choices, selected, disabled))
    if (multiple) 
        selectTag$attribs$multiple <- "multiple"
    res <- div(class = "form-group shiny-input-container", style = if (!is.null(width)) 
        paste0("width: ", validateCssUnit(width), ";"), controlLabel(inputId, 
                                                                     label), div(selectTag))
    if (!selectize) 
        return(res)
    selectizeIt(inputId, res, NULL, nonempty = !multiple && !("" %in% 
                                                                  choices))
}

selectOptions2 <- function (choices, selected = NULL, disabled = NULL) 
{
    html <- mapply(choices, names(choices), FUN = function(choice, 
                                                           label) {
        if (is.list(choice)) {
            sprintf("<optgroup label=\"%s\">\n%s\n</optgroup>", 
                    htmlEscape(label, TRUE), selectOptions2(choice, 
                                                           selected, disabled))
        }
        else {
            sprintf("<option value=\"%s\"%s>%s</option>", 
                    htmlEscape(choice, TRUE), 
                    if (choice %in% selected) " selected" else if (choice %in% disabled) " disabled" else "", 
                    htmlEscape(label))
        }
    })
    HTML(paste(html, collapse = "\n"))
}


#' PIVOT help modules, UI
#'
#' @export
pivot_help_UI <- function(id, title, label = NULL, icn = "question-circle", type = "button", tooltip = T, style = NULL){
    ns<- NS(id)
    if(tooltip) {
        tip <- shinyBS::bsTooltip(
            ns("pivot_help"),
            title = title,
            options = list(container = "body")
        )
    } else {
        tip <- NULL
    }
    if(type == "button") {
        btn <-  actionButton(ns("pivot_help"), label = label, icon = icon(icn), style = style)
    } else {
        btn <- actionLink(ns("pivot_help"), label = label, icon = icon(icn), style = style)
    }
    tagList(
        btn,
        tip
    )
}

#' PIVOT help modules, server
#'
#' @export
pivot_help <- function (input, output, session, title, content, size = "m") {
    observeEvent(input$pivot_help, {
        showModal(modalDialog(
            title = title,
            size = size,
            content,
            easyClose = TRUE
        ))
    })
}


