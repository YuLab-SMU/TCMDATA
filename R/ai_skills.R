#' Skill directory management
#'
#' Functions to manage a user-local skills directory for customizing analysis
#' preferences and adding new skills.
#'
#' @name skill-management
#' @rdname skill-management
NULL

#' Get path to a bundled aisdk skill
#'
#' Returns the absolute path to a skill bundled with the \code{aisdk}
#' package, such as \code{"skill-creator"}. This is useful when combining
#' TCMDATA's own skills with upstream aisdk skills without copying them into
#' the TCMDATA package.
#'
#' @param name Character. Skill name in \code{aisdk/inst/skills}.
#'   Default \code{"skill-creator"}.
#'
#' @return The absolute path to the aisdk skill directory.
#'
#' @examples
#' \dontrun{
#'   skill_creator <- tcm_aisdk_skill()
#'   agent <- create_tcm_task_agent(
#'     skills = c(tcm_skill_dir(), skill_creator)
#'   )
#' }
#' @export
tcm_aisdk_skill <- function(name = "skill-creator") {
  .check_aisdk()

  path <- system.file(file.path("skills", name), package = "aisdk")
  if (!nzchar(path) || !dir.exists(path)) {
    stop("Cannot find aisdk skill: ", name, call. = FALSE)
  }

  normalizePath(path, mustWork = TRUE)
}

#' Initialize a local skills directory
#'
#' Copies all bundled skills from the TCMDATA package to a local directory.
#' Once initialized, the agent will use the local directory instead of the
#' package defaults, allowing full customization.
#'
#' @param path Character. Directory to create. Default \code{"tcm_skills"}.
#' @param overwrite Logical. Overwrite existing directory? Default FALSE.
#'
#' @return The absolute path to the skills directory (invisibly).
#'
#' @examples
#' \dontrun{
#'   # Copy all skills to ./tcm_skills/ for customization
#'   tcm_init_skills()
#'
#'   # Then edit the preferences:
#'   # file.edit("tcm_skills/analysis-preferences/SKILL.md")
#'
#'   # Or place new custom skills under tcm_skills/
#'   # aisdk's skill-creator remains available by default
#' }
#' @export
tcm_init_skills <- function(path = "tcm_skills", overwrite = FALSE) {
  path <- normalizePath(path, mustWork = FALSE)

  if (dir.exists(path) && !overwrite) {
    message("Skills directory already exists: ", path)
    message("Use overwrite = TRUE to replace, or tcm_use_skills() to activate it.")
    tcm_use_skills(path)
    return(invisible(path))
  }

  pkg_skills <- system.file("skills", package = "TCMDATA")
  if (!nzchar(pkg_skills) || !dir.exists(pkg_skills)) {
    stop("Cannot find bundled skills in TCMDATA package.")
  }

  # Copy the entire skills directory
  if (dir.exists(path)) unlink(path, recursive = TRUE)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  skill_dirs <- list.dirs(pkg_skills, full.names = TRUE, recursive = FALSE)
  for (sd in skill_dirs) {
    target <- file.path(path, basename(sd))
    dir.create(target, recursive = TRUE, showWarnings = FALSE)
    files <- list.files(sd, recursive = TRUE, full.names = TRUE)
    for (f in files) {
      rel <- sub(paste0(sd, "/"), "", f, fixed = TRUE)
      dest <- file.path(target, rel)
      dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
      file.copy(f, dest, overwrite = TRUE)
    }
  }

  tcm_use_skills(path)

  skill_names <- basename(skill_dirs)
  message("Initialized skills directory: ", path)
  message("Copied skills: ", paste(skill_names, collapse = ", "))
  message("Edit any SKILL.md to customize. New skills can be added as subdirectories.")

  invisible(path)
}


#' Set the active skills directory
#'
#' Points the agent to an existing skills directory. All subsequent calls to
#' \code{\link{tcm_agent}}, \code{\link{tcm_chat}}, and
#' \code{\link{create_tcm_task_agent}} will use skills from this directory.
#'
#' @param path Character. Path to an existing skills directory containing
#'   skill subdirectories with SKILL.md files.
#'
#' @return The absolute path (invisibly).
#'
#' @examples
#' \dontrun{
#'   tcm_use_skills("my_project/skills")
#' }
#' @export
tcm_use_skills <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  if (!dir.exists(path)) {
    stop("Directory does not exist: ", path,
         "\nUse tcm_init_skills() to create one.")
  }
  options(tcmdata.user_skills = path)
  invisible(path)
}


#' Get the active skills directory
#'
#' Returns the current skills directory path. If a user-local directory has been
#' set via \code{\link{tcm_init_skills}} or \code{\link{tcm_use_skills}}, that
#' path is returned. Otherwise, returns the package default.
#'
#' @return Character path to the active skills directory, or NULL if no skills
#'   are available.
#'
#' @examples
#' \dontrun{
#'   tcm_skill_dir()
#' }
#' @export
tcm_skill_dir <- function() {
  user_dir <- getOption("tcmdata.user_skills")
  if (!is.null(user_dir) && dir.exists(user_dir)) {
    return(user_dir)
  }
  pkg_skills <- system.file("skills", package = "TCMDATA")
  if (nzchar(pkg_skills) && dir.exists(pkg_skills)) {
    return(pkg_skills)
  }
  NULL
}


#' Reset skills to package defaults
#'
#' Clears the user skills directory setting so the agent uses the bundled
#' package skills again.
#'
#' @return NULL (invisibly).
#'
#' @examples
#' \dontrun{
#'   tcm_reset_skills()
#' }
#' @export
tcm_reset_skills <- function() {
  options(tcmdata.user_skills = NULL)
  message("Skills reset to package defaults.")
  invisible(NULL)
}
