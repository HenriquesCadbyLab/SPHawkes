#' Animate previously evaluated intensity snapshots
#'
#' Creates a GIF from the snapshots returned by background_trigger_grid().
#' Each frame displays one evaluated model time. No intermediate intensity
#' surfaces are generated, and the model is not refitted or reevaluated.
#'
#' All frames share the same colour scale.
#'
#' @param grid A list returned by background_trigger_grid().
#' @param component One of "total", "background", or "triggering".
#' @param times Optional distinct evaluated model times to animate.
#'   If NULL, uses all times in grid$days. Selected times are sorted.
#' @param filename Output GIF filename. Defaults to a new file in the R
#'   session's temporary directory, so nothing is written elsewhere unless you
#'   choose a location; the path is returned.
#' @param fps Frames per second.
#' @param width,height Output dimensions in pixels.
#' @param cell_size Optional tile size passed to plot_intensity_snapshot().
#' @param overwrite Whether to replace an existing output file.
#'
#' @return The output filename, invisibly.
#' @examples
#' \donttest{
#' # A quick fit to one year of the bundled simulated data: a coarse grid, a
#' # coarse triggering-kernel grid and two EM iterations keep it to a few
#' # seconds. Real analyses need finer settings and more iterations.
#' cache_dir <- tempfile("spHawkes_marks_")
#' fit <- fit_semihawkes(
#'   sim_events, sim_region,
#'   date_col  = "date",
#'   end_date  = "2022-01-01",
#'   cell_size = 1000,
#'   cache_dir = cache_dir,
#'   trigger_spatial_args = list(x_range = 2000, x_step = 250),
#'   trigger_window_dist  = 2000,
#'   n_iter  = 2,
#'   verbose = FALSE,
#'   seed    = 1
#' )
#' grid <- background_trigger_grid(fit, sim_region, cell_size = 500,
#'                                 days = seq(100, 130, by = 10))
#' if (requireNamespace("gifski", quietly = TRUE)) {
#'   # Written to a temporary file unless `filename` is given
#'   gif <- animate_intensity_snapshots(grid, component = "triggering",
#'                                      fps = 2, width = 400, height = 300)
#'   file.exists(gif)
#' }
#' }
#' @export
animate_intensity_snapshots <- function(
    grid,
    component = c("total", "background", "triggering"),
    times = NULL,
    filename = tempfile("intensity_", fileext = ".gif"),
    fps = 5,
    width = 900,
    height = 700,
    cell_size = NULL,
    overwrite = FALSE) {

  component <- match.arg(component)

  for (pkg in c("ggplot2", "gifski")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        "Package '", pkg, "' is required for animation.",
        call. = FALSE
      )
    }
  }

  if (!is.list(grid) ||
      !all(c("days", "x", "y", "inside", "snapshots") %in% names(grid))) {
    stop(
      "'grid' must be an output of background_trigger_grid().",
      call. = FALSE
    )
  }

  if (!is.numeric(grid$days) || !length(grid$days) ||
      any(!is.finite(grid$days)) || anyDuplicated(grid$days) ||
      length(grid$snapshots) != length(grid$days)) {
    stop("The grid contains invalid times or snapshots.", call. = FALSE)
  }

  if (is.null(times)) {
    times <- grid$days
  }

  if (!is.numeric(times) || !length(times) ||
      any(!is.finite(times)) || anyDuplicated(times)) {
    stop("'times' must contain distinct finite numeric values.",
         call. = FALSE)
  }

  times <- sort(times)
  indices <- match(times, grid$days)

  if (anyNA(indices)) {
    stop("All requested times must be available in 'grid$days'.",
         call. = FALSE)
  }

  positive_scalar <- function(z, integer = FALSE) {
    is.numeric(z) &&
      length(z) == 1L &&
      is.finite(z) &&
      z > 0 &&
      (!integer || z == floor(z))
  }

  if (!positive_scalar(fps)) {
    stop("'fps' must be positive.", call. = FALSE)
  }

  if (!positive_scalar(width, TRUE) ||
      !positive_scalar(height, TRUE)) {
    stop("'width' and 'height' must be positive integers.",
         call. = FALSE)
  }

  if (!is.logical(overwrite) || length(overwrite) != 1L ||
      is.na(overwrite)) {
    stop("'overwrite' must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.character(filename) || length(filename) != 1L ||
      is.na(filename) ||
      !grepl("\\.gif$", filename, ignore.case = TRUE)) {
    stop("'filename' must be a GIF filename.", call. = FALSE)
  }

  if (!dir.exists(dirname(filename))) {
    stop("The output directory does not exist.", call. = FALSE)
  }

  if (file.exists(filename) && !overwrite) {
    stop(
      "The output file exists. Use 'overwrite = TRUE' to replace it.",
      call. = FALSE
    )
  }

  # Validate the plotting inputs before creating temporary files.
  first_plot <- plot_intensity_snapshot(
    grid = grid,
    time = times[1L],
    component = component,
    cell_size = cell_size
  )

  expected_dim <- c(length(grid$x), length(grid$y))

  frame_maxima <- vapply(
    indices,
    function(i) {
      snapshot <- grid$snapshots[[i]]

      values <- switch(
        component,
        background = snapshot$bg,
        triggering = snapshot$tr,
        total = snapshot$bg + snapshot$tr
      )

      if (!is.matrix(values) || !is.numeric(values) ||
          !identical(dim(values), expected_dim)) {
        stop("A snapshot contains invalid intensity matrices.",
             call. = FALSE)
      }

      values <- as.vector(values)[grid$inside]

      if (any(!is.finite(values)) || any(values < 0)) {
        stop("Intensity values must be finite and non-negative.",
             call. = FALSE)
      }

      max(values)
    },
    numeric(1)
  )

  upper_limit <- max(frame_maxima)

  if (upper_limit == 0) {
    upper_limit <- 1
  }

  frame_dir <- tempfile("spHawkes-frames-")
  dir.create(frame_dir)

  on.exit(
    unlink(frame_dir, recursive = TRUE),
    add = TRUE
  )

  png_files <- file.path(
    frame_dir,
    sprintf("frame_%06d.png", seq_along(times))
  )

  for (i in seq_along(times)) {
    p <- if (i == 1L) {
      first_plot
    } else {
      plot_intensity_snapshot(
        grid = grid,
        time = times[i],
        component = component,
        cell_size = cell_size
      )
    }

    # Replace the frame-specific scale with the shared animation scale.
    p$scales$scales <- Filter(
      function(scale) !("fill" %in% scale$aesthetics),
      p$scales$scales
    )

    p <- p +
      ggplot2::scale_fill_viridis_c(
        limits = c(0, upper_limit),
        name = "Intensity"
      )

    ggplot2::ggsave(
      filename = png_files[i],
      plot = p,
      width = width / 100,
      height = height / 100,
      units = "in",
      dpi = 100,
      bg = "white"
    )
  }

  gifski::gifski(
    png_files = png_files,
    gif_file = filename,
    width = width,
    height = height,
    delay = 1 / fps,
    loop = TRUE,
    progress = interactive()
  )

  invisible(filename)
}
