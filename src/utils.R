#' Fix and simplify polygon geometry
#'
#' Cleans invalid polygon geometries by repairing topology issues,
#' removing collapsed geometries, selecting the largest polygon,
#' and stripping interior holes.
#'
#' The geometry is:
#' \itemize{
#'   \item transformed to EPSG:3857 for robust geometric operations,
#'   \item snapped to a fixed precision grid,
#'   \item repaired using \code{st_make_valid()},
#'   \item split into polygon components,
#'   \item reduced to the largest polygon by area,
#'   \item stripped of interior holes,
#' }
#'
#' @param g An `sf` or `sfc` polygon geometry object.
#'
#' @return
#' An `sfc_POLYGON` object containing a single valid polygon geometry
#' with holes removed.
#'
#' @details
#' All functions used are from the sf package
#' 
fix_geometry <- function(g) {
  g_parts <- (
    g
    %>% st_transform(3857)
    %>% st_set_precision(100)
    %>% st_make_valid(geos_keep_collapsed = FALSE)
    %>% st_transform(4326)
    %>% st_cast("POLYGON")
  )
  g_clean <- g_parts[which.max(st_area(g_parts))]
  coords <- st_coordinates(g_clean)
  val_coords <- coords[coords[, "L1"] == 1, ]
  xy <- val_coords[, c("X", "Y")]
  g_no_holes <- st_sfc(st_polygon(list(xy)), crs = st_crs(g))
  g_no_holes
}