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