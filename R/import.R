#' Import MSigDB XML database file into glacier-specific format
#'
#' @param path path to file
#'
#' @return glacier-specific imported database
#' @keywords internal
#'
#' @importFrom magrittr %>%
#' @examples \dontrun{
#' data <- import_msigdb(msig_path)
#' }
import_msigdb <- function(path) {
  data <- path %>% xml2::read_xml() %>% xml2::xml_children()

  # extract information
  gs_info <- tibble::tibble(
    name = xml2::xml_attr(data, "STANDARD_NAME"),
    info = xml2::xml_attr(data, "DESCRIPTION_BRIEF"),
    desc = xml2::xml_attr(data, "DESCRIPTION_FULL") %>% as.factor(),
    category = stringr::str_c(xml2::xml_attr(data, "CATEGORY_CODE"),
                              xml2::xml_attr(data, "SUB_CATEGORY_CODE"),
                              sep = " ") %>%
      stringr::str_squish() %>%
      as.factor(),
    organism = xml2::xml_attr(data, "ORGANISM") %>% as.factor()
  )

  # extract genes
  gs_genes <- xml2::xml_attr(data, "MEMBERS_SYMBOLIZED") %>%
    stringr::str_split(",") %>%
    purrr::map(~.[. != ""]) %>%
    purrr::set_names(gs_info$name)

  list(gs_genes = gs_genes, gs_info = gs_info)
}
