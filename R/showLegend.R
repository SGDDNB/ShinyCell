#' Shows the legends for single-cell metadata
#'
#' Shows the legends for single-cell metadata based on the shinycell config 
#' data.table. This allows user to visualise the different metadata to be 
#' plotted and make any modifications if necessary. Note that the display name 
#' is shown here instead of the actual name. For more information regarding 
#' display name, see \code{?modMetaName}.
#'
#' @param scConf shinycell config data.table
#' @param fontSize font size of legends. Decrease it if you have too many 
#'   items to display
#'
#' @return gtable plot showing the legends for different metadata
#'
#' @author John F. Ouyang
#'
#' @import data.table ggplot2 RColorBrewer grid gridExtra
#'
#' @examples
#' showLegend(scConf)
#' 
#' # Can also save the legend for plotting later
#' scLegend = showLegend(scConf)
#' leg = do.call(gtable_rbind, scLegend)
#' grid.newpage()
#' grid.draw(leg)
#'
#' @export
showLegend <- function(scConf, fontSize = 14){
  
  # Start making config data.table
  scLegend = list()
  for(iMeta in scConf$ID){
    if(!is.na(scConf[ID == iMeta]$fID)){
      # Ensure categorical
      ggData = data.table(X = 1, Y = 1,
                          col = strsplit(scConf[ID == iMeta]$fID, "\\|")[[1]])
      ggOut = ggplot(ggData, aes(X, Y, color = col)) + 
        geom_point(size = fontSize / 5) +
        scale_color_manual(scConf[ID == iMeta]$UI,
                           values = strsplit(scConf[ID == iMeta]$fCL, "\\|")[[1]],
                           labels = strsplit(scConf[ID == iMeta]$fUI, "\\|")[[1]]) + 
        theme_classic(base_size = fontSize) + theme(legend.position = "bottom") + 
        guides(color = guide_legend(nrow = scConf[ID == iMeta]$fRow))

      # Get legends
      tmp = ggplot_gtable(ggplot_build(ggOut)) 
      ggLeg = which(sapply(tmp$grobs, function(y) y$name) == "guide-box") 
      ggLeg = tmp$grobs[[ggLeg]]
      scLegend[[iMeta]] = ggLeg
    } 
  }
  
  # Make legend for continuous
  ggData = data.table(X = 1, Y = 1, col = 1:100)
  ggOut = ggplot(ggData, aes(X, Y, color = col)) + 
    geom_point(size = fontSize / 5) +
    scale_color_gradientn(paste0(scConf[is.na(fID)]$UI, collapse = "\n"),
                          colours = c("grey85", brewer.pal(9, "OrRd")),
                          breaks = c(1,100), labels = c("low","high")) +
    theme_classic(base_size = fontSize) + theme(legend.position = "bottom") +
    guides(color = guide_colorbar(barwidth = 15, barheight = 1))
  tmp = ggplot_gtable(ggplot_build(ggOut)) 
  ggLeg = which(sapply(tmp$grobs, function(y) y$name) == "guide-box") 
  ggLeg = tmp$grobs[[ggLeg]]
  scLegend[["continuous"]] = ggLeg

  # Plot all legends
  legMulti = do.call(gtable_rbind, scLegend)
  grid.newpage()
  grid.draw(legMulti)

  return(scLegend)
}


