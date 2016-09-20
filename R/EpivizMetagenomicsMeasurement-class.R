#' Metagenomic Measurement object
#' 
#' @importClassesFrom epivizrData EpivizMeasurement
#' @docType class
#' 
#' @export

EpivizMetagenomicsMeasurement <- setClass("EpivizMetagenomicsMeasurement",
                                      contains="EpivizMeasurement"
)

setMethod("as.list", signature(x="EpivizMetagenomicsMeasurement"),
          function(x) {
            nms <- slotNames("EpivizMetagenomicsMeasurement")
            out <- lapply(nms, function(slot_name) slot(x, slot_name))
            names(out) <- nms
            length(out)
          }          
)