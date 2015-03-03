MetavizDeviceManager <- setRefClass("MetavizDeviceManager",
  fields=list(
    .epivizMgr="ANY" # EpivizDeviceManager
  ),
  methods=list(
    initialize=function(epivizMgr, ...) {
      .epivizMgr <<- epivizMgr
    },

    epivizDeviceManager=function() { .epivizMgr },



    # Wrap
#     show=function(...) { .epivizMgr$show(...) },
#     registerType=function(...) { .epivizMgr$registerType(...) },
#     registerAction=function(...) { .epivizMgr$registerAction(...) },
#     getSeqInfos=function(...) { .epivizMgr$getSeqInfos() },
#     standalone=function(...) { .epivizMgr$standalone() },
#     daemonized=function(...) { .epivizMgr$daemonized() },
#     waitToClearRequests=function(...) { .epivizMgr$waitToClearRequests(...) },
    addMeasurements=function(...) { .epivizMgr$addMeasurements(...) },
#     clearDatasourceGroupCache=function(...) { .epivizMgr$.clearDatasourceGroupCache(...) },
#     updateMeasurements=function(...) {},
#     getMsObject=function(...) { .epivizMgr$.getMsObject(msObjId) },
#     rmMeasurements=function(...) {},
#     rmAllMeasurements=function(...) {},
#     listMeasurements=function(...) {},
#     getMeasurements=function(...) {},
#     getMeasurementType=function(...) {},
    .findDatasource=function(...) { .epivizMgr$.findDataSource(...) },
#     getRows=function(...) {},
#     getValues=function(chr, start, end, datasource, measurement) {},
#     addChart=function(...) {},
#     .getChartObject=function(...) {},
#     rmChart=function(...) {},
#     rmAllCharts=function(...) {},
#     listCharts=function(...) {},
#     addDevice=function(...) {},
#     rmDevice=function(...) {},
#     rmAllDevices=function(...) {},
#     clearDeviceList=function(...) {},
#     updateDevice=function(...) {},
#     listDevices=function(...) {},
#     bindToServer=function(...) {},
#     isClosed=function(...) {},
#     openBrowser=function(...) {},
    service=function(...) { .epivizMgr$service(...) },
#     stopService=function(...) {},
    startServer=function(...) { .epivizMgr$startServer(...) },
    stopServer=function(...) { .epivizMgr$stopServer(...) }
#     refresh=function(...) {},
#     navigate=function(...) {},
#     getCurrentLocation=function(...) {},
#     slideshow=function(...) {},
#     addSeqinfo=function(...) {},
#     rmSeqinfo=function(...) {},
#     handle=function(...) {}
#
#     blockChart=function(...) {},
#     lineChart=function(...) {},
#     scatterChart=function(...) {},
#     heatmapChart=function(...) {},
#     genesChart=function(...) {},

  )
)
