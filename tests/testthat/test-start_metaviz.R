context("create MetavizApp class")

test_that("startMetaviz creates a MetavizApp Object", {
  app <- startMetaviz(non_interactive=TRUE)
  expect_is(app, "MetavizApp")
  
  expect_is(app$server, "EpivizServer")
  expect_is(app$chart_mgr, "EpivizChartMgr")
  expect_is(app$data_mgr, "EpivizDataMgr")
  
  expect_true(app$server$is_closed())
})
