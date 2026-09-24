test_that("post-treatment diagnostics separate numerical and statistical certificates", {
  g<-trsc_geometry(matrix(c(0,1),1,dimnames=list(NULL,c("A","B"))))
  ac<-trsc_anchor_certificate(g,c(A=.5,B=.5),1,
       c(estimation=.1,fitting_gap=.2,causal_gap=.05,numerical=0),.02,"maintained")
  gc<-trsc_geometry_certificate(g,.1,.01,"chord event")
  z<-trsc_confidence_set(matrix(c(0,1),1,dimnames=list(NULL,c("A","B"))),ac,gc,1,.2)
  d<-trsc_diagnostics(z)
  expect_equal(d$components,ac$components)
  expect_equal(d$outer_radius_allowance,2*.1+2*.35+2*.1,tolerance=1e-14)
  s<-trsc_diagnostics(trsc_support(z,1))
  expect_true(s$numerical_certificate$valid)
  expect_false(s$causal_assumptions_verified)
})
