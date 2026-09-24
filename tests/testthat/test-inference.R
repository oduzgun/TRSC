cert_fixture <- function(q=1,components=c(estimation=.05,fitting_gap=.02,causal_gap=0,numerical=0)) {
  g<-trsc_geometry(matrix(c(0,1),1,dimnames=list(NULL,c("A","B"))))
  ac<-trsc_anchor_certificate(g,c(A=.5,B=.5),q,components,.02,"maintained causal-anchor bounds",
       previous=if(q==2)c(A=.6,B=.4) else NULL,W=matrix(c(2,.5,.5,1),2))
  gc<-trsc_geometry_certificate(g,.01,.01,"maintained absolute chord certificate")
  list(g=g,ac=ac,gc=gc,X=rbind(c(A=0,B=1),c(A=0,B=1)))
}

test_that("Theorem 11.1 retains good-event witnesses with explicit allowances", {
  f<-cert_fixture(1)
  cs<-trsc_confidence_set(f$X,f$ac,f$gc,q=1,R=.2)
  pop<-trsc_counterfactual_set(f$g,f$X,c(A=.45,B=.55),1,.2,type="population")
  expect_equal(cs$energy_allowance,sqrt(2)*.01,tolerance=1e-14)
  expect_gte(trsc_confidence_support(cs,c(1,1))$upper,trsc_support(pop,c(1,1))$upper-1e-7)
  expect_equal(cs$coverage$failure_probability,.03)
  b<-trsc_bounds(cs);expect_true(all(b$lower<=b$upper))
  expect_true(trsc_confidence_support(cs,c(1,1))$certificate$valid)
})

test_that("q2 keeps joint level-velocity information on one event", {
  f<-cert_fixture(2)
  cs<-trsc_confidence_set(f$X,f$ac,f$gc,q=2,R=.1)
  expect_equal(cs$energy_allowance,2*sqrt(2)*.01,tolerance=1e-14)
  expect_equal(cs$coverage$failure_probability,.03)
  expect_true(any(f$ac$W[row(f$ac$W)!=col(f$ac$W)]!=0))
  expect_equal(trsc_anchor_propagation(f$ac,1:3)$level_plus_velocity,
    f$ac$radius0+(1:3)*f$ac$radius_velocity,tolerance=1e-12)
  expect_true(is.finite(trsc_confidence_support(cs,c(1,1))$upper))
})

test_that("confidence inference rejects missing, diagnostic and conditional certificates", {
  f<-cert_fixture()
  expect_error(trsc_confidence_set(f$X,NULL,f$gc,1,.1),"anchor_certificate")
  expect_error(trsc_confidence_set(f$X,f$ac,f$gc,1,.1,conditional=TRUE),"conditional")
  d<-trsc_geometry_certificate(f$g,.1,.01,"rank is only diagnostic",provenance="diagnostic_only")
  expect_error(trsc_confidence_set(f$X,f$ac,d,1,.1),"diagnostic")
  g2<-trsc_geometry(matrix(c(0,1,2),1,dimnames=list(NULL,c("A","B","C"))))
  d2<-trsc_geometry_certificate(g2,.1,.01,"other donor geometry")
  expect_error(trsc_confidence_set(f$X,f$ac,d2,1,.1),"geometry")
})

test_that("common-library events are counted once and no q is selected causally", {
  f<-cert_fixture();a<-f$ac;b<-f$ac
  a$candidate_id<-"m1";b$candidate_id<-"m2"
  a$covered_candidates<-b$covered_candidates<-c("m1","m2")
  cs<-trsc_confidence_set(f$X,a,f$gc,1,.2,library=list(a,b))
  expect_equal(cs$coverage$failure_probability,.03)
  expect_false(cs$coverage$structural_order_validated)
  b$event_id<-"second-path-event"
  expect_equal(trsc_confidence_set(f$X,a,f$gc,1,.2,library=list(a,b))$coverage$failure_probability,.05)
})

test_that("mean-effect translation has correct signs and probability accounting", {
  f<-cert_fixture();z<-trsc_confidence_set(f$X,f$ac,f$gc,1,.1)
  e<-trsc_effect_interval(z,c(1,1),treated_mean=3)
  expect_equal(e$interval,c(lower=3-trsc_support(z,c(1,1))$upper,upper=3+trsc_support(z,c(-1,-1))$upper))
  ee<-trsc_effect_interval(z,c(1,1),treated_mean_interval=c(2,4),mean_failure_probability=.01)
  expect_equal(ee$failure_probability,.04)
  expect_error(trsc_effect_interval(z,c(1,1),treated_mean=3,mean_type="realized"),"conditional mean")
})

test_that("deployment modulus sees exact and near-null ambiguity", {
  g<-trsc_geometry(matrix(c(0,0),1,dimnames=list(NULL,c("A","B"))))
  expect_equal(trsc_deployment_modulus(g,c(A=0,B=2),0)$upper,2,tolerance=1e-7)
  g<-trsc_geometry(matrix(c(0,1e-5),1,dimnames=list(NULL,c("A","B"))))
  expect_equal(trsc_deployment_modulus(g,c(A=0,B=1),1e-5)$upper,1,tolerance=1e-7)
})
