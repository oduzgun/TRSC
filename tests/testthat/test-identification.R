line_geometry <- function() trsc_geometry(matrix(c(0, 1), 1, dimnames = list(NULL, c("A", "B"))))

test_that("Theorem 8.1 retains evaluation fibers and compatible lifts", {
  g <- trsc_geometry(matrix(c(0, 0, 1), 1, dimnames = list(NULL, c("A", "B", "C"))))
  e <- trsc_evaluation_set(g, c(A = 0, B = 2, C = 7), representative = c(A=1,B=0,C=0))
  expect_equal(e$interval, c(lower=0, upper=2), tolerance=1e-7)
  expect_equal(trsc_evaluation_width(e), 2, tolerance=1e-7)
  expect_false(trsc_point_identified(e))
  singleton <- trsc_evaluation_set(g, c(A=1,B=1,C=3), representative=c(A=1,B=0,C=0))
  expect_true(trsc_point_identified(singleton, tolerance=1e-6))
  expect_length(e$compatible_lifts, 2)
})

test_that("q1 persistence fixes profiles but leaves null evaluation freedom", {
  g <- trsc_geometry(matrix(c(0,0,1),1,dimnames=list(NULL,c("A","B","C"))))
  z <- trsc_counterfactual_set(g, rbind(c(A=0,B=2,C=5),c(A=1,B=3,C=7)),
                                anchor=c(A=1,B=0,C=0),q=1,R=0)
  expect_equal(trsc_support(z,c(1,1))$upper,5,tolerance=1e-6)
  expect_equal(trsc_support(z,c(-1,-1))$upper,-1,tolerance=1e-6)
})

test_that("q2 affine continuation and minimum curvature reproduce analytic cases", {
  g <- line_geometry(); X <- matrix(rep(c(0,1),3),3,byrow=TRUE,dimnames=list(NULL,c("A","B")))
  z <- trsc_counterfactual_set(g,X,anchor=c(A=.6,B=.4),previous=c(A=.7,B=.3),q=2,R=0)
  expect_equal(trsc_support(z,c(0,0,1))$upper,.7,tolerance=2e-6)
  r <- trsc_min_radius(g,anchor=c(A=0,B=1),previous=c(A=1,B=0),q=2,H=2)
  expect_lte(r$lower,1);expect_gte(r$upper,1)
  expect_equal(r$estimate,1,tolerance=1e-6)
  at <- trsc_counterfactual_set(g,X[1:2,,drop=FALSE],c(A=0,B=1),2,1,previous=c(A=1,B=0))
  expect_equal(trsc_support(at,c(1,1))$upper,2,tolerance=1e-6)
  below <- trsc_counterfactual_set(g,X[1:2,,drop=FALSE],c(A=0,B=1),2,.5,previous=c(A=1,B=0))
  expect_identical(trsc_support(below,c(1,1))$status,"structural_infeasible")
})

test_that("dynamic supports are monotone and bounded by donor-hull supports", {
  g <- line_geometry(); X <- rbind(c(A=0,B=1),c(A=0,B=1))
  z <- trsc_counterfactual_set(g,X,c(A=1,B=0),1,.1)
  curve <- trsc_sensitivity(z,c(0,.1,.5,2),a=c(1,1))
  expect_true(all(diff(curve$upper)>=-1e-7))
  expect_lte(max(curve$upper),2)
  expect_lt(curve$upper[2],2)
})

test_that("auxiliary refinement uses one future witness before projection", {
  g <- trsc_geometry(matrix(c(0,0),1,dimnames=list(NULL,c("A","B"))))
  z <- trsc_counterfactual_set(g,matrix(c(0,2),1,dimnames=list(NULL,c("A","B"))),c(A=1,B=0),1,0)
  aux <- trsc_auxiliary_refine(z,list(matrix(c(0,1),1)),list(.25),source="untreated composition measurement")
  expect_equal(trsc_support(aux,1)$upper,.5,tolerance=1e-6)
  expect_equal(-trsc_support(aux,-1)$upper,.5,tolerance=1e-6)
  bad <- trsc_auxiliary_refine(z,list(matrix(c(0,1),1)),list(2),source="incompatible exact measurement")
  expect_identical(trsc_support(bad,1)$status,"structural_infeasible")
})

test_that("support agrees with independent enumeration and donor permutation", {
  g<-line_geometry();X<-rbind(c(A=0,B=1),c(A=1,B=-1))
  z<-trsc_counterfactual_set(g,X,c(A=1,B=0),1,.7)
  grid<-expand.grid(b1=seq(0,1,by=.002),b2=seq(0,1,by=.002))
  ok<-grid$b1^2+(grid$b2-grid$b1)^2<=.7^2
  target<-grid$b1+1-2*grid$b2
  s<-trsc_support(z,c(1,1))
  expect_gte(s$upper,max(target[ok]))
  expect_lt(s$upper-max(target[ok]),.01)
  gp<-trsc_geometry(matrix(c(1,0),1,dimnames=list(NULL,c("B","A"))))
  zp<-trsc_counterfactual_set(gp,X[,c("B","A")],c(A=1,B=0),1,.7)
  expect_equal(trsc_support(zp,c(1,1))$upper,s$upper,tolerance=1e-6)
})

test_that("duplicating all augmented donor characteristics preserves supports", {
  g<-line_geometry();X<-rbind(c(A=0,B=1),c(A=0,B=2))
  z<-trsc_counterfactual_set(g,X,c(A=1,B=0),1,.4)
  gd<-trsc_geometry(matrix(c(0,1,0),1,dimnames=list(NULL,c("A","B","A2"))))
  zd<-trsc_counterfactual_set(gd,cbind(X,A2=X[,"A"]),c(A=1,B=0,A2=0),1,.4)
  expect_equal(trsc_support(z,c(1,1))$upper,trsc_support(zd,c(1,1))$upper,tolerance=2e-6)
})
