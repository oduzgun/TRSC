stress_geometry <- function() trsc_geometry(matrix(c(0,1),1,dimnames=list(NULL,c("A","B"))))
stress_cert_fixture <- function() {
  g<-stress_geometry()
  ac<-trsc_anchor_certificate(g,c(A=.5,B=.5),2,
      c(estimation=.05,fitting_gap=.02,causal_gap=0,numerical=0),.02,"maintained bounds",
      previous=c(A=.6,B=.4),W=matrix(c(2,.5,.5,1),2))
  gc<-trsc_geometry_certificate(g,.01,.01,"maintained chord bound")
  list(g=g,ac=ac,gc=gc,X=rbind(c(A=0,B=1),c(A=0,B=1)))
}

test_that("S03/S04/S12/S14: singular J>T and noise-free nonunique fits are valid", {
  t<-1:4;units<-c("Target",paste0("D",1:6))
  p<-expand.grid(unit=units,period=t,stringsAsFactors=FALSE);p$y<-0
  d<-trsc_data(p,"unit","period","y","Target",units[-1],5,t,"y")
  g<-trsc_geometry(matrix(0,2,6,dimnames=list(NULL,units[-1])))
  fit<-trsc_fit(d,g,kappa=1,q=1)
  expect_s3_class(fit,"trsc_fit");expect_equal(fit$objective$total,0,tolerance=1e-10)
  expect_gt(fit$diagnostics$simplex_null_dimension_empirical,0)
  expect_equal(sum(stats::fitted(fit)^2),0,tolerance=1e-10)
})

test_that("S03/S05/S09/S15: forward-seminorm inference keeps a wide valid set", {
  g<-trsc_geometry(matrix(c(0,0,0),1,dimnames=list(NULL,c("A","B","C"))))
  a<-trsc_anchor_certificate(g,c(A=1,B=0,C=0),1,
       c(estimation=0,fitting_gap=0,causal_gap=0,numerical=0),0,"exact boundary")
  gc<-trsc_geometry_certificate(g,1,0,"declared population chord allowance")
  X<-matrix(c(-1e6,1e6,1e-6),1,dimnames=list(NULL,c("A","B","C")))
  z<-trsc_confidence_set(X,a,gc,1,1e8)
  expect_equal(trsc_support(z,1)$upper,1e6,tolerance=1e-3)
  expect_equal(-trsc_support(z,-1)$upper,-1e6,tolerance=1e-3)
})

test_that("S17: equivalent raw row order does not change pre-treatment fit", {
  panel<-expand.grid(unit=c("Target","A","B"),period=1:7,stringsAsFactors=FALSE)
  panel$outcome<-c(rbind(.4*(1:7)+.6*(7:1),1:7,7:1))
  make<-function(x)trsc_data(x,"unit","period","outcome","Target",c("A","B"),8,1:7,"outcome")
  p<-make(panel);pp<-make(panel[rev(seq_len(nrow(panel))),])
  g<-trsc_geometry(p$blocks$Y0)
  f<-trsc_fit(p,g,1,1);ff<-trsc_fit(pp,g,1,1)
  expect_equal(stats::fitted(f),stats::fitted(ff),tolerance=1e-10)
})

test_that("Lemma M.2 endpoint interpolation holds for vector-valued sequences", {
  set.seed(912)
  for(k in 2:20) {
    E<-matrix(rnorm(k*4),k,4)
    lhs<-sum((E[k,]-E[k-1,])^2)
    curv<-if(k>2)sum((.trsc_difference_matrix(k,2)%*%E)^2) else 0
    rhs<-32/k^3*sum(E^2)+2*k*curv
    expect_lte(lhs,rhs+1e-12)
  }
})

test_that("Proposition M.1 sandwich and q2 minimum-radius inclusion hold", {
  f<-stress_cert_fixture();f$gc$epsilon<-0
  cs<-trsc_confidence_set(f$X,f$ac,f$gc,2,.3)
  ac<-f$ac;t0<-2*ac$radius0;tv<-2*ac$radius_velocity
  omega<-sqrt((t0+tv)^2+t0^2)
  pop<-trsc_counterfactual_set(f$g,f$X,c(A=.5,B=.5),2,.3,previous=c(A=.6,B=.4))
  broad<-pop;broad$R<-.3+omega
  for(a in list(c(1,1),c(-1,1),c(0,-1))) {
    expect_gte(trsc_support(cs,a)$upper,trsc_support(pop,a)$upper-1e-6)
    expect_lte(trsc_support(cs,a)$upper,trsc_support(broad,a)$upper+1e-6)
  }
  ac<-trsc_anchor_certificate(f$g,c(A=0,B=1),2,
      c(estimation=0,fitting_gap=0,causal_gap=0,numerical=0),0,"exact population anchors",
      previous=c(A=1,B=0))
  exact<-trsc_confidence_set(f$X,ac,f$gc,2,1)
  expect_gte(trsc_support(exact,c(1,1))$upper,2)
})

test_that("risk-to-anchor Markov conversion uses one path allocation", {
  g<-stress_geometry();ell<-rbind(c(0,0,1),c(0,1,0))
  a<-trsc_anchor_from_risk(g,c(A=.5,B=.5),2,.01,.02,.005,ell,
       fitting_gap=.1,causal_gap=.02,numerical=0,source="known envelope and transfers",
       previous=c(A=.6,B=.4))
  expect_equal(a$failure_probability,.025)
  expect_equal(a$common_path_event$squared_path_radius,1.5)
  expect_true(a$common_path_event$level_and_velocity_share_event)
  L<-rbind(ell[1,],ell[1,]-ell[2,])
  expected<-sqrt(1.5)*svd(L,nu=0,nv=0)$d[1]
  expect_gte(a$components[["estimation"]],expected-1e-14)
  expect_equal(a$components[["estimation"]],expected,tolerance=1e-12)
})
