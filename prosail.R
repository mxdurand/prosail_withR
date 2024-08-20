### Required functions

# jfunc
jfunc1 <- function(k, l, t) {
  del <- (k - l) * t
  out <- numeric(length(l))
  iii <- abs(del) > 1e-3
  out[iii] <- (exp(-l[iii] * t) - exp(-k * t)) / (k - l[iii])
  out[!iii] <- 0.5 * t * (exp(-k * t) + exp(-l[!iii] * t)) *
    (1 - (del[!iii]^2) / 12)
  out
}

jfunc2 <- function(k, l, t) {
  (1 - exp(-(k + l) * t)) / (k + l)
}

# Volume scattering function
volscatt <- function(solar_zenith, instrument_zenith, azimuth, leaf_angle)
{
  d2r <- pi / 180
  
  # Sun-sensor geoms
  costs <- cos(d2r * solar_zenith)
  sints <- sin(d2r * solar_zenith)
  costo <- cos(d2r * instrument_zenith)
  sinto <- sin(d2r * instrument_zenith)
  cospsi <- cos(d2r * azimuth)
  psirad <- d2r * azimuth
  
  # Leaf geoms
  costl <- cos(d2r * leaf_angle)
  sintl <- sin(d2r * leaf_angle)
  
  cs <- costl * costs
  co <- costl * costo
  ss <- sintl * sints
  so <- sintl * sinto
  
  # Transition angles (beta) for solar (bts) and instrument (bto) directions
  cosbts <- rep_len(5, length(cs))
  i <- abs(ss) > 1e-6
  if (any(i)) cosbts[i] <- -cs[i] / ss[i]
  
  cosbto <- rep_len(5, length(co))
  i <- abs(so) > 1e-6
  if (any(i)) cosbto[i] <- -co[i] / so[i]
  
  bts <- rep_len(pi, length(cosbts))
  ds <- rep_len(cs, length(cosbts))
  i <- abs(cosbts) < 1
  if (any(i)) {
    bts[i] <- acos(cosbts[i])
    ds[i] <- ss[i]
  }
  
  chi_s <- 2 / pi * ((bts - 0.5 * pi) * cs + sin(bts) * ss)
  
  bto <- rep_len(0, length(cosbto))
  doo <- rep_len(-co, length(cosbto))
  i <- abs(cosbto) < 1
  if (any(i)) {
    bto[i] <- acos(cosbto[i])
    doo[i] <- so[i]
  }
  j <- !i & instrument_zenith < 90
  if (any(j)) {
    bto[j] <- pi
    doo[j] <- co[j]
  }
  
  chi_o <- 2 / pi * ((bto - 0.5 * pi) * co + sin(bto) * so)
  
  # Auxiliary azimuth angles bt1, bt2, bt3, for bidirectional scattering
  btran1 <- abs(bts - bto)
  btran2 <- pi - abs(bts + bto - pi)
  
  bt1 <- rep_len(psirad, length(btran1))
  bt2 <- btran1
  bt3 <- btran2
  
  i <- psirad >= btran1
  bt1[i] <- btran1[i]
  bt2[i] <- psirad
  bt3[i] <- btran2[i]
  j <- psirad >= btran2
  bt2[j] <- btran2[j]
  bt3[j] <- psirad
  
  t1 <- 2 * cs * co + ss * so * cospsi
  t2 <- rep_len(0, length(t1))
  i <- bt2 > 0
  t2[i] <- sin(bt2[i]) * (2 * ds[i] * doo[i] + ss[i] * so[i] * cos(bt1[i]) * cos(bt3[i]))
  
  denom <- 2 * pi^2
  frho <- ((pi - bt2) * t1 + t2) / denom
  ftau <- (-bt2 * t1 + t2) / denom
  frho[frho < 0] <- 0
  ftau[ftau < 0] <- 0
  
  return(list(chi_s = chi_s, chi_o = chi_o, frho = frho, ftau = ftau))
}


# Leaf angle distribution (translated from matlab)
dcum <- function(a, b, theta)
{
  d2r <- pi / 180
  thetarad <- theta * d2r
  x0 <- x <- 2 * thetarad
  t <- 1e-6
  dx_abs <- 1
  while (any(dx_abs > t)) {
    y <- a * sin(x) + 0.5 * b * sin(2 * x)
    dx <- 0.5 * (y - x + x0)
    x <- x + dx
    dx_abs <- abs(dx)
  }
  2 * (y + thetarad) / pi
}

campbell <- function(ala)
{
  tx1=c(10,20,30,40,50,60,70,80,82,84,86,88,90)
  tx2=c(0,10,20,30,40,50,60,70,80,82,84,86,88)
  
  litab = (tx2 + tx1) / 2
  n = length(litab)
  tl1 = tx1 * (pi / 180)
  tl2 = tx2 * (pi / 180)
  excent = exp(-1.6184e-5 * ala^3 + 2.1145e-3 * ala^2 - 1.2390e-1 * ala + 3.2491)
  sum0 = 0
  freq = c()
  
  for (i in 1:n)
  {
    x1  = excent / (sqrt(1 + excent^2 * tan(tl1[i])^2))
    x2  = excent / (sqrt(1 + excent^2 * tan(tl2[i])^2))
    if(excent == 1){
      freq[i] = abs(cos(tl1[i]) - cos(tl2[i]))
    } else {
      alpha = excent / sqrt(abs(1 - excent^2))
      alpha2 = alpha^2
      x12 = x1^2
      x22 = x2^2
      
      if(excent > 1){
        alpx1 = sqrt(alpha2 + x12)
        alpx2 = sqrt(alpha2 + x22)
        dum   = x1 * alpx1 + alpha2 * log(x1 + alpx1)
        freq[i] = abs(dum - (x2 * alpx2 + alpha2 * log(x2 + alpx2)))
      } else {
        almx1 = sqrt(alpha2 - x12)
        almx2 = sqrt(alpha2 - x22)
        dum   = x1 * almx1 + alpha2 * asin(x1 / alpha)
        freq[i] = abs(dum - (x2 * almx2 + alpha2 * asin(x2 / alpha)))
      }
    }
  }
  sum0 = sum(freq)
  freq0 = freq / sum0
  return(data.frame("litab" = litab, "freq0" = freq0))
}

# 4SAIL
sail4 <- function(leaf_refl, leaf_trans, soil_refl, LAI, hot_spot = 0, 
                  solar_zenith = 0, instrument_zenith = 0, azimuth = 0,
                  LIDFa = -0.35, LIDFb = -0.15, lidf_type = 1)
{
  stopifnot(LAI >= 0, all(leaf_refl >= 0), all(leaf_trans >= 0), all(soil_refl >= 0), length(leaf_trans) %in% c(1, length(leaf_refl)), length(soil_refl) %in% c(1, length(leaf_refl)))
  
  if (LAI == 0) {
    return(list(bhr = soil_refl, dhr = soil_refl, hdr = soil_refl, bdr = soil_refl))
  }
  
  d2r <- pi/180
  cts <- cos(d2r * solar_zenith)
  cto <- cos(d2r * instrument_zenith)
  ctscto <- cts * cto
  tants <- tan(d2r * solar_zenith)
  tanto <- tan(d2r * instrument_zenith)
  cospsi <- cos(d2r * azimuth)
  dso <- sqrt(tants^2 + tanto^2 - 2 * tants * tanto * cospsi)
  
  if(lidf_type == 1){
    litab <- c(seq(5, 75, 10), seq(81, 89, 2))
    dcum_in <- c(seq(10, 80, 10), seq(82, 88, 2))
    F_lidf <- c(dcum(LIDFa, LIDFb, dcum_in), 1)
    lidf <- F_lidf
    for (i in seq(length(litab), 2)) {
      lidf[i] <- lidf[i] - lidf[i - 1]
    }
  } else if(lidf_type == 2){
    LAD <- campbell(LIDFa)
    lidf <- LAD$freq0
    litab <- LAD$litab
  }
  
  
  ctl <- cos(d2r * litab)
  vs <- volscatt(solar_zenith, instrument_zenith, azimuth, litab)
  
  ksli <- vs$chi_s/cts
  koli <- vs$chi_o/cto
  sobli <- vs$frho * pi/ctscto
  sofli <- vs$ftau * pi/ctscto
  bfli <- ctl * ctl
  
  ks <- sum(ksli * lidf)
  ko <- sum(koli * lidf)
  bf <- sum(bfli * lidf)
  sob <- sum(sobli * lidf)
  sof <- sum(sofli * lidf)
  
  sdb <- 0.5 * (ks + bf)
  sdf <- 0.5 * (ks - bf)
  dob <- 0.5 * (ko + bf)
  dof <- 0.5 * (ko - bf)
  ddb <- 0.5 * (1 + bf)
  ddf <- 0.5 * (1 - bf)
  
  sigb <- (ddb * leaf_refl) + (ddf * leaf_trans)
  sigf <- (ddf * leaf_refl) + (ddb * leaf_trans)
  att <- 1 - sigf
  
  m2 <- (att + sigb) * (att - sigb)
  m2[m2 < 0] <- 0
  m <- sqrt(m2)
  
  sb <- sdb * leaf_refl + sdf * leaf_trans
  sf <- sdf * leaf_refl + sdb * leaf_trans
  vb <- dob * leaf_refl + dof * leaf_trans
  vf <- dof * leaf_refl + dob * leaf_trans
  w <- sob * leaf_refl + sof * leaf_trans
  
  e1 <- exp(-m * LAI)
  e2 <- e1 * e1
  
  rinf <- (att - m)/sigb
  rinf2 <- rinf * rinf
  re <- rinf * e1
  denom <- 1 - rinf2 * e2
  
  J1ks <- jfunc1(ks, m, LAI)
  J2ks <- jfunc2(ks, m, LAI)
  J1ko <- jfunc1(ko, m, LAI)
  J2ko <- jfunc2(ko, m, LAI)
  
  Ps <- (sf + sb * rinf) * J1ks
  Qs <- (sf * rinf + sb) * J2ks
  Pv <- (vf + vb * rinf) * J1ko
  Qv <- (vf * rinf + vb) * J2ko
  
  rdd <- rinf * (1 - e2)/denom
  tdd <- (1 - rinf2) * e1/denom
  tsd <- (Ps - re * Qs)/denom
  rsd <- (Qs - re * Ps)/denom
  tdo <- (Pv - re * Qv)/denom
  rdo <- (Qv - re * Pv)/denom
  
  tss <- exp(-ks * LAI)
  too <- exp(-ko * LAI)
  
  z <- jfunc2(ks, ko, LAI)
  
  g1 <- (z - J1ks * too)/(ko + m)
  g2 <- (z - J1ko * tss)/(ks + m)
  
  Tv1 <- (vf * rinf + vb) * g1
  Tv2 <- (vf + vb * rinf) * g2
  
  T1 <- Tv1 * (sf + sb * rinf)
  T2 <- Tv2 * (sf * rinf + sb)
  T3 <- (rdo * Qs + tdo * Ps) * rinf
  rsod <- (T1 + T2 - T3)/(1 - rinf2)
  
  alf <- 1e+06
  if(hot_spot > 0) {alf <- (dso/hot_spot) * 2/(ks + ko)}
  if(alf > 200) {alf <- 200}
  if (alf == 0) {
    tstoo <- tss
    sumint <- (1 - tss)/(ks * LAI)
  } else {
    fhot <- LAI * sqrt(ks * ko)
    x1 = 0
    y1 = 0
    f1 = 1
    fint = (1 - exp(-alf)) * 0.05
    sumint = 0
    for (i in seq(1, 20))
    {
      if (i < 20) {
        x2 <- -log(1 - i * fint)/alf
      } else {
        x2 <- 1
      }
      y2 <- -(ko + ks) * LAI * x2 + fhot * (1 - exp(-alf * 
                                                      x2))/alf
      f2 <- exp(y2)
      sumint <- sumint + (f2 - f1) * (x2 - x1)/(y2 - y1)
      x1 <- x2
      y1 <- y2
      f1 <- f2
    }
    tsstoo <- f1
  }
  
  rsos <- w * LAI * sumint
  rso <- rsos + rsod
  dn <- 1 - soil_refl * rdd
  bhr <- rdd + tdd * soil_refl * tdd / dn
  dhr <- rsd + (tsd + tss) * soil_refl * tdd / dn
  hdr <- rdo + tdd * soil_refl * (tdo + too) / dn
  rsodt <- rsod + ((tss + tsd) * tdo + (tsd + tss * soil_refl * rdd) * too) * soil_refl / dn
  rsost <- rsos + tsstoo * soil_refl
  bdr <- rsost + rsodt
  
  return(list(bhr = bhr, dhr = dhr, hdr = hdr, bdr = bdr))
}

