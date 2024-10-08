! ----------------------------------------------------------------------------
SUBROUTINE g_est(y_rank, z_sort, tz, n, p, K, delta0, pro0, ret)
! ----------------------------------------------------------------------------
  IMPLICIT NONE
  ! - - - - input - - - - 
  INTEGER :: y_rank (n, p), z_sort (n, p), tz (K), n, p, K
  DOUBLE PRECISION :: delta0, pro0 (K)
  ! - - - - output - - - - 
  DOUBLE PRECISION :: ret (n, p)
  ! - - - - local arg - - - - 
  INTEGER :: i, j, l, ct (n, K)
  DOUBLE PRECISION :: tz_inv (K), prob1 (K, n, p)
  DOUBLE PRECISION :: one_minus_delta0, qval
  qval = 0.0D0
  one_minus_delta0 = 1.0D0 - delta0
  DO i = 1, K
    IF (tz(i) > 0) THEN
      tz_inv(i) = 1.0D0 / tz(i)
    ELSE
      tz_inv(i) = 0.0D0
    ENDIF
  ENDDO
  DO j = 1, p
    ct(1, :) = 0
    ct(1, z_sort(1, j)) = 1
    DO i = 2, n
      ct(i, :) = ct(i-1, :)
      ct(i, z_sort(i, j)) = ct(i, z_sort(i, j)) + 1
    ENDDO
    ! CALL INTPR("ct", -1, ct, n*K)
    ! CALL DBLEPR("prod", -1, ct(y_rank(i, j), :) * tz_inv, K)
    DO i = 1, n
      prob1(:, i, j) = ct(y_rank(i, j), :) * tz_inv
    ENDDO
  ENDDO
  DO j = 1, p
    DO i = 1, n
      DO l = 1, K
        IF (prob1(l, i, j) .LE. delta0) THEN
          prob1(l, i, j) = delta0
        ELSEIF (prob1(l, i, j) > one_minus_delta0) THEN
          prob1(l, i, j) = one_minus_delta0
        ENDIF
        CALL qnorm(prob1(l, i, j), qval)
        prob1(l, i, j) = qval
      ENDDO
      ret(i, j) = Dot_product(prob1(:, i, j), pro0)
    ENDDO
  ENDDO
END SUBROUTINE g_est

! ----------------------------------------------------------------------------
SUBROUTINE alpha_fn(g, mu, pro, gam, n, p, K, ret)
! ----------------------------------------------------------------------------
  IMPLICIT NONE
  ! - - - - input - - - - 
  INTEGER :: n, p, K
  DOUBLE PRECISION :: g (n, p), mu (K, p), pro (K), gam (K-1, p)
  ! - - - - output - - - - 
  DOUBLE PRECISION :: ret (n, K)
  ! - - - - local arg - - - - 
  INTEGER :: i, j
  DOUBLE PRECISION :: mu_sum (p, K-1), log_pro (K-1)
  log_pro = Log(pro(2:K)) 
  ret(:, 1) = - Log(pro(1))
  DO j = 1, K-1
    mu_sum(:, j) = (mu(j+1, :) + mu(1, :)) * 0.5D0
    DO i = 1, n
      ret(i, j+1) = Dot_product(gam(j, :), g(i, :) - mu_sum(:, j))
    ENDDO
  ENDDO
  ! CALL DBLEPR("mu1", -1, mu(1, :) , p)
  ! CALL DBLEPR("mu_sum1", -1, (mu(1+1, :) + mu(1, :)) * 0.5D0, p)
  ! CALL DBLEPR("mu_sum", -1, mu_sum, p*(K-1))
  DO i = 1, n
    ret(i, 2:K) = ret(i, 2:K) - log_pro
  ENDDO
END SUBROUTINE alpha_fn


! ----------------------------------------------------------------------------
SUBROUTINE z_fn(g, mu, pro, gam, n, p, K, ret)
! ----------------------------------------------------------------------------
  IMPLICIT NONE
  ! - - - - input - - - - 
  INTEGER :: n, p, K
  DOUBLE PRECISION :: g (n, p), mu (K, p), pro (K), gam (K-1, p)
  ! - - - - output - - - - 
  DOUBLE PRECISION :: ret (n, K)
  ! - - - - local arg - - - - 
  INTEGER :: i, j
  DOUBLE PRECISION :: mu_sum (p, K-1), log_pro (K-1)
  log_pro = Log(pro(2:K)) 
  ret(:, 1) = - Log(pro(1))
  DO j = 1, K-1
    mu_sum(:, j) = (mu(j+1, :) + mu(1, :)) * 0.5D0
    DO i = 1, n
      ret(i, j+1) = Dot_product(gam(j, :), g(i, :) - mu_sum(:, j))
    ENDDO
  ENDDO
  ! CALL DBLEPR("mu1", -1, mu(1, :) , p)
  ! CALL DBLEPR("mu_sum1", -1, (mu(1+1, :) + mu(1, :)) * 0.5D0, p)
  ! CALL DBLEPR("mu_sum", -1, mu_sum, p*(K-1))
  DO i = 1, n
    ret(i, 2:K) = ret(i, 2:K) - log_pro
  ENDDO
END SUBROUTINE z_fn

! ----------------------------------------------------------------------------
SUBROUTINE sigma_est(mu, g, alpha, n, p, K, ninv, ret)
! ----------------------------------------------------------------------------
  IMPLICIT NONE
  ! - - - - input - - - - 
  INTEGER :: n, p, K
  DOUBLE PRECISION :: mu (K, p), g (n, p), alpha (n, K), ninv
  ! - - - - output - - - - 
  DOUBLE PRECISION :: ret (p, p)
  ! - - - - local arg - - - - 
  INTEGER :: i, j, ii, jj
  DOUBLE PRECISION :: sig_tmp (p, p), dif (p)
  ret = 0.0D0
  sig_tmp = 0.0D0

  DO j = 1, K
    DO i = 1, n
      dif = g(i, :) - mu(j, :)
      DO jj = 1, p
        DO ii = 1, p
          IF (ii .GE. jj) THEN
            sig_tmp(ii, jj) = dif(ii) * dif(jj) 
          ELSE
            sig_tmp(ii, jj) = sig_tmp(jj, ii)
          ENDIF
        ENDDO
      ENDDO
      ret = ret + alpha(i, j) * ninv * sig_tmp
    ENDDO
  ENDDO
END SUBROUTINE sigma_est

! ! ----------------------------------------------------------------------------
! SUBROUTINE sigma_est(mu, g, alpha, n, p, K, ninv, ret)
! ! ----------------------------------------------------------------------------
!   IMPLICIT NONE
!   ! - - - - input - - - - 
!   INTEGER :: n, p, K
!   DOUBLE PRECISION :: mu (K, p), g (n, p), alpha (n, K), ninv
!   ! - - - - output - - - - 
!   DOUBLE PRECISION :: ret (p, p)
!   ! - - - - local arg - - - - 
!   INTEGER :: i, j, ii, jj
!   DOUBLE PRECISION :: sig_tmp (p, p), dif (p)
!   ret = 0.0D0
!   sig_tmp = 0.0D0

!   DO j = 1, K
!     DO i = 1, n
!       dif = g(i, :) - mu(j, :)
!       DO jj = 1, p
!         DO ii = 1, p
!           IF (ii .GE. jj) THEN
!             sig_tmp(ii, jj) = dif(ii) * dif(jj) 
!           ELSE
!             sig_tmp(ii, jj) = sig_tmp(jj, ii)
!           ENDIF
!         ENDDO
!       ENDDO
!       ret = ret + alpha(i, j) * ninv * sig_tmp
!     ENDDO
!   ENDDO
! END SUBROUTINE sigma_est


! ----------------------------------------------------------------------------
SUBROUTINE loglik(y, n, p, K, mu, ginv, det_val, ret)
! ----------------------------------------------------------------------------
  IMPLICIT NONE
  ! - - - - input - - - - 
  INTEGER :: n, p, K
  DOUBLE PRECISION :: y (n, p), mu (K, p), ginv (p, p), det_val
  ! - - - - output - - - - 
  DOUBLE PRECISION :: ret (n, K)
  ! - - - - local arg - - - - 
  INTEGER :: j, i
  DOUBLE PRECISION :: y_minus_mu (p)
  DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932D0
  
  DO i = 1, K
    DO j = 1, n
      y_minus_mu = y(j, :) - mu(i, :)
      ret(j, i) = - Dot_product(Matmul(ginv, y_minus_mu), y_minus_mu) * 0.5D0
    ENDDO
  ENDDO
  ret = ret - Log(2 * pi) * p * 0.5D0 - det_val
END SUBROUTINE loglik

! ----------------------------------------------------------------------------
SUBROUTINE qnorm(pval, qval)
! ----------------------------------------------------------------------------
  IMPLICIT NONE
  DOUBLE PRECISION :: pval, qval, r, qtmp
  qtmp = pval - 0.5D0
  IF (Abs(qtmp) .LE. 0.425) THEN
    r = .180625 - qtmp * qtmp
    qval = qtmp * (((((((r * 2509.0809287301226727 + &
      & 33430.575583588128105) * r + 67265.770927008700853) * r + &
      &  45921.953931549871457) * r + 13731.693765509461125) * r + &
      & 1971.5909503065514427) * r + 133.14166789178437745) * r + &
      & 3.387132872796366608) &
      &  / (((((((r * 5226.495278852854561 + &
      & 28729.085735721942674) * r + 39307.89580009271061) * r + &
      & 21213.794301586595867) * r + 5394.1960214247511077) * r + &
      & 687.1870074920579083) * r + 42.313330701600911252) * r + 1.0D0)
  ELSE 
    IF (qtmp > 0) THEN
      r = 1.0D0 - pval
    ELSE 
      r = pval
    ENDIF
    r = Sqrt(-Log(r))

    IF (r .LE. 5.0D0) THEN
      r = r - 1.6D0
      qval = (((((((r * 7.7454501427834140764D-4 + &
        & 0.0227238449892691845833) * r + 0.24178072517745061177) * &
        & r + 1.27045825245236838258) * r + &
        & 3.64784832476320460504) * r + 5.7694972214606914055) * &
        & r + 4.6303378461565452959) * r + &
        & 1.42343711074968357734) &
        & / (((((((r * &
        & 1.05075007164441684324D-9 + 5.475938084995344946D-4) * &
        & r + 0.0151986665636164571966) * r + &
        & 0.14810397642748007459) * r + 0.68976733498510000455) * &
        & r + 1.6763848301838038494) * r + &
        & 2.05319162663775882187) * r + 1.0D0)
    ELSE 
      r = r - 5.0D0
      qval = (((((((r * 2.01033439929228813265D-7 + &
        & 2.71155556874348757815D-5) * r + &
        & 0.0012426609473880784386) * r + 0.026532189526576123093) * &
        & r + 0.29656057182850489123) * r + &
        & 1.7848265399172913358) * r + 5.4637849111641143699) * &
        & r + 6.6579046435011037772) &
        & / (((((((r * &
        & 2.04426310338993978564D-15 + 1.4215117583164458887D-7)* &
        & r + 1.8463183175100546818D-5) * r + &
        & 7.868691311456132591D-4) * r + 0.0148753612908506148525) &
        & * r + 0.13692988092273580531) * r + &
        & 0.59983220655588793769) * r + 1.0D0)
    ENDIF
    IF (qtmp < 0.0D0) THEN
      qval = -qval
    ENDIF
  ENDIF

END SUBROUTINE qnorm



! copied from msda package
! --------------------------------------------------
SUBROUTINE msda(obj,nk,nvars,sigma,delta,pf,dfmax,pmax,nlam,flmin,ulam,&
  eps,maxit,sml,verbose,nalam,theta,m,ntheta,alam,npass,jerr)
! --------------------------------------------------
    IMPLICIT NONE
    ! - - - arg types - - -
    DOUBLE PRECISION, PARAMETER :: big=9.9E30
    DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
    INTEGER, PARAMETER :: mnlam = 6
    INTEGER::mnl
    INTEGER::nk
    INTEGER::nvars
    INTEGER::dfmax
    INTEGER::pmax
    INTEGER::nlam
    INTEGER::nalam
    INTEGER::npass
    INTEGER::jerr
    INTEGER::verbose
    INTEGER::maxit
    INTEGER::m(pmax)
    INTEGER::ntheta(nlam)
    DOUBLE PRECISION::flmin
    DOUBLE PRECISION::eps
    DOUBLE PRECISION::sml
    DOUBLE PRECISION::sigma(nvars,nvars)
    DOUBLE PRECISION::delta(nk,nvars)
    DOUBLE PRECISION::pf(nvars)
    DOUBLE PRECISION::ulam(nlam)
    DOUBLE PRECISION::theta(nk,pmax,nlam)
    DOUBLE PRECISION::alam(nlam)
    DOUBLE PRECISION::obj(nlam)
    ! - - - local declarations - - -
    INTEGER::mm(nvars)
    INTEGER::k
    INTEGER::j
    INTEGER::jj
    INTEGER::l
    INTEGER::vrg
    INTEGER::ni
    INTEGER::me
    DOUBLE PRECISION::dif
    DOUBLE PRECISION::v
    DOUBLE PRECISION::al
    DOUBLE PRECISION::alf
    DOUBLE PRECISION::unorm
    DOUBLE PRECISION::thetanew(nk,nvars)
    DOUBLE PRECISION::thetaold(nk,nvars)
    DOUBLE PRECISION::r(nk,nvars)
    DOUBLE PRECISION::ab(nk,nvars)
    DOUBLE PRECISION::d(nk)
    DOUBLE PRECISION::theta_sum(nk)
    DOUBLE PRECISION::thetatmp(nk)
    DOUBLE PRECISION::u(nk)
    DOUBLE PRECISION::loss_diff
    DOUBLE PRECISION::penalty_diff
    DOUBLE PRECISION::dev
    DOUBLE PRECISION::dev_tmp
    DOUBLE PRECISION::tmp1_new
    DOUBLE PRECISION::tmp2_new
    DOUBLE PRECISION::dev_new
    DOUBLE PRECISION::dev1_new
    DOUBLE PRECISION::dev2_new
    DOUBLE PRECISION::dev3_new
    DOUBLE PRECISION::tmp1_old
    DOUBLE PRECISION::tmp2_old
    DOUBLE PRECISION::dev_old
    DOUBLE PRECISION::dev1_old
    DOUBLE PRECISION::dev2_old
    DOUBLE PRECISION::dev3_old
! - - - begin - - -
    IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
    ENDIF
    pf=max(0.0D0,pf)
! - - - some initial setup - - -
    mnl = Min (mnlam, nlam)
    r = delta
    thetanew=0.0D0
    thetaold=0.0D0
    dev=0.0D0
    m=0
    mm=0
    npass=0
    ni=npass
! --------- lambda loop ----------------------------
    IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf=flmin**(1.0D0/(nlam-1.0D0))
    ENDIF
    DO l=1,nlam
        IF(flmin>=1.0D0) THEN
            al=ulam(l)
        ELSE
            IF(l > 2) THEN
                al=al*alf
            ELSE IF(l==1) THEN
                al=big
            ELSE IF(l==2) THEN
                al=0.0D0
                DO j = 1,nvars
                    IF(pf(j)>0.0D0) THEN
                            u = delta(:,j)
                            v = sqrt(dot_product(u,u))
                            al=max(al,v/pf(j))
                    ENDIF
                END DO
                al=al*alf
            ENDIF
        ENDIF
! --------- outer loop ----------------------------
        DO
            IF(ni>0) thetaold(:,m(1:ni))=thetanew(:,m(1:ni))
! --middle loop-------------------------------------
            DO
                npass=npass+1
                dif=0.0D0
                dev_tmp = dev
                DO k=1,nvars
                    thetatmp=thetanew(:,k)
                    u = r(:,k)/sigma(k,k) + thetatmp
                    unorm = sqrt(dot_product(u,u))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! check this?                    
                    v = unorm-al*pf(k)/sigma(k,k)
                    IF(v > 0.0D0) THEN
                        thetanew(:,k) = v*u/unorm
                    ELSE
                        thetanew(:,k)=0.0D0
                    ENDIF
                    d=thetanew(:,k)-thetatmp
                    ! CALL DBLEPR("d", -1, maxval(abs(d)), 1)
                    theta_sum=thetanew(:,k)+thetatmp
                    IF(any(d/=0.0D0)) THEN
                        dif=max(dif,maxval(abs(d)))
                        loss_diff = sum(d*(0.5*theta_sum-r(:,k)-sigma(k,k)*thetatmp))
                        penalty_diff = al*pf(k)*(sqrt(dot_product(thetanew(:,k),thetanew(:,k))) &
                        - sqrt(dot_product(thetatmp,thetatmp)))
                        dev = dev + loss_diff + penalty_diff
                        ab=spread(sigma(:,k),dim=1,ncopies=nk)*spread(d,dim=2,ncopies=nvars)
                        r=r-ab
                        IF(mm(k)==0) THEN
                            ni=ni+1
                            IF(ni>pmax) EXIT
                                mm(k)=ni
                                m(ni)=k
                            ENDIF
                        ENDIF
                ENDDO
                IF(abs((dev-dev_tmp)/dev)<sml) EXIT
                IF(ni>pmax) EXIT
                IF(dif<eps) EXIT
                 ! CALL DBLEPR("dif_middle", -1, dif, 1)
                IF(npass > maxit) THEN
                    jerr=-l
                    RETURN
               ENDIF
! --inner loop----------------------
                DO
                    npass=npass+1
                    dif=0.0D0
                    dev_tmp = dev
                    DO j=1,ni
                        k=m(j)
                        thetatmp=thetanew(:,k)
                        u = r(:,k)/sigma(k,k) + thetatmp
                        unorm = sqrt(dot_product(u,u))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! check this?
                        v = unorm-al*pf(k)/sigma(k,k)
                        IF(v > 0.0D0) THEN
                            thetanew(:,k) = v*u/unorm
                        ELSE
                            thetanew(:,k)=0.0D0
                        ENDIF
                        d=thetanew(:,k)-thetatmp
                    ! CALL DBLEPR("d", -1, maxval(abs(d)), 1)
                        theta_sum=thetanew(:,k)+thetatmp
                        IF(any(d/=0.0D0)) THEN
                            dif=max(dif,maxval(abs(d)))
                            loss_diff = sum(d*(0.5*theta_sum-r(:,k)-sigma(k,k)*thetatmp))
                            penalty_diff = al*pf(k)*(sqrt(dot_product(thetanew(:,k),thetanew(:,k))) &
                            - sqrt(dot_product(thetatmp,thetatmp)))
                            dev = dev + loss_diff + penalty_diff
                            ab=spread(sigma(:,k),dim=1,ncopies=nk)*spread(d,dim=2,ncopies=nvars)
                            r=r-ab
                        ENDIF
                    ENDDO
                    IF(abs((dev-dev_tmp)/dev)<sml) EXIT
                    IF(dif<eps) EXIT
                 ! CALL DBLEPR("dif_inner", -1, dif, 1)

                    IF(npass > maxit) THEN
                        jerr=-l
                        RETURN
                    ENDIF
                ENDDO
            ENDDO
            IF(ni>pmax) EXIT
!--- this is the final check ------------------------
            vrg=1
            DO j=1,ni
                IF(maxval(abs(thetanew(:,m(j))-thetaold(:,m(j))))>=eps) THEN
                    vrg=0
                    EXIT
                ENDIF
            ENDDO
            IF(vrg==1) EXIT
            ! test deviance loop
            dev1_new = 0.0
            dev2_new = 0.0
            dev1_old = 0.0
            dev2_old = 0.0
            DO jj = 1,nk
                tmp1_new = dot_product(MATMUL(thetanew(jj,m(1:ni)),sigma(m(1:ni),m(1:ni))),thetanew(jj,m(1:ni)))
                tmp2_new = dot_product(thetanew(jj,m(1:ni)),delta(jj,m(1:ni)))
                dev1_new = dev1_new + tmp1_new
                dev2_new = dev2_new + tmp2_new
                tmp1_old = dot_product(MATMUL(thetaold(jj,m(1:ni)),sigma(m(1:ni),m(1:ni))),thetaold(jj,m(1:ni)))
                tmp2_old = dot_product(thetaold(jj,m(1:ni)),delta(jj,m(1:ni)))
                dev1_old = dev1_old + tmp1_old
                dev2_old = dev2_old + tmp2_old
            ENDDO
            dev3_new = al * sum(pf(m(1:ni)) * sqrt(sum(thetanew(:,m(1:ni)) * thetanew(:,m(1:ni)), DIM = 1)))
            dev3_old = al * sum(pf(m(1:ni)) * sqrt(sum(thetaold(:,m(1:ni)) * thetaold(:,m(1:ni)), DIM = 1)))
            dev_new = 0.5 * dev1_new - dev2_new + dev3_new
            dev_old = 0.5 * dev1_old - dev2_old + dev3_old
            IF(verbose==1) THEN
                CALL intpr('Current Lambda',-1,l,1)
                CALL dblepr('Obj-func Jump',-1,abs((dev_new-dev_old)/dev_new),1)
            ENDIF
            IF(abs((dev_new-dev_old)/dev_new)<sml) EXIT
            ! test deviance loop end
        ENDDO
!--- final update variable save results------------
        IF(ni>pmax) THEN
            jerr=-10000-l
            EXIT
        ENDIF
        IF(ni>0) theta(:,1:ni,l)=thetanew(:,m(1:ni))
        me = count(maxval(abs(theta(:,1:ni,l)),dim=1)/=0.0D0)
        IF(me>dfmax) THEN
      jerr=-20000-l
      EXIT
    ENDIF
        obj(l) = dev_new
        ntheta(l)=ni
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        IF (flmin >= 1.0D0) CYCLE
    ENDDO
    RETURN
END SUBROUTINE msda


! designed for msda with nonconvex penalties
!  changed pf (nvars) to pfmat (nvars, nlam)
! --------------------------------------------------
SUBROUTINE msda_ncx(obj,nk,nvars,sigma,delta,pfmat,dfmax,pmax,nlam,flmin,ulam,&
  eps,maxit,sml,verbose,nalam,theta,m,ntheta,alam,npass,jerr)
! --------------------------------------------------
    IMPLICIT NONE
    ! - - - arg types - - -
    DOUBLE PRECISION, PARAMETER :: big=9.9E30
    DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
    INTEGER, PARAMETER :: mnlam = 6
    INTEGER::mnl
    INTEGER::nk
    INTEGER::nvars
    INTEGER::dfmax
    INTEGER::pmax
    INTEGER::nlam
    INTEGER::nalam
    INTEGER::npass
    INTEGER::jerr
    INTEGER::verbose
    INTEGER::maxit
    INTEGER::m(pmax)
    INTEGER::ntheta(nlam)
    DOUBLE PRECISION::flmin
    DOUBLE PRECISION::eps
    DOUBLE PRECISION::sml
    DOUBLE PRECISION::sigma(nvars,nvars)
    DOUBLE PRECISION::delta(nk,nvars)
    DOUBLE PRECISION::pfmat(nvars, nlam)
    DOUBLE PRECISION::ulam(nlam)
    DOUBLE PRECISION::theta(nk,pmax,nlam)
    DOUBLE PRECISION::alam(nlam)
    DOUBLE PRECISION::obj(nlam)
    ! - - - local declarations - - -
    INTEGER::mm(nvars)
    INTEGER::k
    INTEGER::j
    INTEGER::jj
    INTEGER::l
    INTEGER::vrg
    INTEGER::ni
    INTEGER::me
    DOUBLE PRECISION::dif
    DOUBLE PRECISION::v
    DOUBLE PRECISION::al
    DOUBLE PRECISION::alf
    DOUBLE PRECISION::unorm
    DOUBLE PRECISION::thetanew(nk,nvars)
    DOUBLE PRECISION::thetaold(nk,nvars)
    DOUBLE PRECISION::r(nk,nvars)
    DOUBLE PRECISION::ab(nk,nvars)
    DOUBLE PRECISION::d(nk)
    DOUBLE PRECISION::theta_sum(nk)
    DOUBLE PRECISION::thetatmp(nk)
    DOUBLE PRECISION::u(nk)
    DOUBLE PRECISION::loss_diff
    DOUBLE PRECISION::penalty_diff
    DOUBLE PRECISION::dev
    DOUBLE PRECISION::dev_tmp
    DOUBLE PRECISION::tmp1_new
    DOUBLE PRECISION::tmp2_new
    DOUBLE PRECISION::dev_new
    DOUBLE PRECISION::dev1_new
    DOUBLE PRECISION::dev2_new
    DOUBLE PRECISION::dev3_new
    DOUBLE PRECISION::tmp1_old
    DOUBLE PRECISION::tmp2_old
    DOUBLE PRECISION::dev_old
    DOUBLE PRECISION::dev1_old
    DOUBLE PRECISION::dev2_old
    DOUBLE PRECISION::dev3_old
! - - - begin - - -
    IF(maxval(pfmat) <= 0.0D0) THEN
        jerr=10000
        RETURN
    ENDIF
    pfmat=max(0.0D0, pfmat)
! - - - some initial setup - - -
    mnl = Min (mnlam, nlam)
    r = delta
    thetanew=0.0D0
    thetaold=0.0D0
    dev=0.0D0
    m=0
    mm=0
    npass=0
    ni=npass
! --------- lambda loop ----------------------------
    IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf=flmin**(1.0D0/(nlam-1.0D0))
    ENDIF
    DO l=1,nlam
        IF(flmin>=1.0D0) THEN
            al=ulam(l)
        ELSE
            IF(l > 2) THEN
                al=al*alf
            ELSE IF(l==1) THEN
                al=big
            ELSE IF(l==2) THEN
                al=0.0D0
                DO j = 1,nvars
                    IF(pfmat(j, l)>0.0D0) THEN
                            u = delta(:,j)
                            v = sqrt(dot_product(u,u))
                            al=max(al,v/pfmat(j, l))
                    ENDIF
                END DO
                al=al*alf
            ENDIF
        ENDIF
! --------- outer loop ----------------------------
        DO
            IF(ni>0) thetaold(:,m(1:ni))=thetanew(:,m(1:ni))
! --middle loop-------------------------------------
            DO
                npass=npass+1
                dif=0.0D0
                dev_tmp = dev
                DO k=1,nvars
                    thetatmp=thetanew(:,k)
                    u = r(:,k)/sigma(k,k) + thetatmp
                    unorm = sqrt(dot_product(u,u))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! check this?                    
                    v = unorm-al*pfmat(k, l)/sigma(k,k)
                    IF(v > 0.0D0) THEN
                        thetanew(:,k) = v*u/unorm
                    ELSE
                        thetanew(:,k)=0.0D0
                    ENDIF
                    d=thetanew(:,k)-thetatmp
                    ! CALL DBLEPR("d", -1, maxval(abs(d)), 1)
                    theta_sum=thetanew(:,k)+thetatmp
                    IF(any(d/=0.0D0)) THEN
                        dif=max(dif,maxval(abs(d)))
                        loss_diff = sum(d*(0.5*theta_sum-r(:,k)-sigma(k,k)*thetatmp))
                        penalty_diff = al*pfmat(k, l)*(sqrt(dot_product(thetanew(:,k),thetanew(:,k))) &
                        - sqrt(dot_product(thetatmp,thetatmp)))
                        dev = dev + loss_diff + penalty_diff
                        ab=spread(sigma(:,k),dim=1,ncopies=nk)*spread(d,dim=2,ncopies=nvars)
                        r=r-ab
                        IF(mm(k)==0) THEN
                            ni=ni+1
                            IF(ni>pmax) EXIT
                                mm(k)=ni
                                m(ni)=k
                            ENDIF
                        ENDIF
                ENDDO
                IF(abs((dev-dev_tmp)/dev)<sml) EXIT
                IF(ni>pmax) EXIT
                IF(dif<eps) EXIT
                 ! CALL DBLEPR("dif_middle", -1, dif, 1)
                IF(npass > maxit) THEN
                    jerr=-l
                    RETURN
               ENDIF
! --inner loop----------------------
                DO
                    npass=npass+1
                    dif=0.0D0
                    dev_tmp = dev
                    DO j=1,ni
                        k=m(j)
                        thetatmp=thetanew(:,k)
                        u = r(:,k)/sigma(k,k) + thetatmp
                        unorm = sqrt(dot_product(u,u))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! check this?
                        v = unorm-al*pfmat(k, l)/sigma(k,k)
                        IF(v > 0.0D0) THEN
                            thetanew(:,k) = v*u/unorm
                        ELSE
                            thetanew(:,k)=0.0D0
                        ENDIF
                        d=thetanew(:,k)-thetatmp
                    ! CALL DBLEPR("d", -1, maxval(abs(d)), 1)
                        theta_sum=thetanew(:,k)+thetatmp
                        IF(any(d/=0.0D0)) THEN
                            dif=max(dif,maxval(abs(d)))
                            loss_diff = sum(d*(0.5*theta_sum-r(:,k)-sigma(k,k)*thetatmp))
                            penalty_diff = al*pfmat(k, l)*(sqrt(dot_product(thetanew(:,k),thetanew(:,k))) &
                            - sqrt(dot_product(thetatmp,thetatmp)))
                            dev = dev + loss_diff + penalty_diff
                            ab=spread(sigma(:,k),dim=1,ncopies=nk)*spread(d,dim=2,ncopies=nvars)
                            r=r-ab
                        ENDIF
                    ENDDO
                    IF(abs((dev-dev_tmp)/dev)<sml) EXIT
                    IF(dif<eps) EXIT
                 ! CALL DBLEPR("dif_inner", -1, dif, 1)

                    IF(npass > maxit) THEN
                        jerr=-l
                        RETURN
                    ENDIF
                ENDDO
            ENDDO
            IF(ni>pmax) EXIT
!--- this is the final check ------------------------
            vrg=1
            DO j=1,ni
                IF(maxval(abs(thetanew(:,m(j))-thetaold(:,m(j))))>=eps) THEN
                    vrg=0
                    EXIT
                ENDIF
            ENDDO
            IF(vrg==1) EXIT
            ! test deviance loop
            dev1_new = 0.0
            dev2_new = 0.0
            dev1_old = 0.0
            dev2_old = 0.0
            DO jj = 1,nk
                tmp1_new = dot_product(MATMUL(thetanew(jj,m(1:ni)),sigma(m(1:ni),m(1:ni))),thetanew(jj,m(1:ni)))
                tmp2_new = dot_product(thetanew(jj,m(1:ni)),delta(jj,m(1:ni)))
                dev1_new = dev1_new + tmp1_new
                dev2_new = dev2_new + tmp2_new
                tmp1_old = dot_product(MATMUL(thetaold(jj,m(1:ni)),sigma(m(1:ni),m(1:ni))),thetaold(jj,m(1:ni)))
                tmp2_old = dot_product(thetaold(jj,m(1:ni)),delta(jj,m(1:ni)))
                dev1_old = dev1_old + tmp1_old
                dev2_old = dev2_old + tmp2_old
            ENDDO
            dev3_new = al * sum(pfmat(m(1:ni), l) * sqrt(sum(thetanew(:,m(1:ni)) * thetanew(:,m(1:ni)), DIM = 1)))
            dev3_old = al * sum(pfmat(m(1:ni), l) * sqrt(sum(thetaold(:,m(1:ni)) * thetaold(:,m(1:ni)), DIM = 1)))
            dev_new = 0.5 * dev1_new - dev2_new + dev3_new
            dev_old = 0.5 * dev1_old - dev2_old + dev3_old
            IF(verbose==1) THEN
                CALL intpr('Current Lambda',-1,l,1)
                CALL dblepr('Obj-func Jump',-1,abs((dev_new-dev_old)/dev_new),1)
            ENDIF
            IF(abs((dev_new-dev_old)/dev_new)<sml) EXIT
            ! test deviance loop end
        ENDDO
!--- final update variable save results------------
        IF(ni>pmax) THEN
            jerr=-10000-l
            EXIT
        ENDIF
        IF(ni>0) theta(:,1:ni,l)=thetanew(:,m(1:ni))
        me = count(maxval(abs(theta(:,1:ni,l)),dim=1)/=0.0D0)
        IF(me>dfmax) THEN
      jerr=-20000-l
      EXIT
    ENDIF
        obj(l) = dev_new
        ntheta(l)=ni
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        IF (flmin >= 1.0D0) CYCLE
    ENDDO
    RETURN
END SUBROUTINE msda_ncx

