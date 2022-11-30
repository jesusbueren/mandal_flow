	SUBROUTINE powell(p,xi,ftol,iter,fret)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	USE nr, ONLY : linmin
	IMPLICIT NONE
	REAL(DP), DIMENSION(7), INTENT(INOUT) :: p
	REAL(DP), DIMENSION(7,7), INTENT(INOUT) :: xi
	INTEGER(I4B), INTENT(OUT) :: iter
	REAL(DP), INTENT(IN) :: ftol
	REAL(DP), INTENT(OUT) :: fret
	INTERFACE
		FUNCTION log_likelihood2(p)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), DIMENSION(7), INTENT(IN) :: p
		REAL(DP) :: log_likelihood2
		END FUNCTION log_likelihood2
	END INTERFACE
	INTEGER(I4B), PARAMETER :: ITMAX=200
	REAL(DP), PARAMETER :: TINY=1.0e-25_dp
	INTEGER(I4B) :: i,ibig,n
	REAL(DP) :: del,fp,fptt,t
	REAL(DP), DIMENSION(size(p)) :: pt,ptt,xit
	n=assert_eq(size(p),size(xi,1),size(xi,2),'powell')
	fret=log_likelihood2(p)
	pt(:)=p(:)
	iter=0
	do
		iter=iter+1
		fp=fret
		ibig=0
		del=0.0
		do i=1,n
			xit(:)=xi(:,i)
			fptt=fret
			call linmin(p,xit,fret)
			if (fptt-fret > del) then
				del=fptt-fret
				ibig=i
			end if
		end do
		if (2.0_dp*(fp-fret) <= ftol*(abs(fp)+abs(fret))+TINY) RETURN
		if (iter == ITMAX) call &
			nrerror('powell exceeding maximum iterations')
		ptt(:)=2.0_dp*p(:)-pt(:)
		xit(:)=p(:)-pt(:)
		pt(:)=p(:)
		fptt=log_likelihood2(ptt)
		if (fptt >= fp) cycle
		t=2.0_dp*(fp-2.0_dp*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
		if (t >= 0.0) cycle
		call linmin(p,xit,fret)
		xi(:,ibig)=xi(:,n)
		xi(:,n)=xit(:)
	end do
	END SUBROUTINE powell
