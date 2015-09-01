smax=17
nr1=0l
nread=4644l
mesh=dblarr(2400)
cs=dblarr(50)
func=dblarr(2,2400)
ksir=dblarr(2400,nread)
etar=dblarr(2400,nread)
lmode=intarr(nread)
nmode=intarr(nread)
numode=intarr(nread)
openr,1,'/scratch/seismo/schou/swapjcd/WB-2.1_models_l5bi.d.15_amde.1',/f77_unformatted
readu,1,nr1,mesh
for nfunc=0,nread-1 do begin
  readu,1,cs,func
  l=nint(cs(17))
  n=nint(cs(18))
  nu=1000*cs(26)
; Already multiplied by sqrt(rho). Multiply by r and correct for JCD scaling.
  ksir(*,nfunc)=func(0,*)/sqrt(mesh)
  etar(*,nfunc)=func(1,*)/sqrt(mesh)
; Fix center
  ksir(0,nfunc)=0
  etar(0,nfunc)=0
  lmode(nfunc)=l
  nmode(nfunc)=n
  numode(nfunc)=nu
end
close,1

lmax=max(lmode)

; Set integration weights
setint,2,mesh,intw

; Set up v_2s+1 from eq. 7
vx=dblarr(smax+1,lmax+1)
for l=0l,lmax do begin
  vx(0,l)=1
  for j=1l,smax do begin
    vx(j,l)=-vx(j-1,l)*(2*j+2D0)*(2*j+1D0)*(l+j+1D0)*(l-j)/(j*(j+1D0)*(2*l+2*j+2D0)*(2*l+2*j+1D0))
  end
end

; v_2s+1 for l->infinity
vinf=dblarr(smax+1)
vinf(0)=1
for j=1l,smax do vinf(j)=-vinf(j-1)*(2*j+2D0)*(2*j+1D0)/(4*j*(j+1D0))

; Calculate factor to correct for different definition of polynomials in
; latitude in Pijpers and calcsplit code
sratio=dblarr(smax+1)
sratio(0)=1
for s=1l,smax do sratio(s)=-sratio(s-1)*(2*s+1D0)*(2*s-1D0)/((4*s+1D0)*(4*s-1D0))

betax=dblarr(smax+1,nfunc)
for nfunc=0,nread-1 do begin
  l=lmode(nfunc)
  n=nmode(nfunc)
  ll=sqrt(l*(l+1D))
; Calculate normalization integral, f1 and f2 from eqs. 5 and 6.
  inl=total((ksir(*,nfunc)^2+etar(*,nfunc)^2)*intw)
  f1=(ksir(*,nfunc)^2-2*ksir(*,nfunc)*etar(*,nfunc)/ll+etar(*,nfunc)^2)/inl
  f2=etar(*,nfunc)^2/ll^2/inl
  for s=0,smax do begin
    fls=(f1-f2*(2*s+2)*(2*s+1)/2)*vx(s,l)
    beta=total(fls*intw)
    cg=total(mesh*fls*intw)/beta
    beta_old=beta/vinf(s) ; Beta asymptotically equal to 1
    beta_js=beta_old*sratio(s) ; Beta from calcsplit code
if ((l eq 300) and (n eq 0)) then print,l,n,s,2*s+1,beta,beta_old,beta_js,cg
    betax(s,nfunc)=beta_js
  end
end

end

