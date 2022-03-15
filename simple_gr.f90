PROGRAM Simple
IMPLICIT NONE
INTEGER :: nx,ny,i,j,it=0,it2=0
REAL(KIND=8) :: lx,ly,dx,dx2,dy,dy2,ri,pr,re,rel,i_pr,i_re,i_dx,i_dy,i_dx2,i_dy2,dt,eps,eps_p,gr,i_gr
REAL(KIND=8) :: vt,vb,a,b,c,AA,ur,ul,BB,temp,al,u,v,err_x,err_y,ut,ub,vr,vl,err_t
REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: vx,vy,vx1,vy1,p,p1,p2,div,T,T1,vxs,vys
nx=64
ny=64
lx=1
ly=1
dx=lx/nx
dy=ly/ny
dx2=dx*dx
dy2=dy*dy
i_dx=1/dx
i_dy=1/dy
i_dx2=1/dx2
i_dy2=1/dy2
dt=0.00005 !timestep
pr=1 !numero prandtl
i_pr=1/pr
gr=10;
i_gr=1/gr;
rel=1.97 !relaxation factor
eps_p=1.0E-11 !precisione solutore poisson
eps=1.0E-10 !criterio per convergenza
ALLOCATE(vx(1:ny+2,1:nx+1))
ALLOCATE(vx1(1:ny+2,1:nx+1))
ALLOCATE(vy(1:ny+1,1:nx+2))
ALLOCATE(vy1(1:ny+1,1:nx+2))
ALLOCATE(p(1:ny+2,1:nx+2))
ALLOCATE(p1(1:ny+2,1:nx+2))
ALLOCATE(p2(1:ny+2,1:nx+2))
ALLOCATE(div(1:ny+1,1:nx+1))
ALLOCATE(T(1:ny+2,1:nx+2))
ALLOCATE(T1(1:ny+2,1:nx+2))
T1=0
T=0
vx=0
vy=0
p=0
p1=0
p2=0
ALLOCATE(vxs(1:ny,1:nx))
ALLOCATE(vys(1:ny,1:nx))

principale: DO
it=it+1
WRITE(*,*) "it n*", it, "it poisson" , it2 ,"differenza velocità", max(err_x,err_y), "differenza temp" , err_t
!componente x velocità
DO i=2,nx
 DO j=2,ny+1 
  
  vt=0.5*(vy(j-1,i)+vy(j-1,i+1))
  vb=0.5*(vy(j,i)+vy(j,i+1))
  ut=0.5*(vx(j-1,i)+vx(j,i))
  ub=0.5*(vx(j+1,i)+vx(j,i))
  ur=0.5*(vx(j,i)+vx(j,i+1))
  ul=0.5*(vx(j,i)+vx(j,i-1))
  
  a=-(ur**2-ul**2)*i_dx-(ub*vb-ut*vt)*i_dy
  !a=-(vx(j,i+1)**2-vx(j,i-1)**2)*0.5*i_dx-(vx(j+1,i)*vb-vx(j-1,i)*vt)*0.5*i_dy !convettivo
  b=(vx(j,i+1)-2*vx(j,i)+vx(j,i-1))*i_dx2 !diffusivo 1
  !b=(vy(j,i+1)-2*vy(j,i)+vy(j,i)-1)*i_dx2; !----ERRORE--------ERRORE--------ERRORE--------ERRORE--------ERRORE----
  c=(vx(j+1,i)-2*vx(j,i)+vx(j-1,i))*i_dy2 !diffusivo 2
  AA=-a+sqrt(i_gr)*(b+c)
  vx1(j,i)=vx(j,i)+dt*(AA-(p(j,i+1)-p(j,i))*i_dx) !velocita vx al passo n+1    
 END DO 
END DO 
!componente y velocità
DO i=2,nx+1
 DO j=2,ny 
  
  ur=0.5*(vx(j,i)+vx(j+1,i))
  ul=0.5*(vx(j,i-1)+vx(j+1,i-1))
  vr=0.5*(vy(j,i+1)+vy(j,i))
  vl=0.5*(vy(j,i-1)+vy(j,i))
  vt=0.5*(vy(j-1,i)+vy(j,i))
  vb=0.5*(vy(j+1,i)+vy(j,i))
  a=-(vb**2-vt**2)*i_dy-(vr*ur-vl*ul)*i_dx !convettivo
  !a=-(vy(j+1,i)**2-vy(j-1,i)**2)*0.5*i_dy-(vy(j,i+1)*ur-vy(j,i-1)*ul)*0.5*i_dx !convettivo
  b=(vy(j,i+1)-2*vy(j,i)+vy(j,i-1))*i_dx2 !diffusivo 1
  c=(vy(j+1,i)-2*vy(j,i)+vy(j-1,i))*i_dy2 !diffusivo 2
  temp=-(T1(j+1,i)+T1(j,i))*0.5 !termine galleggiamento
  BB=-a-temp+sqrt(i_gr)*(b+c) !approssimo gr>>1 gr/re^2
  vy1(j,i)=vy(j,i)+dt*(BB-(p(j+1,i)-p(j,i))*i_dy) !velocita vy al passo n+1


 END DO 
END DO 

!applico le condizioni al contorno
!bc vx
!parete sx e dx impermiabilità
vx1(2:ny+1,1)=0
vx1(2:ny+1,nx+1)=0
!parete top e bottom noslip
vx1(1,1:nx+1)=2*0-vx1(2,1:nx+1)
vx1(ny+2,1:nx+1)=-vx1(ny+1,1:nx+1)

!bc vy
!parete sx e dx noslip
vy1(1:ny+1,1)=-vy1(1:ny+1,2)
vy1(1:ny+1,nx+2)=-vy1(1:ny+1,nx+1)
!parete top e bottom impermiabilità
vy1(1,1:nx+2)=0
vy1(ny+1,1:nx+2)=0

!calcoliamo la divergenza ,come source dell'eq di poisson e poi risolviamo
!tutto
DO i=2,nx+1
    DO j=2,ny+1
    div(j,i)=+(vx1(j,i)-vx1(j,i-1))*i_dx+(vy1(j,i)-vy1(j-1,i))*i_dy
    END DO 
END DO

div=div*dx2*dy2/dt
a=dx2*0.5/(dy2+dx2)
b=dy2*0.5/(dy2+dx2)
c=0.5/(dx2+dy2)
!solver poisson
it2=0
poisson: DO
p2=p1
it2=it2+1
DO i=2,nx+1
 DO j=2,ny+1
      al=a*(p1(j+1,i)+p1(j-1,i))+b*(p1(j,i+1)+p1(j,i-1))-c*div(j,i)
      p1(j,i)=al*rel+(1-rel)*p1(j,i)
 END DO 
END DO
    !bc poisson,gradiente normale alla parete nullo
    p1(1:ny+2,1)=p1(1:ny+2,2) !parete sx
    p1(1:ny+2,nx+2)=p1(1:ny+2,nx+1) !parete dx
    p1(1,1:nx+2)=p1(2,1:nx+2) !parete top
    p1(ny+2,1:nx+2)=p1(ny+1,1:nx+2) !parete bottom

    !aggiungere gli spigoli
    p1(1,1)=p1(2,2)
    
    !p1(1,nx+2)=p1(1,nx+1) ---possibile ERRORE----
    p1(1,nx+2)=p1(2,nx+1) 
    p1(ny+2,nx+2)=p1(ny+1,nx+1)
    p1(ny+2,1)=p1(ny+1,2)
    ! p1(ny+2,1)=p1(ny+1,1) -------possibile ERRORE
    !
    IF (MAXVAL(ABS(p2-p1))<eps_p) EXIT
  END DO poisson
 !correggiamo le velocità 
 DO i=2,nx
  DO j=2,ny+1
   vx1(j,i)=vx1(j,i)-(p1(j,i+1)-p1(j,i))*i_dx*dt
   END DO 
   END DO
   
    DO i=2,nx+1
  DO j=2,ny
   vy1(j,i)=vy1(j,i)-(p1(j+1,i)-p1(j,i))*i_dy*dt
   END DO 
   END DO
   !campo scalare della temperatura 
   DO j=2,ny+1
    DO i=2,nx+1
    !ricavo le velocita
        u=0.5*(vx1(j,i)+vx1(j,i-1))
        v=0.5*(vy1(j,i)+vy1(j-1,i))
        a=0.5*u*i_dx*(T(j,i+1)-T(j,i-1))
        b=0.5*v*i_dy*(T(j+1,i)-T(j-1,i))
        
        c=sqrt(i_Gr)*i_pr*((T(j,i+1)-2*T(j,i)+T(j,i-1))*i_dx2+(T(j+1,i)-2*T(j,i)+T(j-1,i))*i_dy2)
        T1(j,i)=T(j,i)+dt*(-a-b+c)
    END DO
END DO
err_t=maxval(abs(T-T1))
!condizioni al contorno temperatura
T1(2:ny+1,1)=2*0-T1(2:ny+1,2) !parete sx fredda
T1(2:ny+1,nx+2)=2*1-T1(2:ny+1,nx+1) !parete dx sorgente calda
!T1(2:ny+1,nx+2)=2*1-T1(2:ny+1,nx+2) !!----ERRORE--------ERRORE--------ERRORE--------ERRORE----
T1(1,2:nx+1)=T1(2,2:nx+1) !parete top adiabatica
T1(ny+2,2:nx+1)=T1(ny+1,2:nx+1) !parete bottom adiabatica

T=T1;

err_x=maxval(abs(vx-vx1))
err_y=maxval(abs(vy-vy1)) 
vx=vx1
vy=vy1
p=p1+p
p1=p1*0
 
!IF (it>40000) EXIT !esco per numero di iterazioni
!END DO principale
WRITE(*,*) err_x,err_y
IF ((max(err_x,err_y)<eps) .AND. it>2) EXIT !esco per la differenza massima tra step raggiunta
END DO principale 



!torno dalla griglia staggered a una normale
DO j=1,ny
    DO i=1,nx
        vys(j,i)=0.5*(vy1(j,i+1)+vy1(j+1,i+1));
        vxs(j,i)=0.5*(vx1(j+1,i+1)+vx1(j+1,i));
    END DO
END DO
OPEN(10, file="vex.bin", form="unformatted")
WRITE(10) vxs
CLOSE(10)
OPEN(10, file="vey.bin", form="unformatted")
WRITE(10) vys
CLOSE(10)
OPEN(10, file="temp.bin", form="unformatted")
WRITE(10) T(2:ny+1,2:nx+1)
CLOSE(10)
WRITE(*,*) SIZE(vys)

END PROGRAM Simple
