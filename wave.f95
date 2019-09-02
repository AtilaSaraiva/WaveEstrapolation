module difFinitas
    implicit none
    integer,parameter:: dp = kind(0.d0)
    real,parameter:: numPi = 3.14159263
    contains

        subroutine waveEstrap (ordem,Nz,Nx,Nt,dx,dz,dt,wavelet,campoVel,fontes)
            integer,intent(in):: Nx,Nz,Nt,fontes(:,:),ordem
            real,intent(in):: campoVel(Nz,Nx),dx,dt,dz,wavelet(:)
            !real:: P_futuro(Nz,Nx),P_passado(Nz,Nx),P_atual(Nz,Nx)
            real,allocatable:: P_futuro(:,:),P_passado(:,:),P_atual(:,:)
            real,allocatable:: lap(:,:)
            real,allocatable:: campoVelExt(:,:),mascaraAtenuacaoBorda(:,:)
            character(len=30) filename

            real:: waveletExtendida(Nt),prod=1.0
            integer:: i,j,k,L,f1,f2,j1,lFontes,counter,Nb,Nxb,Nzb

            filename = 'snap.ad'

            open(30,file=filename,status='replace',recl=4*Nz*Nx,form='unformatted',access='direct')

            ! Extraindo a quantidade de fontes
            lFontes = ubound(fontes,2)


            Nb = 30
            Nxb = Nx + 2*Nb
            Nzb = Nz + 2*Nb
            allocate(campoVelExt(Nzb,Nxb),mascaraAtenuacaoBorda(Nzb,Nxb))
            allocate(P_futuro(Nzb,Nxb),P_passado(Nzb,Nxb),P_atual(Nzb,Nxb),lap(Nzb,Nxb))

            campoVelExt = extenCampoVel(Nx,Nz,Nb,campoVel)
            mascaraAtenuacaoBorda = mascaraAtenuacao(Nx,Nz,Nb,dx,dz)

            P_passado = 0.0
            P_atual = 0.0
            P_futuro = 0.0


            j=1
            !do k=2,Nt-1
            do k=2,1000
                write(0,*) 'it',k
                call atenuacao(nxb,nzb,nb,P_passado)
                call atenuacao(nxb,nzb,nb,P_atual)

                lap = laplaciano(Nb,ordem,Nxb,Nzb,dx,dz,P_atual)
                P_futuro = (2.0*P_atual - P_passado  + (dt**2.)*(campoVelExt**2.)*lap)! * mascaraAtenuacaoBorda

                do i=1,lFontes
                    f1=fontes(1,i)+Nb
                    f2=fontes(2,i)+Nb
                    P_futuro(f1,f2) = wavelet(k) + P_futuro(f1,f2)
                end do

                if(mod(k,20) == 0)then
                    write(30,rec=j) P_futuro(Nb+1:Nz+Nb,Nb+1:Nx+Nb)
                    j = j + 1
                end if

                P_passado = P_atual
                P_atual = P_futuro
            end do
            close(30)
        end subroutine waveEstrap

        subroutine atenuacao(nxb,nzb,nb,p2)
            integer:: nxb, nzb,nb,i,j,lz,lx
            real, dimension (:,:):: p2(nzb,nxb)
            real:: alpha,alpha1

            lz=nzb
            do i =1,nb
               do j=1,nxb
                   alpha = exp(-(0.008*(nb-i))**2)
                   !alpha = 0.95**((float(Nb-i)/float(Nb))**2.0)
                   p2(i,j)=p2(i,j)*alpha
                   p2(lz,j)=p2(lz,j)*alpha
               end do
               lz=lz-1
            enddo

            lx=nxb
            do i =1,nb
               do j=1,nzb
                   !alpha = 0.95**(((i-1)/Nb)**2.0)
                   alpha = exp(-(0.008*(nb-i))**2)
                   p2(j,i)=p2(j,i)*alpha
                   p2(j,lx)=p2(j,lx)*alpha
               enddo
               lx=lx-1
            enddo

            return
        end subroutine atenuacao

        function mascaraAtenuacao(Nx,Nz,Nb,dx,dz)
            integer,intent(in):: Nx,Nz,Nb
            real,intent(in)::dx,dz
            real:: mascaraAtenuacao(-Nb+1:Nz+Nb,-Nb+1:Nx+Nb),alpha
            integer:: i,j
            real:: x,z

            mascaraAtenuacao(1:Nz,1:Nx) = 1.0

            do j=-Nb+1,0
                x = j*dx
                z = j*dz
                alpha=0.15**(-x/(Nb*dx))
                mascaraAtenuacao(:,j) = alpha
                mascaraAtenuacao(:,-j+1+Nx) = alpha
                alpha=0.15**(-z/(Nb*dz))
                mascaraAtenuacao(j,1:Nx) = alpha
                mascaraAtenuacao(-j+1+Nx,1:Nx) = alpha
            end do

            open(40,file="mascara.ad",status='replace',recl=4*(Nz+2*Nb),form='unformatted',access='direct')

            do j=-Nb+1,Nx+Nb
                write(40,rec=j+Nb) mascaraAtenuacao(:,j)
            end do
            close(40)
        end function mascaraAtenuacao

        function extenCampoVel (Nx,Nz,Nb,campoVel) result(campoVelExt)
            integer,intent(in):: Nx,Nz,Nb
            real,intent(in):: campoVel(Nz,Nx)
            real:: campoVelExt(-Nb+1:Nz+Nb,-Nb+1:Nx+Nb)
            integer:: i,j

            campoVelExt(1:Nz,1:Nx) = campoVel
            do i=-Nb+1,0
                campoVelExt(i,1:Nx) = campoVel(1,:)
                campoVelExt(1:Nz,i) = campoVel(:,1)
                j = -i+1
                campoVelExt(Nz+j,1:Nx) = campoVel(Nx,:)
                campoVelExt(1:Nz,Nx+j) = campoVel(:,Nz)
            end do

            campoVelExt(-Nb+1:0,-Nb+1:0) = campoVel(1,1)
            campoVelExt(-Nb+1:0,Nx+1:Nx+Nb) = campoVel(1,Nx)
            campoVelExt(Nz+1:Nz+Nb,-Nb+1:0) = campoVel(Nz,1)
            campoVelExt(Nz+1:Nz+Nb,Nx+1:Nx+Nb) = campoVel(Nz,Nx)

        end function extenCampoVel

        function laplaciano(Nb,ordem,Nx,Nz,dx,dz,P)
            integer,intent(in):: Nx,Nz,ordem,Nb
            real,intent(in):: dx,dz,P(Nz,Nx)
            real:: laplaciano(Nx,Nz)
            ! Matrizes auxiliares
            real:: coef(ordem+1)
            real:: Pxx,Pzz
            integer:: i,j,k,in_n,lim_nx,lim_nz
            laplaciano = 0.0

            coef = 0.0
            select case (ordem)

                case(2)
                    coef(1) = 1.
                    coef(2) = -2.
                    coef(3) = 1.
                case(4)
                    coef(1) = (-1.)/12.
                    coef(2) = (4.)/3.
                    coef(3) = (-5.)/2.
                    coef(4) = (4.)/3.
                    coef(5) = (-1.)/12.
                case(6)
                    coef(1) = 1./90.
                    coef(2) = (-3.)/20.
                    coef(3) = 3./2.
                    coef(4) = (-49.)/18.
                    coef(5) = 3./2.
                    coef(6) = (-3.)/20.
                    coef(7) = 1./90.
                case(8)
                    coef(1) = -1./560.
                    coef(2) = 8./315.
                    coef(3) = (-1.)/5.
                    coef(4) = 8./5.
                    coef(5) = (-205.)/72.
                    coef(6) = 8./5.
                    coef(7) = (-1.)/5.
                    coef(8) = 8./315.
                    coef(9) = -1./560.
                case default
                    write(0,*) "Ordem inválida"

            end select

            in_n   = ordem / 2 + 1
            lim_nx = Nx - ordem / 2
            lim_nz = Nz - ordem / 2
            do j=in_n,lim_nx
                do i=in_n,lim_nz
                    Pxx = 0.0
                    Pzz = 0.0

                    do k=1,ordem+1
                        Pxx = Pxx + coef(k)*P(i,j+k-in_n)
                        Pzz = Pzz + coef(k)*P(i+k-in_n,j)
                    end do

    !                if(j<=Nb .or. j>=Nx-Nb+1 .or. i<=Nb .or. i>=Nz-Nb+1)then
    !                    laplaciano(i,j) = Pxx / (dx*10.0)**2. + Pzz / ((dz*10.0)**2.)
    !                else
                        laplaciano(i,j) = Pxx / dx**2. + Pzz / (dz**2.)
    !                end if
                end do
            end do

        end function laplaciano

end module difFinitas

program wave
    use difFinitas
    implicit none
    integer,parameter:: Nx=200,Nz=200,Nt=1500,ordem=8
    real,parameter:: h1=1000,h2=1000
    real,parameter:: dx=10.0,dz=10.0
    real:: dt, campoVel(Nz,Nx),t,fmax,tAux
    real:: wavelet(Nt)
    integer:: i,j,h1A,h2A,fontes(2,1),k,counter

    h1A = Nz/3.0
    h2A = h1A
    campoVel(1:h1A,1:Nx) = 2000.0
    campoVel(h1A:h2A+h1A,1:Nx) = 3500.0
    campoVel(h2A+h1A:Nz,1:Nx) = 4000.0

    fontes(:,1) = [Nz/2,Nx/2]

    dt=0.001

    open(50,file="pulso.ad",status='old',recl=4*Nt,form='unformatted',access='direct')
    !open(40,file="campovel.ad",status='replace',recl=4*Nz,form='unformatted',access='direct')

    read(50,rec=1) wavelet

    call waveEstrap (ordem,Nz,Nx,Nt,dx,dz,dt,wavelet,campoVel,fontes)
    close(50)
    !close(40)

end program wave
