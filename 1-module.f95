module difFinitas
    use omp_lib
    !=============================================!
    ! Módulo desenvolvido para efetuar a modelagem
    ! direta da propagação de uma onda ao resolver
    ! a equação da onda acústica num campo de vel-
    ! ocidade específico.
    !
    ! Autor: Átila Saraiva Quintela Soares
    ! Email: atilasaraiva@gmail.com
    !
    ! Contribuição: Protásio
    !=============================================!
    implicit none
    integer,parameter:: dp = kind(0.d0)
    real,parameter:: numPi = 3.14159263
    contains

        subroutine waveEstrap (ordem,Nz,Nx,Nt,Nb,igzbeg,igxbeg,nr,jgx,dx,dz,dt,wavelet,campoVel,fontes,snaps,dat)
            !=============================================!
            ! Subrotina para fazer a propagação de uma
            ! onda num campo de velocidade.
            !=============================================!

            ! Entrada e Saída
            integer,intent(in) :: Nx,Nz,Nt,fontes(:,:),ordem,Nb,nr,igxbeg,igzbeg
            real,intent(in)    :: campoVel(Nz,Nx),dx,dt,dz,wavelet(:),jgx
            real,intent(out):: snaps(Nz,Nx,Nt/20-1),dat(nt,nr)

            ! Variáveis Auxiliares
            real,allocatable   :: P_futuro(:,:),P_passado(:,:),P_atual(:,:) ! Campos de pressão
            real,allocatable   :: lap(:,:)                                  ! Laplaciano
            real,allocatable   :: campoVelExt(:,:)                          ! Campo de velocidade com borda de atenuacao
            real               :: coefAtenuacao(Nb)                         ! Coefientes de atenuacao
            character(len=30)  :: filename                                  ! Nome do Arquivo
            real               :: prod=1.0                                  ! Aux de produto
            integer            :: i,j,k                                     ! Contadores
            integer            :: f1,f2,lFontes                             ! Relacionados a posição e num de fontes
            integer            :: Nxb,Nzb                                   ! Dimensões com a borda
            integer:: passo

            ! Criando arquivo para guardar os snaps
            filename = 'snap.ad'
            !open(30,file=filename,status='replace',recl=4*Nz*Nx,form='unformatted',access='direct')

            ! Extraindo a quantidade de fontes
            lFontes = ubound(fontes,2)


            ! Alocando arrays
            Nxb = Nx + 2*Nb
            Nzb = Nz + 2*Nb
            allocate(campoVelExt(Nzb,Nxb))
            allocate(P_futuro(Nzb,Nxb),P_passado(Nzb,Nxb),P_atual(Nzb,Nxb),lap(Nzb,Nxb))

            ! Extendendo o campo de velocidade com a borda de atenuação
            campoVelExt = extenCampoVel(Nx,Nz,Nb,campoVel)

            ! Zerando os campos
            P_passado = 0.0
            P_atual = 0.0
            P_futuro = 0.0

            ! Calculando os coeficientes de Atenuação para serem multiplicados na Borda
            coefAtenuacao = coeficientesDeAtenuacao(Nb)

            j=1
            do k=2,Nt-1
                if(mod(k,100) == 0) write(0,*) 'it',k

                ! Atenuando os campos passado e presente
                call atenuacao(nxb,nzb,nb,coefAtenuacao,P_passado)
                call atenuacao(nxb,nzb,nb,coefAtenuacao,P_atual)

                ! Calculando o laplaciano
                lap = laplaciano(Nb,ordem,Nxb,Nzb,dx,dz,P_atual)

                ! Resolvendo a equação da onda por diferenças finitas
                P_futuro = (2.0*P_atual - P_passado  + (dt**2.)*(campoVelExt**2.)*lap)

                ! Inserindo a wavelet nas posições de fonte
                do i=1,lFontes
                    f1=fontes(1,i)+Nb
                    f2=fontes(2,i)+Nb
                    P_futuro(f1,f2) = wavelet(k) + P_futuro(f1,f2)
                end do

                ! Se a iteração for múltiplo de 20, escrever o campo no arquivo
                if(mod(k,20) == 0)then
                    !write(30,rec=j) P_futuro(Nb+1:Nz+Nb,Nb+1:Nx+Nb)
                    snaps(:,:,j) = P_futuro(Nb+1:Nz+Nb,Nb+1:Nx+Nb)
                    j = j + 1
                end if

                passo = int(jgx / dx)
                do i=1,nr
                    dat(k,i) = P_futuro(Nb+igzbeg,Nb+igxbeg+i*passo-1)
                end do
                P_passado = P_atual
                P_atual = P_futuro
            end do
            !close(30)
        end subroutine waveEstrap

        function coeficientesDeAtenuacao(Nb)
            !=============================================!
            ! Função para calcular os coeficientes de
            ! atenuação da onda na borda.
            !=============================================!

            ! Entrada e Saída
            integer,intent(in):: Nb
            real:: coeficientesDeAtenuacao(Nb)
            ! Auxiliares
            integer:: i

            ! Calculando os coeficientes
            do i =1,nb
               coeficientesDeAtenuacao(i) = exp(-(0.008*(nb-i))**2)
            enddo

            return
        end function coeficientesDeAtenuacao

        subroutine atenuacao(nxb,nzb,nb,coef,p2)
            !=============================================!
            ! Subrotina para a atenuação de borda da onda
            !=============================================!

            ! Entrada e Saída
            integer:: nxb, nzb,nb
            real,intent(in):: coef(Nb)
            real, dimension (:,:):: p2(nzb,nxb)
            ! Auxiliares
            integer:: i,j,nz,nx,i2,j2

            nz=nzb-2*nb
            nx=nxb-2*nb
            !$omp parallel
            !$omp do
            do i=1,nxb
                do j=1,nb
                    p2(i,j) = p2(i,j) * coef(j)
                    p2(i,nx+nb+j-1) = p2(i,nx+nb+j-1) * coef(nb-j+1)
                    i2=j
                    j2=i
                    p2(i2,j2) = p2(i2,j2) * coef(i2)
                    p2(nz+nb+i2-1,j2) = p2(nz+nb+i2-1,j2) * coef(nb-i2+1)
                end do
            end do
            !$omp end do
            !$omp end parallel

            return
        end subroutine atenuacao

        function extenCampoVel (Nx,Nz,Nb,campoVel) result(campoVelExt)
            !=============================================!
            ! Função para extensão do campo de velocidade
            ! para que seja possível se criar uma borda de
            ! atenuação.
            !=============================================!
            ! Entrada e Saída
            integer,intent(in):: Nx,Nz,Nb
            real,intent(in):: campoVel(Nz,Nx)
            real:: campoVelExt(-Nb+1:Nz+Nb,-Nb+1:Nx+Nb)
            ! Auxiliares
            integer:: i,j

            ! Copiando o campo original no centro do extendido
            campoVelExt(1:Nz,1:Nx) = campoVel

            ! Copiando os valores do contorno do campo de vel original
            ! na borda do campo extendido
            !$omp parallel
            !$omp do
            do i=-Nb+1,0
                campoVelExt(i,1:Nx) = campoVel(1,:)
                campoVelExt(1:Nz,i) = campoVel(:,1)
                j = -i+1
                campoVelExt(Nz+j,1:Nx) = campoVel(Nx,:)
                campoVelExt(1:Nz,Nx+j) = campoVel(:,Nz)
            end do
            !$omp end do nowait
            !$omp end parallel

            ! Fazendo o mesmo para os valores nos 4 cantos do campo de
            ! velocidade original.
            campoVelExt(-Nb+1:0,-Nb+1:0) = campoVel(1,1)
            campoVelExt(-Nb+1:0,Nx+1:Nx+Nb) = campoVel(1,Nx)
            campoVelExt(Nz+1:Nz+Nb,-Nb+1:0) = campoVel(Nz,1)
            campoVelExt(Nz+1:Nz+Nb,Nx+1:Nx+Nb) = campoVel(Nz,Nx)

        end function extenCampoVel

        function laplaciano(Nb,ordem,Nx,Nz,dx,dz,P)
            !=============================================!
            ! Subrotina para cálculo do laplaciano de um
            ! campo de pressão, com ordens de 2,4,6 e 8.
            !=============================================!

            ! Entrada e Saída
            integer,intent(in):: Nx,Nz,ordem,Nb
            real,intent(in):: dx,dz,P(Nz,Nx)
            real:: laplaciano(Nx,Nz)
            ! Matrizes auxiliares
            real:: coef(ordem+1)
            real:: Pxx,Pzz
            integer:: i,j,k,in_n,lim_nx,lim_nz
            laplaciano = 0.0

            ! Calculo dos coeficientes do operador de
            ! diferenças finitas
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
                    return
            end select

            ! Definindo os limites em que o operador vai trabalhar
            in_n   = ordem / 2 + 1
            lim_nx = Nx - ordem / 2
            lim_nz = Nz - ordem / 2

            ! Calculando o laplaciano

            !$omp parallel shared(coef,P,laplaciano)
            !$omp do reduction(+:Pxx) reduction(+:Pzz) schedule(dynamic,1) private(i,j)
            do j=in_n,lim_nx
                do i=in_n,lim_nz
                    Pxx = 0.0
                    Pzz = 0.0

                    do k=1,ordem+1
                        ! Derivada em x
                        Pxx = Pxx + coef(k)*P(i,j+k-in_n)
                        ! Derivada em z
                        Pzz = Pzz + coef(k)*P(i+k-in_n,j)
                    end do

                    laplaciano(i,j) = Pxx / dx**2. + Pzz / (dz**2.)
                end do
            end do
            !$omp end do nowait
            !$omp barrier
            !$omp end parallel

            return
        end function laplaciano

end module difFinitas
