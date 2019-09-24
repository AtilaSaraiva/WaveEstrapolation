module difFinitas
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

        subroutine waveEstrap (ordem,Nz,Nx,Nt,Nb,dx,dz,dt,wavelet,campoVel,fontes,snaps)
            !$acc routine(waveEstrap)
            !=============================================!
            ! Subrotina para fazer a propagação de uma
            ! onda num campo de velocidade.
            !=============================================!

            ! Entrada e Saída
            integer,intent(in) :: Nx,Nz,Nt,fontes(:,:),ordem,Nb
            real,intent(in)    :: campoVel(Nz,Nx),dx,dt,dz,wavelet(:)
            real,intent(out):: snaps(Nz,Nx,Nt/20-1)

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
                P_futuro = (2.0*P_atual - P_passado  + (dt**2.)*(campoVelExt**2.)*lap)! * mascaraAtenuacaoBorda

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
            integer:: i,j,lz,lx

            lz=nzb
            lx=nxb
            !$acc kernels
            do i =1,nb
                do j=1,nxb
                    p2(i,j)=p2(i,j)*coef(i)
                    p2(lz,j)=p2(lz,j)*coef(i)
                end do
                do j=1,nzb
                    p2(j,i)=p2(j,i)*coef(i)
                    p2(j,lx)=p2(j,lx)*coef(i)
                enddo
                lx=lx-1
                lz=lz-1
            enddo
            !$acc end kernels

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
            !$acc kernels
            do i=-Nb+1,0
                campoVelExt(i,1:Nx) = campoVel(1,:)
                campoVelExt(1:Nz,i) = campoVel(:,1)
                j = -i+1
                campoVelExt(Nz+j,1:Nx) = campoVel(Nx,:)
                campoVelExt(1:Nz,Nx+j) = campoVel(:,Nz)
            end do

            ! Fazendo o mesmo para os valores nos 4 cantos do campo de
            ! velocidade original.
            campoVelExt(-Nb+1:0,-Nb+1:0) = campoVel(1,1)
            campoVelExt(-Nb+1:0,Nx+1:Nx+Nb) = campoVel(1,Nx)
            campoVelExt(Nz+1:Nz+Nb,-Nb+1:0) = campoVel(Nz,1)
            campoVelExt(Nz+1:Nz+Nb,Nx+1:Nx+Nb) = campoVel(Nz,Nx)
            !$acc end kernels

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

            !$acc parallel loop
            do j=in_n,lim_nx
                !$acc loop
                do i=in_n,lim_nz
                    Pxx = 0.0
                    Pzz = 0.0

                    !$acc loop reduction(+:Pxx) reduction(+:Pzz)
                    do k=1,ordem+1
                        ! Derivada em x
                        Pxx = Pxx + coef(k)*P(i,j+k-in_n)
                        ! Derivada em z
                        Pzz = Pzz + coef(k)*P(i+k-in_n,j)
                    end do

                    laplaciano(i,j) = Pxx / dx**2. + Pzz / (dz**2.)
                end do
            end do

            return
        end function laplaciano

end module difFinitas

program wave
    use difFinitas
    use rsf
    implicit none
    logical:: verb
    ! Arquivos de entrada
    type(file) :: FcampoVel,Fpulso,Fsnaps
    type(axa) :: at,az,ax
    ! Parâmetros do modelo
    integer:: nx,nz,nt,nb,ordem
    real:: dt,dx,dz,start,finish
    ! Campo de velocidade
    real,allocatable:: campoVel(:,:)
    ! Wavelet
    real,allocatable:: pulso(:)
    real,allocatable:: snaps(:,:,:)
    ! Auxiliares
    integer:: i,j,h1A,h2A,fontes(2,1),k,counter


    ! Inicializando o madagascar
    call sf_init()
    ! Variavel de verbose
    call from_par("verb",verb,.false.)

    ! Definindo variáveis que vão conter os arquivos
    FcampoVel = rsf_input("vel")
    Fpulso = rsf_input("wav")
    Fsnaps = rsf_output("out")

    ! Retirando do header as informações de geometria
    call from_par(FcampoVel,"n1",az%n)
    call from_par(FcampoVel,"n2",ax%n)
    call from_par(FcampoVel,"d1",az%d)
    call from_par(FcampoVel,"d2",ax%d)
 !  call iaxa(FcampoVel,az,1)
 !  call iaxa(FcampoVel,ax,2)
    call iaxa(Fpulso,at,1)

 !  ! Definindo a geometria do output
    call to_par (Fsnaps,"d1",az%d)
    call to_par (Fsnaps,"d2",ax%d)
    call to_par (Fsnaps,"d3",at%d*20)
    call to_par (Fsnaps,"n1",az%n)
    call to_par (Fsnaps,"n2",ax%n)
    call to_par (Fsnaps,"n3",at%n/20-1)
    call to_par (Fsnaps,"o1",0)
    call to_par (Fsnaps,"o2",0)
    call to_par (Fsnaps,"o3",0)
!   call oaxa(Fsnaps,az,1)
!   call oaxa(Fsnaps,ax,2)
!   call oaxa(Fsnaps,at,3)

 !  ! Alocando variáveis e lendo
    allocate(campoVel(az%n,ax%n))
    campoVel=0.
    call rsf_read(FcampoVel,campoVel)

    allocate(pulso(at%n))
    pulso=0.
    call rsf_read(Fpulso,pulso)

    allocate(snaps(az%n,ax%n,at%n/20-1))

    ! Retirando da variável de geometria as informações de geometria
    dt = at%d
    dz = az%d
    dx = ax%d

    nt = at%n
    nz = az%n
    nx = ax%n

    ! Tamanho de borda
    nb = 0.2 * nx

    ! Ordem do operador laplaciano
    ordem = 8

    ! Posição da fonte
    fontes(:,1) = [1,nx/2]

    ! Chamada da subrotina de propagação da onda
    call cpu_time(start)
    call waveEstrap (ordem,nz,nx,nt,nb,dx,dz,dt,pulso,campoVel,fontes,snaps)
    call cpu_time(finish)

    write(0,*) "Tempo total: ",finish-start

    ! Escrevendo o arquivo de output
    call rsf_write(Fsnaps,snaps)

    ! Saino
    call exit(0)
end program wave

