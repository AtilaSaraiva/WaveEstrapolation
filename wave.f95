program wave
    use difFinitas
    use rsf
    implicit none
    logical:: verb
    real:: sx,sz,gxbeg,gzbeg,jgx
    integer:: isx,isz,igxbeg,igzbeg,o1,o2,nr
    ! Arquivos de entrada
    type(file) :: FcampoVel,Fpulso,Fsnaps,Fdata
    type(axa) :: at,az,ax
    ! Parâmetros do modelo
    integer:: nx,nz,nt,nb,ordem
    real:: dt,dx,dz
    double precision::start,finish
    ! Campo de velocidade
    real,allocatable:: campoVel(:,:)
    ! Wavelet
    real,allocatable:: pulso(:)
    real,allocatable:: snaps(:,:,:),dat(:,:)
    ! Auxiliares
    integer:: i,j,h1A,h2A,fontes(2,1),k,counter


    ! Inicializando o madagascar
    call sf_init()
    ! Variavel de verbose
    call from_par("verb",verb,.false.)

    ! Definindo variáveis que vão conter os arquivos
    FcampoVel = rsf_input("vel")
    Fpulso = rsf_input("wav")
    Fdata = rsf_output("out")
    Fsnaps = rsf_output("snaps")

    ! Retirando do header as informações de geometria
    call from_par(FcampoVel,"n1",az%n)
    call from_par(FcampoVel,"n2",ax%n)
    call from_par(FcampoVel,"d1",az%d)
    call from_par(FcampoVel,"d2",ax%d)
    call from_par(FcampoVel,"o1",o1)
    call from_par(FcampoVel,"o2",o2)
 !  call iaxa(FcampoVel,az,1)
 !  call iaxa(FcampoVel,ax,2)
    call iaxa(Fpulso,at,1)
    call from_par('sx',sx)
    call from_par('sz',sz)
    call from_par('gxbeg',gxbeg)
    call from_par('gzbeg',gzbeg)
    call from_par('jgx',jgx)
    call from_par('nr',nr)

    !igxbeg = int((gxbeg - o2 + 1) / ax%d)
    isx = int((sx - o2) / ax%d) + 1
    isz = int((sz - o1) / az%d) + 1
    igxbeg = int((gxbeg - o2)/ ax%d) + 1
    igzbeg = int((gzbeg - o1) / az%d) + 1

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
    call to_par (Fdata,"d1",at%d)
    call to_par (Fdata,"d2",jgx)
    call to_par (Fdata,"o2",gxbeg)
    call to_par (Fdata,"o1",gzbeg)
    call to_par (Fdata,"n1",at%n)
    call to_par (Fdata,"n2",nr)

 !  ! Alocando variáveis e lendo
    allocate(campoVel(az%n,ax%n))
    campoVel=0.
    call rsf_read(FcampoVel,campoVel)

    allocate(pulso(at%n))
    pulso=0.
    call rsf_read(Fpulso,pulso)

    allocate(snaps(az%n,ax%n,at%n/20-1))
    allocate(dat(at%n,nr))

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
    fontes(:,1) = [isz,isx]


    ! Chamada da subrotina de propagação da onda
    start = omp_get_wtime()
    call waveEstrap (ordem,nz,nx,nt,nb,igzbeg,igxbeg,nr,jgx,dx,dz,dt,pulso,campoVel,fontes,snaps,dat)
    finish = omp_get_wtime()

    write(0,*) "Tempo total: ",finish-start

    ! Escrevendo o arquivo de output
    call rsf_write(Fsnaps,snaps)
    call rsf_write(Fdata,dat)

    ! Saino
    call exit(0)
end program wave

