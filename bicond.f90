program bicond
    implicit none
    integer:: n=1,i, nx=10, ny=10, j=1 !nx= numero de celulas em x, ny=numero de celulas em y, j= contador
    double precision:: alpha=200.0d0, cfl=0.5d0!(constante condutividade)
    double precision:: dt, dx,dy, x, y, z, time !passo to, espaço, z=dt*alpha/dx**2
    double precision:: top, bot, wes, eas, tol=1.0d-6 !t no top, baixo, oeste, leste, toleranca consdição parada
    double precision:: max_val1=0, max_val2=1 !valores maximos do vetor
    double precision, allocatable, dimension(:)::t, ti! t n+1, t n e inicial 
    double precision, allocatable, dimension(:,:):: a
    x=10.0d0
    y=x
    dx=x/(real(nx-1))
    dy=dx
    dt=cfl*((dx**2.0d0)/alpha)
    top=600.0d0
    bot=top
    wes=top
    eas=top
    allocate(a(nx*ny,nx*ny))
    write(*,*) 'a'
    allocate(t(nx*ny))
    write(*,*) 't'
    allocate(ti(nx*ny))
    write(*,*)'ti'
    ti=0.0d0
    z=(dt*alpha)/(dx**2.0d0)
    a=0
    i=1
    !#####################
    ! escrita da matriz!
    do while (i<=nx*ny)
        a(i,i)=z-4.0d0
        if (i<nx*ny) then
            a(i,i+1)=z
        end if
        if (i>=2) then
            a(i,i-1)=z
        end if
        if (i<nx*(ny-1)) then
            a(i,i+nx)=z
        end if
        if (i>nx) then
            a(i,i-nx)=z
        end if
        i=i+1
    end do
    !########################
    !gauss-seidel
    i=1
    do while (time<2.0d0)
        max_val1=0
        max_val2=1
        j=1
        do while(abs(max_val2-max_val1)>tol)
            max_val2=max_val1
            do i=1,nx*ny
                if(i==1)then !celula topo esquerdo
                    t(i)=(ti(i)-top-wes-z*(t(i+1)+t(i+nx)))/(z-4.0d0)
                
                else if(i==nx) then !celula top direito
                    t(i)=(ti(i)-top-eas-z*(t(i-1)+t(i+nx)))/(z-4.0d0)
                
                else if (i==(nx*(ny-1)+1)) then !celula bot esquerdo
                    t(i)=(ti(i)-bot-wes-z*(t(i+1)+t(i-nx)))/(z-4.0d0)
                
                else if (i==(nx*ny)) then !celula bot direito
                    t(i)=(ti(i)-bot-eas-z*(t(i-1)+t(i-nx)))/(z-4.0d0)
                
                else if (i>=2 .and. i<=nx-1) then !celulas top ñ cantos
                    t(i)=(ti(i)-top-z*(t(i-1)+t(i+1)+t(i+nx)))/(z-4.0d0)
                
                else if (i==(n*nx)+1 .and. n<=ny-2) then !celulas west
                    t(i)=(ti(i)-wes-z*(t(i+1)+t(i+nx)+t(i-nx)))/(z-4.0d0)
                    n=n+1
                else if (i==n*nx .and. n<=ny-1) then !celulas east
                    t(i)=(ti(i)-eas-z*(t(i-1)+t(i+nx)+t(i-nx)))/(z-4.0d0)
                else if (i>=nx*(ny-1)+1 .and. i<=nx*ny-1) then !celulas bot
                    t(i)=(ti(i)-bot-z*(t(i-1)+t(i+1)+t(i-nx)))/(z-4.0d0)
                else
                    t(i)=(ti(i)-z*(t(i-1)+t(i+1)+t(i-nx)+t(i+nx)))/(z-4.0d0)
                end if
            end do
            max_val1=(maxval(abs(t)))
            write(*,*) j,'          ', abs(max_val2-max_val1), time
            j=j+1
            
        end do
        ti=t
        time=time+dt 
    end do 
    write(*,*) t


end program
