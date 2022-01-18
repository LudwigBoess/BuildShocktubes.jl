function get_bfield_from_angle(theta::Real;
                            reduction_scale::Real=1.0e-10)

    shock_normal = [1.0, 0.0, 0.0]
    # convert to radians
    theta *= π/180.0

    bfld = zeros(3)

    bfld[1] = sin(90.0 * π/180.0 - theta)
    bfld[2] = sin(theta)

    bfld ./= sqrt(bfld[1]^2 + bfld[2]^2 + bfld[3]^2)

    #bfld += shock_normal

    bfld .*= reduction_scale

    return bfld

end


function read_Bfield(filename, Npart)
    f = open(filename, "r")
    bfld = read!(f, Matrix{Float32}(undef, 3, Npart))
    close(f)

    return bfld
end



function findii(x::Real, xtab::Array{<:Real})

    n = length(xtab)

    imin  = -1
    imax  = -1
    dmin1 = maximum(xtab) - minimum(xtab)
    dmin2 = maximum(xtab) - minimum(xtab)

    for i = 2:n
        d = x - xtab[i]
        if (0.0 < d < dmin1)
            imin = i
            dmin1 = d
        end

        if ( d < 0.0  && -d < dmin2 )
            imax = i
            dmin2 = -d
        end
    end
    return imin, imax
end

function interpolate_components(Bh::Array{<:Real}, 
                                dx::Real, dy::Real, dz::Real,
                                ixmin::Integer, ixmax::Integer,
                                iymin::Integer, iymax::Integer,
                                izmin::Integer, izmax::Integer )

    fx1y1z1 = Bh[ixmin, iymin, izmin]
    fx1y1z2 = Bh[ixmin, iymin, izmax]
    fx1y2z1 = Bh[ixmin, iymax, izmin]
    fx1y2z2 = Bh[ixmin, iymax, izmax]
    fx2y1z1 = Bh[ixmax, iymin, izmin]
    fx2y1z2 = Bh[ixmax, iymin, izmax]
    fx2y2z1 = Bh[ixmax, iymax, izmin]
    fx2y2z2 = Bh[ixmax, iymax, izmax]

    b = ( 1.0 - dx ) * ( 1.0 - dy )* ( 1.0 - dz ) * fx1y1z1 +
        ( 1.0 - dx ) * ( 1.0 - dy )*         dz   * fx1y1z2 +
        ( 1.0 - dx ) *         dy  * ( 1.0 - dz ) * fx1y2z1 +
        ( 1.0 - dx ) *         dy  *         dz   * fx1y2z2 +
                dx   * ( 1.0 - dy) * ( 1.0 - dz)  * fx2y1z1 +
                dx   * ( 1.0 - dy) *         dz   * fx2y1z2 +
                dx   *         dy  * ( 1.0 - dz)  * fx2y2z1 +
                dx   *         dy  *         dz   * fx2y2z2

    return b
end

function setup_turb_B(pos::Array{<:Real}, 
                    npart::Integer, B0::Real)

    r = zeros(npart)

    for i = 1:length(r)
        r[i] = norm(pos[i,:])
    end

    k = 1.0./r

    # set grid
    # NFFT  = 117
    # NFFT2 = 58

    NFFT  = 233
    NFFT2 = 116

    minn = maximum(abs.(k)) / minimum(abs.(k))
    mink = minimum(abs.(k))

    k0 = maximum(abs.(k))

    α = 5.0/3.0

    A0 = k0^(-α) * B0^2

    α2 = 0.5*α

    nnn = (NFFT2+1)^3 * 4

    val = randn(nnn)
    ϕ2  = rand(nnn) .* 2π

    kx = zeros(Int64, NFFT)
    ky = zeros(Int64, NFFT)
    kz = zeros(Int64, NFFT)

    xp = zeros(Float64, NFFT)
    yp = zeros(Float64, NFFT)
    zp = zeros(Float64, NFFT)

    Bhx = zeros(NFFT, NFFT, NFFT)
    Bhy = zeros(NFFT, NFFT, NFFT)
    Bhz = zeros(NFFT, NFFT, NFFT)


    ii = 1

    @showprogress "Constructing grid..." for iz=0:NFFT2

        niz = (iz > 0) ? (NFFT-iz) : 0
        kz[iz+1]  =  iz
        kz[niz+1] = -iz

        for iy = 0:NFFT2

            niy = ( iy > 0 ) ? (NFFT - iy) : 0
            ky[iy+1]  =  iy
            ky[niy+1] = -iy

            for ix = 0:NFFT2

                nix = ( ix > 0 ) ? (NFFT - ix) : 0
                kx[ix+1]  =  ix
                kx[nix+1] = -ix

                k2 = ( kx[ix+1]^2 + ky[iy+1]^2 + kz[iz+1]^2 ) * mink^2
                PSP = ( k2 != 0.0 ) ? sqrt(A0 * k2^α2 ) :  0.0

                # first quadrant

                Bh = PSP * val[ii]
                if  kz[iz+1] == 0
                    ϕ1 = π/4.0
                else
                    ϕ1 = atan( - (kx[ix+1]*cos(ϕ2[ii])
                            +  ky[iy+1]*sin(ϕ2[ii]))
                            /  kz[iz+1] )
                end

                Bhx[ix+1, iy+1, iz+1]    = Bh*cos(ϕ1)*cos(ϕ2[ii])
                Bhy[ix+1, iy+1, iz+1]    = Bh*cos(ϕ1)*sin(ϕ2[ii])
                Bhz[ix+1, iy+1, iz+1]    = Bh*sin(ϕ1)

                Bhx[nix+1, niy+1, niz+1] = Bh*cos(ϕ1)*cos(ϕ2[ii])
                Bhy[nix+1, niy+1, niz+1] = Bh*cos(ϕ1)*sin(ϕ2[ii])
                Bhz[nix+1, niy+1, niz+1] = Bh*sin(ϕ1)


                # second quadrant
                ii += 1
                Bh = PSP * val[ii]
                if  kz[iz+1] == 0
                    ϕ1 = π/4.0
                else
                    ϕ1 = atan( - (kx[ix+1]*cos(ϕ2[ii])
                            +  ky[iy+1]*sin(ϕ2[ii]))
                            /  kz[iz+1] )
                end

                Bhx[nix+1, iy+1, iz+1]   = Bh*cos(ϕ1)*cos(ϕ2[ii])
                Bhy[nix+1, iy+1, iz+1]   = Bh*cos(ϕ1)*sin(ϕ2[ii])
                Bhz[nix+1, iy+1, iz+1]   = Bh*sin(ϕ1)

                Bhx[ix+1, niy+1, niz+1]  = Bh*cos(ϕ1)*cos(ϕ2[ii])
                Bhy[ix+1, niy+1, niz+1]  = Bh*cos(ϕ1)*sin(ϕ2[ii])
                Bhz[ix+1, niy+1, niz+1]  = Bh*sin(ϕ1)

                # third quadrant
                ii += 1
                Bh = PSP * val[ii]
                if  kz[iz+1] == 0
                    ϕ1 = π/4.0
                else
                    ϕ1 = atan( - (kx[ix+1]*cos(ϕ2[ii])
                            +  ky[iy+1]*sin(ϕ2[ii]))
                            /  kz[iz+1] )
                end

                Bhx[ix+1, niy+1, iz+1]   = Bh*cos(ϕ1)*cos(ϕ2[ii])
                Bhy[ix+1, niy+1, iz+1]   = Bh*cos(ϕ1)*sin(ϕ2[ii])
                Bhz[ix+1, niy+1, iz+1]   = Bh*sin(ϕ1)

                Bhx[nix+1, iy+1, niz+1]  = Bh*cos(ϕ1)*cos(ϕ2[ii])
                Bhy[nix+1, iy+1, niz+1]  = Bh*cos(ϕ1)*sin(ϕ2[ii])
                Bhz[nix+1, iy+1, niz+1]  = Bh*sin(ϕ1)

                # fourth quadrant
                ii += 1
                Bh = PSP * val[ii]
                if  kz[iz+1] == 0
                    ϕ1 = π/4.0
                else
                    ϕ1 = atan( - (kx[ix+1]*cos(ϕ2[ii])
                            +  ky[iy+1]*sin(ϕ2[ii]))
                            /  kz[iz+1] )
                end

                Bhx[nix+1, niy+1, iz+1]  = Bh*cos(ϕ1)*cos(ϕ2[ii])
                Bhy[nix+1, niy+1, iz+1]  = Bh*cos(ϕ1)*sin(ϕ2[ii])
                Bhz[nix+1, niy+1, iz+1]  = Bh*sin(ϕ1)

                Bhx[ix+1, iy+1, niz+1]   = Bh*cos(ϕ1)*cos(ϕ2[ii])
                Bhy[ix+1, iy+1, niz+1]   = Bh*cos(ϕ1)*sin(ϕ2[ii])
                Bhz[ix+1, iy+1, niz+1]   = Bh*sin(ϕ1)


                ii += 1
            end

        end

    end

    #Bhx
    # Bhx_dummy = Bhx
    # Bhx = Bhx_dummy

    for iz = 1:NFFT
        xp[iz] = 1.0/(kx[iz] * mink)
        yp[iz] = 1.0/(ky[iz] * mink)
        zp[iz] = 1.0/(kz[iz] * mink)
    end

    @info "FFTW"

    Bhx = fft(Bhx)
    Bfehler1 = sum(abs.(imag.(Bhx)))
    Bhx = real.(Bhx)
    Bfehler2 = sum(abs.(Bhx))

    Bhy = fft(Bhy)
    #Bfehler1 = sum(abs.(imag.(Bhx)))
    Bhy = real.(Bhy)
    #Bfehler2 = sum(abs.(Bhx))

    Bhz = fft(Bhz)
    #Bfehler1 = sum(abs.(imag.(Bhx)))
    Bhz = real.(Bhz)
    #Bfehler2 = sum(abs.(Bhx))

    anzz=NFFT^3
    bmean  = sqrt((sum(Bhx)/anzz)^2+(sum(Bhy)/anzz)^2+(sum(Bhz)/anzz)^2)
    bmean2 = sum(sqrt.(Bhx.^2+Bhy.^2+Bhz.^2))/anzz

    println("Gitter: <|B|>= $bmean2, |<B>|= $bmean")

    bx = zeros(npart)
    by = zeros(npart)
    bz = zeros(npart)

    @showprogress "Interpolating B..." for i = 1:npart

        xvec = pos[i,:]

        ixmin, ixmax = findii(pos[i,1], xp)
        iymin, iymax = findii(pos[i,2], yp)
        izmin, izmax = findii(pos[i,3], zp)

        dx = ( pos[i,1] - xp[ixmin] ) / ( xp[ixmax] - xp[ixmin] )
        dy = ( pos[i,2] - yp[ixmin] ) / ( yp[ixmax] - yp[ixmin] )
        dz = ( pos[i,3] - zp[ixmin] ) / ( zp[ixmax] - zp[ixmin] )


        bx[i] = interpolate_components(Bhx, dx, dy, dz, 
                                    ixmin, ixmax,
                                    iymin, iymax,
                                    izmin, izmax)

        # y component
        by[i] = interpolate_components(Bhy, dx, dy, dz, 
                                    ixmin, ixmax,
                                    iymin, iymax,
                                    izmin, izmax)

        # z component
        bz[i] = interpolate_components(Bhz, dx, dy, dz, 
                                    ixmin, ixmax,
                                    iymin, iymax,
                                    izmin, izmax)
    end

    return [bx by bz]
end