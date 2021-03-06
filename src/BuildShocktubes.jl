module BuildShocktubes

    using Printf
    using LinearAlgebra: norm
    using FFTW
    using ProgressMeter
    using GadgetIO

    export ShockParameters,
            P_to_U,
            get_bfield_from_angle,
            setup_shocktube

    struct ShockParameters

        glass_file::String
        output_file::String
        n_blocks::Int64
        v::Array{Float64,2}
        B::Array{Float64,2}
        U::Vector{Float64}
        B0::Float64
        density_step::Bool
        turb::Bool

        function ShockParameters(glass_file::String="", output_file::String="",
                                U::Vector{Float64}=zeros(2),
                                B::Array{Float64,2}=zeros(3,2),
                                v::Array{Float64,2}=zeros(3,2);
                                turb::Bool=false,
                                B0::Float64=0.0,
                                density_step::Bool=true,
                                n_blocks::Int64=70)

            new(glass_file, output_file, n_blocks,
                v, B, U, B0, density_step, turb)

        end
    end

    function P_to_U(P::Real, rho::Real)
        return P/((2.0/3.0) * rho)
    end

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

    function getLargeBox(x_in::Array{<:Real}, hsml_in::Array{<:Real}=[0.0])

        n_part = size(x_in, ndims(x_in))
        
        x = Matrix{eltype(x_in[1])}(undef, 3, 8n_part)

        count = 0

        for ix = 0:1, iy = 0:1, iz = 0:1
            for i = 1:n_part
                x[1,count*n_part+i] = x_in[1,i] + ix
                x[2,count*n_part+i] = x_in[2,i] + iy
                x[3,count*n_part+i] = x_in[3,i] + iz
            end
            count += 1
        end

        x = 0.5 .* x

        if hsml_in != [0.0]
            h = Vector{eltype(hsml_in[1])}(undef, 8n_part)
            for i = 0:7
                h[i*n_part+1:(i+1)*n_part] = 0.5 .* hsml_in
            end
            return x, h
        else
            return x
        end
    end

    function buildTube(x0::Array{<:Real}, n_blocks::Integer, hsml0=0, offset=0)

        n_part = size(x0, ndims(x0))
        x = Matrix{eltype(x0[1])}(undef, 3, n_blocks*n_part)

        for i = 0:n_blocks-1
            x[:,i*n_part+1:(i+1)*n_part] = i .* [1.0, 0.0, 0.0] .+ x0 .+ [offset, 0.0, 0.0]
        end

        if hsml0 != 0
            hsml = Vector{eltype(hsml0[1])}(undef, n_blocks*n_part)
            for i = 0:n_blocks-1
                hsml[i*n_part+1:(i+1)*n_part] = hsml0
            end

            return Float32.(x), Float32.(hsml)
        end

        return Float32.(x)
    end

    function build_B_tube(B_in::Array{<:Real}, n_blocks::Integer)

        n_part = size(B_in,2)
        B = Array{eltype(B_in[1]),2}(undef, 3, n_part)

        for i = 0:n_blocks-1
            B[:,i*n_part+1:(i+1)*n_part] .= B_in
        end

        return Float32.(B)
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

    function setup_shocktube(par::ShockParameters; arepo::Bool=false)

        println("reading glass file")
        pos_info = InfoLine("POS", Float32, 3, [1,0,0,0,0,0])
        hsml_info = InfoLine("HSML", Float32, 1, [1,0,0,0,0,0])

        h = head_to_obj(par.glass_file)

        pos  = Float32(1.0/h.boxsize) .* read_block(par.glass_file, "POS", info=pos_info, parttype=0)
        hsml = Float32(1.0/h.boxsize) .* read_block(par.glass_file, "HSML", info=hsml_info, parttype=0)

        println("read!")

        println("setting up tube")
        x_large, hsml_l = getLargeBox(pos, hsml)

        x_small = pos

        n_blocks = par.n_blocks

        # build tubes
        println("Building x tubes")
        if par.density_step
            m = 1.0/size(x_large,2)
            x_left, hsml_left = buildTube(x_large, n_blocks, hsml_l)
        else
            m = 1.0/size(x_small,2)
            x_left, hsml_left = buildTube(x_small, n_blocks, hsml)
        end
        x_right, hsml_right = buildTube(x_small, n_blocks, hsml, n_blocks)

        x        = [x_left x_right]
        hsml_out = [hsml_left; hsml_right]
        N = size(x, ndims(x))

        left_part = findall(x[1,:] .<= Float64(par.n_blocks))
        right_part = findall(x[1,:] .> Float64(par.n_blocks))

        # set up random magnetic field
        if par.turb

            n_large = size(x_large,2)
            B_large = setup_turb_B(x_large, n_large, par.B0)

            n_small = size(x_small,2)
            B_small = setup_turb_B(x_small, n_small, par.B0)

            @info "Building B tube"
            
            B_left = build_B_tube(B_large, n_blocks)
            B_right = build_B_tube(B_small, n_blocks)
            B = [B_left B_right]
            N_B = size(B,ndims(B))

            # println("Setting up turb B")
            # B = setup_turb_B(x, N, par.B0)
            # N_B = length(B[:,1])

            # if N != N_B
            #     error("Error in tube building x!\nN = $N\nN_B = $N_B")
            # end

            B = Float32.(B)


        else
            B = Array{Float32,2}(undef, 3, N)

            B[1, left_part ] .= Float32(par.B[1,1])
            B[1, right_part] .= Float32(par.B[1,2])

            B[2, left_part ] .= Float32(par.B[2,1])
            B[2, right_part] .= Float32(par.B[2,2])

            B[3, left_part ] .= Float32(par.B[3,1])
            B[3, right_part] .= Float32(par.B[3,2])
        end

        println("done")

        println("Checking uniqueness of positions.")
        unique_check = (size(unique(x, dims=2),2) == N) ? true : false

        if unique_check == false
            return "ERROR! Overlapping particles!"
        end


        println(minimum(x[1,:]), " ", maximum(x[1,:]))

        println("Number of particles = ", N)

        println("Assigning shock parameters")

        head = head_to_obj(par.glass_file)
        head.boxsize = 100000.0
        head.npart[1] = N
        head.nall[1] = N
        head.massarr[1] = m

        vel = Array{Float32,2}(undef, 3, N)

        vel[1, left_part ] .= Float32(par.v[1,1])
        vel[1, right_part] .= Float32(par.v[1,2])

        vel[2, left_part ] .= Float32(par.v[2,1])
        vel[2, right_part] .= Float32(par.v[2,2])

        vel[3, left_part ] .= Float32(par.v[3,1])
        vel[3, right_part] .= Float32(par.v[3,2])

        id = UInt32.(collect(0:N-1))

        U = Vector{Float32}(undef, N)

        U[left_part] .= Float32(par.U[1])
        U[right_part] .= Float32(par.U[2])

        println("done")

        println("writing ic file")

        if !arepo

            f = open(par.output_file, "w")
            write_header( f, head)
            write_block(  f, x,        "POS")
            write_block(  f, vel,      "VEL")
            write_block(  f, id,       "ID")
            write_block(  f, U,        "U")
            write_block(  f, hsml_out, "HSML")
            write_block(  f, B,        "BFLD")
            close(f)

        else

            ρ = Vector{Float32}(undef, N)
            ρ[left_part] .= Float32(1.0)
            if par.density_step 
                ρ[right_part] .= Float32(0.125)
            else
                ρ[right_part] .= Float32(1.0)
            end

            f = open(par.output_file, "w")
            write_header( f, head)
            write_block(  f, x,   "POS")
            write_block(  f, vel, "VEL")
            write_block(  f, id,  "ID")
            write_block(  f, U,   "U")
            write_block(  f, ρ,   "MASS")
            write_block(  f, B,   "BFLD")
            close(f)

        end

    end

end # module
