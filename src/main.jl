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
        if par.B_filename == ""
            B = Array{Float32,2}(undef, 3, N)

            B[1, left_part ] .= Float32(par.B[1,1])
            B[1, right_part] .= Float32(par.B[1,2])

            B[2, left_part ] .= Float32(par.B[2,1])
            B[2, right_part] .= Float32(par.B[2,2])

            B[3, left_part ] .= Float32(par.B[3,1])
            B[3, right_part] .= Float32(par.B[3,2])
        else
            println("reading Bfld data")
            n_small = size(x_small,2)
            n_large = size(x_large,2)

            B_small = read_Bfield(par.B_filename * ".dat", n_small)
            B_large = read_Bfield(par.B_filename * "_large.dat", n_large)

            B_left = build_B_tube(B_large, n_blocks)
            B_right = build_B_tube(B_small, n_blocks)

            B = [B_left B_right]

            N_B = size(B,ndims(B))

            B = Float32.(B)
        end
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

    if par.CRp || par.CRe

        println("setting up CR spectral")

        if par.density_step
            rho_left  = 1.0
            rho_right = 0.125
        else
            rho_left  = 1.0
            rho_right = 1.0
        end
        
        # if electrons required inject them at a given
        # electron to proton ration Kep
        if par.CRe
            Kep = par.Kep
        else
            Kep = 0.0
        end

        # total pressure given in parameter file
        P_tot_left  = U_to_P(par.U[1], rho_left  )
        P_tot_right = U_to_P(par.U[2], rho_right )

        # compute thermal pressure component from Xcr = Pcr / Pth 
        #  Ptot = Pth + Pcr = Pth + Xcr * Pth
        #       = Pth * ( 1.0 + Xcr )
        P_th_left  = P_tot_left  / ( 1.0 + par.Xcr[1] )
        P_th_right = P_tot_right / ( 1.0 + par.Xcr[2] )

        U[left_part]  .= Float32( P_to_U( P_th_left,  rho_left  ) )
        U[right_part] .= Float32( P_to_U( P_th_right, rho_right ) )

        P_cr_left      = P_tot_left  - P_th_left
        P_cr_right     = P_tot_right - P_th_right

        # boundaries of CR spectrum
        bounds = set_bounds(par.p_lim, par.Nbins)

        # Protons
        pNorm          = Matrix{Float32}(undef, par.Nbins, N)
        pNorm         .= -Inf 
        pSlope         = par.cr_slope .* ones(Float32, par.Nbins, N)
        pCut           = par.p_lim[2] .* ones(Float32, N)

        if par.CRp

            # left part 
            init_norm_left = find_init_norm(( 1.0 - Kep ) * P_cr_left, par.cr_slope, par.p_lim[1], par.p_lim[2])
            norm_left      = init_powerlaw(init_norm_left, par.cr_slope, bounds)
            for i ∈ left_part
                pNorm[:,i] = norm_left
            end

            # right part 
            init_norm_right = find_init_norm(( 1.0 - Kep ) * P_cr_right, par.cr_slope, par.p_lim[1], par.p_lim[2])
            norm_right      = init_powerlaw(init_norm_right, par.cr_slope, bounds)
            for i ∈ right_part
                pNorm[:,i] = norm_right
            end

        end

        println("protons done!")

        # Electrons
        eNorm          = Matrix{Float32}(undef, par.Nbins, N)
        eNorm         .= -Inf 
        eSlope         = par.cr_slope .* ones(Float32, par.Nbins, N)
        eCut           = par.p_lim[2] .* ones(Float32, N)

        if par.CRe

            # left part 
            init_norm_left = find_init_norm( Kep * P_cr_left, par.cr_slope, par.p_lim[1], par.p_lim[2])
            norm_left      = init_powerlaw(init_norm_left, par.cr_slope, bounds)
            for i ∈ left_part
                eNorm[:,i] = norm_left
            end

            # right part 
            init_norm_right = find_init_norm( Kep * P_cr_right, par.cr_slope, par.p_lim[1], par.p_lim[2])
            norm_right      = init_powerlaw(init_norm_right, par.cr_slope, bounds)
            for i ∈ right_part
                eNorm[:,i] = norm_right
            end

        end

        println("electrons done!")

    else
        U[left_part]  .= Float32(par.U[1])
        U[right_part] .= Float32(par.U[2])
    end

    println("done")

    println("writing ic file")

    if !arepo

        # write Gadget file
        f = open(par.output_file, "w")
        write_header( f, head)
        write_block(  f, x,        "POS")
        write_block(  f, vel,      "VEL")
        write_block(  f, id,       "ID")
        write_block(  f, U,        "U")
        write_block(  f, hsml_out, "HSML")
        write_block(  f, B,        "BFLD")
        if par.CRp || par.CRe
            write_block( f, pNorm,  "CRpN")
            write_block( f, eNorm,  "CReN")

            write_block( f, pSlope, "CRpS")
            write_block( f, eSlope, "CReS")

            write_block( f, pCut,   "CRpC")
            write_block( f, eCut,   "CReC")
        end
        close(f)

    else
        # write arepo file 
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