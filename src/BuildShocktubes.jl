module BuildShocktubes

    using Printf
    using LinearAlgebra: norm
    using FFTW
    using ProgressMeter
    using GadgetIO

    export ShockParameters,
            P_to_U, U_to_P,
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
        CRp::Bool
        CRe::Bool
        Xcr::Vector{Float64}
        Kep::Float64
        cr_slope::Float64
        p_lim::Vector{Float64}
        Nbins::Int64
        B_filename::String

        function ShockParameters(glass_file::String="", output_file::String="",
                                U::Vector{Float64}=zeros(2),
                                B::Array{Float64,2}=zeros(3,2),
                                v::Array{Float64,2}=zeros(3,2);
                                turb::Bool=false,
                                B0::Float64=0.0,
                                CRp::Bool=false,
                                CRe::Bool=false,
                                Xcr::Vector{Float64}=zeros(2),
                                Kep::Float64=0.01,
                                cr_slope::Float64=4.5,
                                p_lim::Vector{Float64}=[1.0, 1.e6],
                                Nbins::Int64=48,
                                density_step::Bool=true,
                                n_blocks::Int64=70,
                                B_filename::String="")

            new(glass_file, output_file, n_blocks,
                v, B, U, B0, density_step, turb,
                CRp, CRe, Xcr, Kep, cr_slope, p_lim, Nbins, B_filename)

        end
    end

    function P_to_U(P::Real, rho::Real)
        return P/((2.0/3.0) * rho)
    end

    function U_to_P(U::Real, rho::Real)
        return (2.0/3.0) * rho * U
    end

    include("build_tube.jl")
    include("magnetic_field.jl")
    include("CRs.jl")
    include("main.jl")
    

end # module
