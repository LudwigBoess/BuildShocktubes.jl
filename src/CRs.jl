

"""
    find_init_norm(pressure::T, slope::T, bound_low::T, bound_up::T) where T

Find norm of first bin for a given total pressure.
"""
function find_init_norm(pressure::Real, slope::Real, bound_low::Real, bound_up::Real)

    c_light    = 2.9979e10
    slope_soft = 1.e-6
    
    norm_cnst = 3.0*pressure / ( 4π * c_light * bound_low^4 )
    init_norm = norm_cnst * ( 4.0 - slope ) / ( (bound_up/bound_low)^(4.0 - slope) - 1.0 )

    if ( 4.0 - slope_soft ) < slope < ( 4.0 + slope_soft )

        slope_var = (slope - 4.0)/slope_soft
        init_norm2 = norm_cnst / log(bound_up/bound_low)

        if slope_var != 0.0
            return init_norm * slope_var + init_norm2 * ( 1.0 - slope_var )
        else
            return init_norm2
        end
    end

    return init_norm
end


"""
    init_powerlaw( CR::AbstractCRSpectrum, 
                   init_norm::T, init_slope::T, bounds::Array{T}, density::T) where T

Initialize a single powerlaw starting from an initial norm `init_norm`.
"""
function init_powerlaw( init_norm::Real, init_slope::Real, bounds::Vector{<:Real})

    Nbins = length(bounds)-1
    Norm  = Vector{Float64}(undef, Nbins)

    # initialize first bin by hand
    Norm[1]  = init_norm

    # loop over the remaining bins
    @inbounds for Nbin = 2:Nbins
        # Calculate the consecutive norms
        Norm[Nbin] = Norm[Nbin-1] * (bounds[Nbin]/bounds[Nbin-1])^(-init_slope) 
    end

    return Float32.(log10.(Norm))
end

"""
    set_bounds(p_lim::Vector{Float64}, Nbins::Int64)

Set bin boundaries of spectrum.
"""
function set_bounds(p_lim::Vector{Float64}, Nbins::Int64)
    bounds    = Vector{Float64}(undef, Nbins+1)
    bin_width = log10( p_lim[2] / p_lim[2] ) / Nbins

    @inbounds for i = 1:Nbins+1
        bounds[i] = p_lim[2] * 10.0^( (i-1) * bin_width )
    end
    return bounds
end