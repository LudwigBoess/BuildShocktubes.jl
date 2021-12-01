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

