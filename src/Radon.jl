module Radon
export radon!

function radon(im::AbstractArray{T, 2}) where {T<:AbstractFloat}
    Nx, Ny = size(im)
    radonim = similar(im, min(Nx, Ny), max(Nx, Ny))
    radon!(radonim, im)
    radonim
end

function radon!(radonim::AbstractArray{T, 2}, im::AbstractArray{T, 2}) where {T<:AbstractFloat}
    Nφ, Ns = size(radonim)
    Nx, Ny = size(im)
    Ldiag = hypot(Nx, Ny)
    K = round(Int, Ldiag)

    dLz = Ldiag/(K-1)
    dLs = Ldiag/(Ns-1)
    dφ = π/Nφ
    for is = 1:Ns, iφ = 1:Nφ
        φ = (iφ-1)*dφ
        s = -Ldiag/2 + (is-1)*dLs
        tx = cos(φ)
        ty = sin(φ)

        radonim[iφ, is] = zero(T)
        for k=1:K
            z = -Ldiag/2 + (k-1)*dLz
            x = round(Int, Nx/2 + z*ty + s*tx)
            y = round(Int, Ny/2 - z*tx + s*ty)

            if 1 <= x <= Nx && 1 <= y <= Ny
                radonim[iφ, is] += im[x, y]*dLz
            end
        end
    end
    radonim
end

function backproject(radonim::AbstractArray{T, 2}) where {T<:AbstractFloat}
    Nφ, Ns = size(radonim)
    im = similar(radonim, Ns, Ns)
    backproject!(im, radonim)
    im
end

function backproject!(im::AbstractArray{T, 2}, radonim::AbstractArray{T, 2}) where {T<:AbstractFloat}
    Nφ, Ns = size(radonim)
    Nx, Ny = size(im)
    Ldiag = hypot(Nx, Ny)
    K = round(Int, Ldiag)

    dLz = Ldiag/(K-1)
    dLs = Ldiag/(Ns-1)
    dφ = π/Nφ

    im[:] = 0.
    for is = 1:Ns, iφ = 1:Nφ
        φ = (iφ-1)*dφ
        s = -Ldiag/2 + (is-1)*dLs
        tx = cos(φ)
        ty = sin(φ)

        for k=1:K
            z = -Ldiag/2 + (k-1)*dLz
            x = round(Int, Nx/2 + z*ty + s*tx)
            y = round(Int, Ny/2 - z*tx + s*ty)

            if 1 <= x <= Nx && 1 <= y <= Ny
                im[x, y] += radonim[iφ, is]*dLz
            end
        end
    end
    im
end

end
