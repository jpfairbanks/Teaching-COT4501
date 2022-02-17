using LinearAlgebra

𝐯(z) = begin
    v = -z
    v[1] += norm(z)
    return v
end

proj(v, z) = v'z*v / v'v
reflector(v::AbstractVector) = I - 2(v*v')/(v'v)
hh(z::AbstractVector) = reflector(𝐯(z))
angle₁(x) = (rad2deg(acos(x[1]/norm(x))))

P₁ = hh(ones(2))
P₁*ones(2)

P₁ = hh(ones(4))
P₁*ones(4)

𝟙 = ones(2)
x = 𝟙 - proj(𝐯(𝟙), 𝟙)
angle₁(x)

proj(𝐯(ones(2)), ones(2))

y = 2*[cos(deg2rad(60)),sin(deg2rad(60))]
angle₁(y - proj(𝐯(y), y))

"""
qrfact(A)

QR factorization by Householder reflections. Returns Q and R.
"""
function qrfact(A)
    m,n = size(A)
    Qt = Matrix(Diagonal(ones(m)))
    R = float(copy(A))
    for k in 1:n
      z = R[k:m,k]
      v = [ -sign(z[1])*norm(z) - z[1]; -z[2:end] ]
      nrmv = norm(v)
      if nrmv < eps() continue; end  # skip this iteration
      v = v / nrmv;                  # simplifies other formulas
      # Apply the reflection to each relevant column of A and Q
      for j in 1:n
        R[k:m,j] -= v*( 2*(v'*R[k:m,j]) )
      end
      for j in 1:m
        Qt[k:m,j] -= v*( 2*(v'*Qt[k:m,j]) )
      end
    end
    
    return Qt',triu(R)
end

function nearlylindep(a₁::AbstractVector{T}, n::Int, ϵ::Real) where T
    m = length(a₁)
    e(i::Int) = begin # standard basis vector in m-dims
        eᵢ = zeros(T, m)
        eᵢ[i] = 1
        return eᵢ
    end
    A = zeros(T, m, n)
    A[:, 1] = a₁
    for i in 2:n
        aᵢ = sum(A, dims=2)/i + ϵ * e(i)
        A[:, i] = aᵢ
    end
    return A
end

println("i\t∣A-Aₙ∣ \t∣hh(A).Q - hh(Aₙ).Q∣ \t∣qr(A).Q - qr(Aₙ).Q∣  (digits)")
for i in 1:6
    σ = 1e-14
    Abad  = nearlylindep(ones(Float64, 6), 6, 10.0^-i)
    Abadₙ = Abad .+ randn(size(Abad))*σ
    ν  = correctdigits(Abadₙ - Abad)
    νq = correctdigits(qrfact(Abad)[1] - qrfact(Abadₙ)[1])
    Qʰ, Rʰ  = qr(Abad)
    Qʰₙ, Rʰₙ = qr(Abadₙ)
    νqʰ = correctdigits(Qʰ-Qʰₙ)
    println("$i\t$ν\t$νq\t$νqʰ")
end

# read through the docs for qr for info on pivoting and the extended API