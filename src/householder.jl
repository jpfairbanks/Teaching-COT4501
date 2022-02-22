using LinearAlgebra

"""    𝐯(z::Vector)

compute the vector v from the definition of the Householder reflectors.
Type as `\\bfv<TAB>` for "BoldFace v".
"""
𝐯(z) = begin
    v = -z
    v[1] += norm(z)
    return v
end

"""    proj(v::AbstractVector, z::AbstractVector)

The vector that takes `z` to the space orthogonal to v
"""
proj(v::AbstractVector, z::AbstractVector) = v*v'z / v'v

"""    reflector(v::AbstractVector)

Construct the reflector P from the vector v in Householder definition
"""
reflector(v::AbstractVector) = I - 2(v*v')/(v'v)

"""    hh(z::AbstractVector)

use `𝐯` and `reflector` to construct `P` the Householder reflector that will zero out
the entries of `z` except for `z₁`.
"""
hh(z::AbstractVector) = reflector(𝐯(z))

"""    angle₁(x::AbstractVector)

Compute the angle to the e₁ line using the cos⁻¹(x[1]/norm(x)) formula. Only works in 2-D.
"""
angle₁(x) = (rad2deg(acos(x[1]/norm(x)))) * sign(x[2])

P₁ = hh(ones(2))
P₁*ones(2)

P₁ = hh(ones(4))
P₁*ones(4)

SE = [1, -1.0]/2
𝟙 = ones(2)/2
SE'𝟙
proj(SE, 𝟙)
angle₁(𝟙)
x = 𝟙 - proj(SE, 𝟙)
angle₁(x)
angle₁(𝟙 - 2proj(SE, 𝟙))

u = [1, 2.0]
u = u / norm(u)
angle₁(u)
u - proj(𝐯(SE), u)
angle₁(u - proj(SE, u))
norm(u - 2proj(SE, u))
angle₁(u - 2proj(SE, u))

# Exercise 1: Draw out the picture of the examples above in the style of textbook Figure 3.3.

# Exercise 2: Pick two additional vectors for defining your plane of reflection and two additional vectors to reflect.
# Draw out the diagram and check that they angles are correct.

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

using PrettyTables

function compare_factorizer(n, Agenerator, σ, factorizer1=qrfact, factorizer2=qr)
    results = NamedTuple[]
    for i in 1:n
        Abad  = Agenerator(i)
        Abadₙ = Abad .+ randn(size(Abad))*σ
        ν  = correctdigits(Abadₙ - Abad)
        νq = correctdigits(factorizer1(Abad)[1] - factorizer1(Abadₙ)[1])
        Qʰ, Rʰ  = factorizer2(Abad)
        Qʰₙ, Rʰₙ = factorizer2(Abadₙ)
        νqʰ = correctdigits(Qʰ-Qʰₙ)
        push!(results, (i=i,ϵ=10.0^-i,σ=σ, normAbad=norm(Abad), κ=cond(Abad), ν=ν, νq=νq, νqʰ=νqʰ))
    end
    pretty_table(results, header = ["i", "ϵ", "σ", "|A|", "κ", "∣A-Aₙ∣", "∣hh(A).Q - hh(Aₙ).Q∣", "∣qr(A).Q - qr(Aₙ).Q∣ (digits)"])
end


compare_factorizer(7, i->nearlylindep(ones(Float64, 6), 6, 10.0^-i), 1e-14, qrfact, qr)
# read through the docs for qr for info on pivoting and the extended API

function qrfact₂(A)
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
      R[k:m,:] .-= v*( 2*(v'*R[k:m,:]) )
      Qt[k:m,:] .-= v*( 2*(v'*Qt[k:m,:]) )
    end
    return Qt',triu(R)
end

A₃ = nearlylindep(ones(Float64, 6), 6, 10.0^-3)
norm(qrfact₂(A₃)[1] - qrfact(A₃)[1])

# Exercise 3: Benchmark the qrfact and qrfact₂ implementations with `@time` as the problem size increases. Remember to run the code twice and take the last time because of precompilation.
# Plot the runtime and memory allocations as a function of problem size. Which one is faster? Does memory use drive runtime?

@time qrfact(A₃);
@time qrfact₂(A₃);