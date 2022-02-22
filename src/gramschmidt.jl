using LinearAlgebra
using PrettyTables

A = Float64[ 1  2  1;
      0  2  3;
     -1  2 -2
    ]

function gramschmidt(A)
    n = size(A,2)
    Q = zeros(eltype(A), size(A))
    a₁ = A[:,1]
    Q[:, 1] = a₁/norm(a₁)
    for i in 2:n
        q̃ᵢ = A[:,i]
        for j in 1:i
            qⱼ = Q[:, j]
            q̃ᵢ .-= qⱼ'A[:,i]*qⱼ
        end
        if norm(q̃ᵢ) < 1e-8
            error("column $i is not linearly independent of basis")
        end
        Q[:,i] = q̃ᵢ/norm(q̃ᵢ)
    end
    return Q
end

Q = gramschmidt(A)

# check that Q is orthogonal
Q'Q

R = Q'A

# check that A = QR
A - Q*R

norm(A-Q*R)
# check that R is Upper Triangular
R - triu(R)

A′ = [A; 1 2 3]
Q′ = gramschmidt(A′)
Q′'A′ 

A₁ = Float64[ 1  1  1;
      1  1+1e-7  3;
     -1  -1 -2
    ]

Q₁ = gramschmidt(A₁)
A₂ = Float64[ 1  1  1;
1  1+2e-7  3;
-1  -1 -2
]

Q₂ = gramschmidt(A₂)
Q₂ - Q₁

# Exercise 1: Modify the gramschmidt function above to return a Tuple{Matrix, Vector}
# where the first element is the Q matrix as currently returned and the second element
# is the vector of vectors that are not linearly independent of the basis built by the
# gramschmidt process.

# we want a family of matrices that get worse and worse for our gramschmidt function
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

Abad  = nearlylindep(Float64[1,1,1,1], 3, 1e-2)
σ = 1e-14
Abadₙ = Abad .+ randn(size(Abad))*σ
norm(Abadₙ - Abad)
norm(gramschmidt(Abad) - gramschmidt(Abadₙ))

correctdigits(x) = floor(Int, -log10(norm(x, 2)))

results = NamedTuple[]
for i in 1:6
    σ = 1e-14
    Abad  = nearlylindep(ones(Float64, 6), 6, 10.0^-i)
    Abadₙ = Abad .+ randn(size(Abad))*σ
    ν  = correctdigits(Abadₙ - Abad)
    νq = correctdigits(gramschmidt(Abad) - gramschmidt(Abadₙ))
    Qʰ, Rʰ  = qr(Abad)
    Qʰₙ, Rʰₙ = qr(Abadₙ)
    νqʰ = correctdigits(Qʰ-Qʰₙ)
    println("$i\t$ν\t$νq\t$νqʰ")
    push!(results, (i=i,ϵ=10.0^-i,σ=σ, ν=ν, νq=νq, νqʰ=νqʰ))
end

pretty_table(results, header = ["i", "ϵ", "σ", "∣A-Aₙ∣", "∣gs(A).Q - gs(Aₙ).Q∣", "∣qr(A).Q - qr(Aₙ).Q∣ (digits)"])

#Exercise 2: use the following function to identify a different family of matrices that cause
# the gramschmidt process to show poor accuracy due to a loss of orthogonalality.

function compare_factorizer(n, Agenerator, σ, factorizer1=gramschmidt, factorizer2=qr)
    results = NamedTuple[]
    for i in 1:n
        Abad  = Agenerator(i)
        Abadₙ = Abad .+ randn(size(Abad))*σ
        ν  = correctdigits(Abadₙ - Abad)
        νq = correctdigits(factorizer1(Abad) - factorizer1(Abadₙ))
        Qʰ, Rʰ  = factorizer2(Abad)
        Qʰₙ, Rʰₙ = factorizer2(Abadₙ)
        νqʰ = correctdigits(Qʰ-Qʰₙ)
        push!(results, (i=i,ϵ=10.0^-i,σ=σ, normAbad=norm(Abad), κ=cond(Abad), ν=ν, νq=νq, νqʰ=νqʰ))
    end
    pretty_table(results, header = ["i", "ϵ", "σ", "|A|", "κ", "∣A-Aₙ∣", "∣gs(A).Q - gs(Aₙ).Q∣", "∣qr(A).Q - qr(Aₙ).Q∣ (digits)"])
end

compare_factorizer(7, i->nearlylindep(ones(Float64, 6), 6, 10.0^-i), 1e-14, gramschmidt, qr)

# This function makes matrices with a known norm and condition number using the fact
# that you can read the 2-norm and 2-κ off the singular value decomposition of a matrix.
function badsvd(n,k)
    Q = qr(randn(n,k)).Q
    Q*diagm([10.0^-(i-1) for i in 1:k])*Q'
end


A = badsvd(10,10)
norm(A) # ≈ 1
cond(A) # ≈ 1e-9

compare_factorizer(7, i->badsvd(1+i,1+i), 1e-10)


function nearlylindep₂(a₁::AbstractVector{T}, n::Int, ϵ::Real) where T
    m = length(a₁)
    e(i::Int) = begin # standard basis vector in m-dims
        eᵢ = zeros(T, m)
        eᵢ[i] = 1
        return eᵢ
    end
    A = zeros(T, m, n)
    A[:, 1] = a₁
    for i in 2:n
        aᵢ = sum(A, dims=2)/i + ϵ * e(i) + ϵ*e(mod(i+1, 1:n))
        A[:, i] = aᵢ
    end
    return A
end

compare_factorizer(7, i->nearlylindep₂(ones(Float64, 6), 5, 10.0^-i), 1e-14)