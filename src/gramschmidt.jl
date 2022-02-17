using LinearAlgebra
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

println("i\t∣A-Aₙ∣ \t∣gs(A).Q - gs(Aₙ).Q∣ \t∣qr(A).Q - qr(Aₙ).Q∣ (digits)")
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
end