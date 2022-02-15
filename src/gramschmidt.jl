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
