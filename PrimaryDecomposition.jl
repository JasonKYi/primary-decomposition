"""
Primary Decomposition

1. Find minimal polynomial as factors
2. Find the generalised evectors
3. Change of basis
"""

import Base.*, Base.+, Base.-

using LinearAlgebra
using SymEngine

struct LinearMap{T <: Any, N <: Integer} <: AbstractArray{T, N} 
    length::N
    matrepr::Array{T, 2}
end

LinearMap{T}(l::N) where {T <: Number, N <: Integer} = LinearMap(l, Matrix{T}(undef, l, l))
LinearMap(l::N) where N <: Integer = LinearMap{Int64}(l)
LinearMap(M::Array{T, 2}) where T <: Number = LinearMap(size(M)[1], M)

*(x::Basic, A::Diagonal{Bool, Array{Bool, 1}}) = x * Array{Int}(A)
*(x::Basic, M::LinearMap{T, N}) where {T, N} = LinearMap(M.length, x * M.matrepr)
+(x::Basic, A::Array{T, 2}) where T = x * Array(size(A)[1] |> I) + A
-(x::Basic, A::Array{T, 2}) where T = x * Array(size(A)[1] |> I) - A
+(x::Int, A::Array{U, 2}) where U = x * Array(size(A)[1] |> I) + A
-(x::Int, A::Array{U, 2}) where U = x * Array(size(A)[1] |> I) - A

function charpoly(M::LinearMap, var::Basic, simpr::Bool = true)::Basic
    eq = det(var * I(M.length) - M.matrepr)
    simpr ? expand(eq) : eq
end

function terms(expr::Basic)::Array{SubString{String}, 1} # Maybe not so useful
    strexpr = "$expr"
    split(strexpr, r"\+ |\- ")
end

function distcomplex(z1::Vector{Basic}, z2::Vector{Basic})::Real    
    sum(map(abs, z1 - z2))
end

# Using Durand-Kerner method to find all complex roots
function durandkerner(M::LinearMap, var::Basic, maxiterations::Int = 15, 
            initvals::Complex = sqrt(2)im)::Vector{Basic}

    # Initiate the values
    d, iter, expr = M.length, maxiterations, charpoly(M, var)
    println(d)
    sol = map(x -> convert(Basic, x), Array(hcat([[0, initvals^i] for i = 1:d]...)'))

    # Begin the main loop
    while distcomplex(sol[:, end], sol[:, end - 1]) > 10e-6 && iter > 0
        newer = Vector{Basic}()
        for i = 1:d
            last = sol[:, end]
            soli = last[i] - expr(last[i]) / prod([last[i] - last[j] for j = 1:d if j != i])
            push!(newer, soli)
        end
        sol = hcat(sol, newer)
        iter -= 1
    end
    return sol[:, end]
end

function rroot(r::Basic)::Union{Basic, Nothing}
    abs(imag(r)) < 10e-6 ? real(r) : nothing
end

function rroots(r::Vector{Basic})::Vector{Basic}
    filter(x -> typeof(x) == Basic, map(rroot, r))
end

function getint(r::Vector{Basic})::Vector{Basic}
    map(x -> abs(x - round(x)) < 10e-2 ? round(x) : x, r)
end

function getlambda(var::Basic, expr::Basic)
    string("$var", "->", "$expr") |> Meta.parse |> eval
end

function findminipolyaux(M::LinearMap, roots::Vector{Basic}, var::Basic)
    maxdegree = M.length
    numroots = length(roots)
    candidates = Vector()
    for α ∈ Iterators.product([1:maxdegree for _ = 1:numroots]...)
        expr = [(var - roots[i])^(α[i]) for i = 1:numroots] |> prod |> expand
        push!(candidates, [α, getlambda(var, expr)])
    end
    return candidates
end

# function findminipoly(M::LinearMap, aux::Vector)
#     annihilators = Vector{Tuple{Int, Int}}()
#     for candidate ∈ aux
#         println(candidate[2](1))
#     end
# end

function main()
    @vars x
    # Companion matrix of x^4 + x^3 + 2x^2 + 1 which has no real roots
    M = [[0 0 0 -1]; 
         [1 0 0  0]; 
         [0 1 0 -2];
         [0 0 1 -1]]
    # Example matrix from lecture notes
    M2 = [[2 0 0]; [-1 -3 -1]; [-1 4 1]]
    M3 = [[0 -1 0]; [1 0 0]; [0 0 3]]
    
    Mw = LinearMap(M3)
    v = durandkerner(Mw, x) |> rroots |> getint
    # vnew = findminipolyaux(Mw, v, x)
    # p = findminipoly(Mw, vnew)
end

main()