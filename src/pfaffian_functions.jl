function indicator_vector(a, k)
    v = zeros(Int64, k)
    for i in a
        v[i] = 1
    end
    return v
end

function isaDeltaMatroid(Bs)
    for A in Bs, B in Bs
        sAB = symdiff(A,B)
        for a in sAB
            if any([sort!(symdiff(A,[a,b])) in Bs for b in sAB if b != a])
                continue
            else
                println("A = ", A)
                println("B = ", B)
                println("a = ", a)
                println("b = ", b)
                return false                
            end
        end
    end
    return true
end

function set2str(F)
    return join([string(i) for i in F])
end

function coordinate_matrix(Q, n, R, x, xdict)
    ijs = [[i,j] for i in 1:n, j in 1:n if (i<j && setdiff(1:n,[i,j]) in Q) ]
    #ijsn = [[i,j] for i in 1:n, j in 1:n if i<j ]
    S = MatrixSpace(R, n, n)
    X = S()
    for i in 1:n, j in 1:n
        ij=[i,j]
        if ij in ijs
            X[ij[1],ij[2]] =  xdict[ij]
            X[ij[2],ij[1]] =  -xdict[ij]
        end
    end
    return X
end


function pfaffian_coordinate(X, S, n)
    nS = [i for i in 1:n if !(i in S)]
    return pfaffian(X[nS, nS])
end

function even_Delta_matroid_stratum(Q, n, F)
    even = [a for a in collect(powerset(1:n)) if length(a)%2 == 0]
    ijs = [[i,j] for i in 1:n, j in 1:n if (i<j && setdiff(1:n,[i,j]) in Q) ]
    R, x = PolynomialRing(F, :"x" => ijs); 

    xdict = Dict{Vector{Int}, MPolyElem}([ijs[i] => x[i] for i in 1:length(ijs)]);
    
    X = coordinate_matrix(Q, n, R, x, xdict)
    Igens = unique!([pfaffian_coordinate(X, S, n) for S in even if !(S in Q)])
    Sgens = unique!([pfaffian_coordinate(X, S, n) for S in Q])
    return (Igens, Sgens, R, x)
    
end


#function delta_matroid_to_reduced_expression(Q, n, F, k)
#    if !(1:n in Q)
#        return "1:n not in even Delta matroid"
#    end
#    RQ = even_Delta_matroid_stratum(Q, n, F)
#    R = parent(RQ[1][1])
    
#    I = reduce_ideal_full(RQ[1], leq_n_terms(RQ[2],k), R, gens(R), false)
#    return (I[1], I[2])
#end

#function twist(Q, S)
#    return sort!([sort!(symdiff(q,S)) for q in Q])
#end
