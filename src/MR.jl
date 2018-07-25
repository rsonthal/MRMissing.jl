module MR
using Distances, StatsBase

export Isomap, mds, procrustes, IOMR, MR_Missing

function Kmin(D,K)
    (n,n) = size(D)
    for i = 1:n
        M = findKmin(D[i,:], K+1)
        for j = 1:n
            if D[i,j] > M
                D[i,j] = Inf
            end
        end
    end
    
    for i = 1:n
        for j = 1:n
            if D[i,j] != Inf
                D[j,i] = D[i,j]
            end
        end
    end
    
    return D
end  

function findKmin(A,k)
    A = Set(A)
    for i = 1:k-1
        m = minimum(A)
        A = delete!(A,m)
        if isempty(A)
            return m
        end
    end
    
    return minimum(A)
end

function apsp(G)
    (n,n) = size(G)
    U = Inf*ones(size(G))
    for i = 1:n
        U[i,i] = 0
    end
    
    for i = 1:n
        for j = 1:n
            if G[i,j] != Inf
                U[i,j] = G[i,j]
            end
        end
    end
    
    for k = 1:n
        for i=1:n
            for j=1:i-1
                if U[i,j] > U[i,k] + U[k,j]
                    U[i,j] = U[i,k]+U[k,j]
                    U[j,i] = U[i,j]
                end
            end
        end
    end
    
    return U
end

function mds(D, d; center=true,align=true)
    n = size(D, 1)
    
    P = eye(n) - ones(n, n) / n
    K = -0.5*(P*(D.^2)*P') 
    K = (K+K')/2

    e, V = eig(K)
    idx = sortperm(e, rev=true)[1:d] 
    
    if (d > sum(e .> 0))
        println("Using some eigenvaectors for negative eigenvalues")
    end

    X =  V[:,idx]*diagm(sqrt.(e[idx])) 

    if center == true
        X = X .- mean(X,1) 
    end

    if align == true
        for i =1:d
            j = indmax(abs.(X[:,i]))
            X[:,i] .*= sign.(X[j,i])
        end
    end

    return X
end

function procrustes(X, Y)
    X₀ = X - ones(size(X)[1],1)*mean(X, 1)
    Y₀ = Y - ones(size(X)[1],1)*mean(Y, 1) 
    Z = X₀'*Y₀

    U,S,V = svd(Z, thin=true) 
    Q = V*U' 
    α = trace(Z*Q)/trace(Z) 
    Ya = α*(Y₀*Q) .+ mean(X, 1)
   
    return Ya
end

function Isomap(X,k,d)
    D = pairwise(Euclidean(1e-12),X',X')
    Dmin = Kmin(D,k)
    Dmani = apsp(Dmin)
    Xhat = mds(Dmani,d)
    return Xhat
end

function IOMR(Dpert)
    n = size(Dpert)[1]
    
    Dnew = copy(Dpert)
    
    for k = 1:n
        for i = 1:n
            for j = 1:i-1
                if i!=j && i!=k && j!=k
                    if Dnew[i,j] > Dnew[i,k]+Dnew[j,k]
                        Dnew[i,k] = Dnew[i,j]-Dnew[j,k]
                        Dnew[k,i] = Dnew[i,j]-Dnew[j,k]
                    end
                end
            end
        end
    end
    
    return Dnew-Dpert
end

function sumInf(x)
    n = length(x)
    s = 0.0
    for i = 1:n
        if (x[i] < Inf)
            s += x[i]
        end
    end
    
    return s
end

function distance(X)
    (n,m) = size(X)
    D = zeros(n,n)
    for i =1:n
        for j = 1:i-1
            b = (X[i,:]-X[j,:]).^2
            D[i,j] = sumInf(b)^0.5
            D[j,i] = (D[i,j])
        end
    end
    
    return D
end

function MR_Missing(X,Q,k)
    (n,m) = size(Q)
    for i = 1:n
        for j = 1:m
            if Q[i,j] == 0
                X[i,j] = Inf
            end
        end
    end
    
    Dmiss = distance(X)
    Dfixed = Dmiss + IOMR(Dmiss)
    Dmin = Kmin(Dfixed,k)
    
    return Dmin
end

end
