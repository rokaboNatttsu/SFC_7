using StatsBase
using StatsPlots

growth_rate = 0.02
M, T = 100, 100
G, G_sum = zeros(M), zeros(T)
G_lst = sample(1:M, Int64(round(M/10)))
for i = G_lst
    G[i] = 1.0
end

inn, outn = 2, 2
e = 0.1

function govExp(growth_rate, M, Δgm, G, G_lst)
    d = 0.0
    outs = sample(G_lst, outn)
    for i in outs
        d += G[i]
        G[i] = 0.0
    end
    for _ =1:inn
        i = G_lst[1]
        while i in G_lst
            i = sample(1:M)
        end
        push!(G_lst, i)
        G[i] = Δgm
        d -= Δgm
    end
    ns = abs.(1 .+ e*randn(length(G_lst)))
    GS = sum(G)*(1.0 + growth_rate)
    for (i, gi) in enumerate(G_lst)
        G[gi] *= ns[i]
    end
    G .*= GS/sum(G)
end
function run(growth_rate, M, G, G_lst)
    for t = 1:T
        Δgm = 1.0 + 0.01*t
        #M += 1
        #push!(G, 0.0)
        govExp(growth_rate, M, Δgm, G, G_lst)
    end
end
run(growth_rate, M, G, G_lst)