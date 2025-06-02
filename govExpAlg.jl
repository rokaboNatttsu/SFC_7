using StatsBase
using StatsPlots

β = 0.02
M, T = 1000, 100
G = zeros(M)
G_lst = sample(1:M, Int64(round(M/5)), replace=false)
Gm0 = 1.0
for i = G_lst
    G[i] = Gm0
end
G_sum = [sum(G_lst)*(1 + β)^(t-1) for t=1:T]

inn, outn = 10, 10
e = 0.1

function govExp(M, G, G_lst, t)
    # 契約終了処理
    outs = sample(G_lst, outn, replace=false)
    for m in outs
        G[m] = 0.0
        setdiff!(G_lst, [m])
    end
    # 新規契約処理
    for _ =1:inn
        m = G_lst[1]
        while m in G_lst
            m = sample(1:M)
        end
        push!(G_lst, m)
        G[m] = Gm0*(1.0 + β)^(t-1)
    end
    # 契約金額を経路依存的に決める
    ns = abs.(1 .+ e*randn(length(G_lst)))
    for (i, m) in enumerate(G_lst)
        G[m] *= ns[i]
    end
    G .*= G_sum[t]/sum(G)
end
function run(M, G, G_lst)
    for t = 1:T
        #push!(G, 0.0)
        govExp(M, G, G_lst, t)
    end
end
run(M, G, G_lst)
plot(1:length(G_lst), sort(G[G.>0.0], rev=true), xscale=:log10, yscale=:log10, xlabel="Rank", ylabel="Value")