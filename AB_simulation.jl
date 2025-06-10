using StatsPlots

L, O, N, TIME = 10000, 100, 10, 100

C, c, G, g, I. i = zeros(O, L, TIME), zeros(O, L, TIME), zeros(O, TIME), zeros(O, TIME), zeros(O, TIME), zeros(O, TIME)
W, Ti, Ta, Tv, Tc = zeros(L, O, TIME), zeros(L, TIME), zeros(L, TIME), zeros(O, TIME), zeros(O, TIME)
Ph, P, Pf, Pb = zeros(L, O, TIME), zeros(O, TIME), zeros(O, TIME), zeros(N, O, TIME)
S = zeros(L, N, TIME)
NLh, NLf, NLb, NLg = zeros(L, TIME), zeros(O, TIME), zeros(N, TIME), zeros(TIME)
p, pe, pf = zeros(O, TIME), zeros(O, TIME), zeros(N, TIME)
eh, e, eb, Δeh, Δe, Δeb = zeros(O, L, TIME), zeros(O, TIME), zeros(O, N, TIME), zeros(O, L, TIME), zeros(O, TIME), zeros(O, N, TIME)
Eh, E, Eb = zeros(O, L, TIME), zeros(O, TIME), zeros(O, N, TIME)
fh, Δfh = zeros(N, L, TIME), zeros(N, L, TIME)
Fh = zeros(N, L, TIME)
Mh, Mf, M, ΔMh, ΔMf, ΔM = zeros(N, L, TIME), zeros(N, O, TIME), zeros(N, TIME), zeros(N, L, TIME), zeros(N, O, TIME), zeros(N, TIME)
Lh, Lf, L, ΔLh, ΔLf, ΔL = zeros(L, N, TIME), zeros(O, N, TIME), zeros(N, TIME), zeros(L, N, TIME), zeros(O, N, TIME), zeros(N, TIME)
H, ΔH = zeros(N, TIME), zeros(N, TIME)
K, k = zeros(O, TIME), zeros(O, TIME)
DE, DF = zeros(TIME), zeros(TIME)
NWh, NWf, NWb, NWg, NWg = zeros(L, TIME), zeros(O, TIME), zeros(N, TIME), zeros(TIME), zeros(TIME)

rL = 0.02
G0 = 1.0*O

G_calc_item, G_potential = ones(N), ones(N)

function G_func(t)
    Gsum = G0*(1+β1)^(t-1)
    G_calc_item .*= 1.0 .+ β2.*(-1.0 .+ β3*randn(N))
    G_calc_item = max.(1.0, G_calc_item)
    G_potential = max.(0.0, G_calc_item .- 2.0)
    G[:,t] = G_potential.*Gsum/sum(G_potential)
end

function c_func(t)
    trans_or_stay = rand(O) .< α3*(p[:,t].-mean(p[:,t]))/p[:,t]
    r = zeros(O)
    r[1] = sum(c[:,l,t-1])+α4*sum(c[:,:,t-1])
    for o=2:O
        r[o] = sum(c[:,l,t-1]) + α4*sum(c[:,:,t-1]) + r[o-1]
    end
    r /= r[end]
    Cs = α1*(sum(W[:,:,t-1], dims=1)-Ta[:,t]-Ti[:,t]-rL*sum(Lh[:,:,t-1], dims=2)+sum(Ph[:,:,t-1], dims=1)+sum(S[:,:,t-1], dims=1))+α2*(sum(Eh[:,:,t-1], dims=1)+sum(Mh[:,:,t-1], dims=1)-sum(Lh[:,:,t-1], dims=1))
    for (o, trans) in enumerate(trans_or_stay)
        if ! trans
            for l=1:L
                if c[o,l,t-1] > 0.0
                    c[o,l,t] = Cs[l]/p[o,t]
                    break
                end
            end
        else
            tmp = rand()
            for new_o=1:O
                if tmp <= r[new_o]
                    c[new_o,l,t] = Cs[l]/p[new_o,t]
                    break
                end
            end
        end
    end
end

function one_season(TIMERANGE)
    for t=TIMERANGE
        p[:,t] = λp*(1+ν1+ν2*sum(Lf[:,:,t-1], dims=2)./(sum(C[:,:,t-1], dims=2)+G[:,t-1]))*(sum(W[:,:,t-1], dims=1)+Tv[:,t-1]+Tc[:,t-1]+δ*k[:,t-1])/(uT*γ1*k[:,t-1])+(1-λp)*ν3*(p[:,t-1].-mean(p[:,t-1]))
        Ti[:,t] = τ1*(W[:,o,t-1]+P[:,o,t-1]+S[:,n,t-1])
        Ta[:,t] = τ2*(sum(Eh[:,:,t-1], dims=1)+sum(Mh[:,:,t-1], dims=1)-sum(Lh[:,:,t-1], dims=2))
        Tv[:,t] = τ3*(sum(C[:,:,t-1], dims=2)+I[:,t-1]+G[:,t-1])
        Tc[:,t] = τ4*(sum(C[:,:,t-1], dims=2)+G[:,t-1]+I[:,t-1]-sum(W[:,:,t-1], dims=1)-Tv[:,t-1])
        G_func(t)
        g[:,t] = G[:,t]./p[:,t]
        c_func(t)
        for o=1:O
            C[o,:,t] = p[o,t]*c[o,:,t]
        end
        A[:,t] = A[:,t-1]*(1+μ1.+μ2*i[:,t-1]./k[:,t-1])
        u[:,t] = (sum(c[:,:,t], dims=2)+i[:,t]+g[:,t])./(γ1*k[:,t-1])
        i[:,t] = 
    end
end