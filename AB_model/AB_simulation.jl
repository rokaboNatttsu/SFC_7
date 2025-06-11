using StatsPlots
using StatsBase

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

EMP = zeros(Int64, L,TIME)
G_calc_item, G_potential = ones(N), ones(N)
Cgrs = zeros(O)

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
            new_o = searchsortedlast(r, rand())
            c[new_o,l,t] = Cs[l]/p[new_o,t]
        end
    end
end

function w_and_W_func(t)
    for l=1:L
        if EMP[l,t-1] > 0
            if rand() < ζ3
                continue
            else
                EMP[l,t] = EMP[l,t-1]
                w[l,t] = (1-ζ1)*w[l,t-1]+ζ1*w[l,t-1]*(1+ζ2*abs(randn()))
                W[l,EMP[l,t],t] = w[l,t]
            end
        else
            offers = [Int64(max(0, (u[o,t-1]*k[o,t-1]-A[o,t]*sum(w[:,o,t-1]>0))/A[o,t])) for o=1:O]
            prob = deepcopy(offers)
            for o=2:O
                prob[o,t] += prob[o,t-1]
            end
            prob /= prob[end]
            appli = [[] for _ = 1:O]
            for l=1:L
                if EMP[l,t-1] == 0
                    o = searchsortedlast(prob, rand())
                    push!(appli[o], l)
                end
            end
            for o=1:O
                if offers[o] > 0 & length(appli[o]) > 0
                    much = collect(Set(sample(appli[o], offers[o])))
                    for x=1:length(much)
                        l = much[x]
                        EMP[l,t] = o
                        w[l,t] = (1-ζ1)*w[l,t-1]+ζ1*w[l,t-1]*(1+ζ2*abs(randn()))
                        W[l,EMP[l,t],t] = w[l,t]
                    end
                end
            end
        end
    end
    for l=1:L
        if EMP[l,t]==0
            w[l,t] = (1-ζ1)*w[l,t-1]+ζ1*w[l,t-1]*(1-ζ3*abs(randn()))
        end
    end
end

function Lh_func(t)
    Lhs = ϵ1*NLh[:,t] + ϵ2*sum(C[:,:,t], dims=1)
    for l=1:L
        alreadyexist = false
        for n=1:N
            if Lh[l,n,t-1] > 0.0
                Lh[l,n,t] = Lhs[l]
                alreadyexist = true
                break
            end
        end
        if ! alreadyexist
            prob = zeros(N)
            prob[1] = NWb[1,t]
            for n=2:N
                prob[n] = NWb[n,t] + prob[t-1]
            end
            prob /= prob[end]
            n = searchsortedlast(prob, rand())
            Lh[l,n,t] = Lhs[n]
        end
    end
end

function ΔLf_and_Lf_func(t)
    ΔLfs = max.(-sum(Lf[:,:,t-1], dims=2), (λ3+λ4*((P[:,t]-Pf[:,t])./(sum(Eh[:,:,t-1], dims=2)+sum(Eb[:,:,t-1], dims=2))-rL)).*(I[:,t]+sum(W[:,:,t], dims=1)+Tv[:,t]+Tc[:,t]+rL*sum(Lf[:,:,t-1], dims=2)-ϕ*sum(Mf[:,:,t-1], dims=1)))
    for o=1:O
        alreadyexist = false
        for n=1:N
            if Lf[o,n,t-1] > 0.0
                Lf[o,n,t] = Lf[o,n,t-1] + ΔLfs[o]
                alreadyexist = true
                break
            end
        end
        if ! alreadyexist
            prob = zeros(N)
            prob[1] = NWb[1,t]
            for n=2:N
                prob[n] = NWb[n,t] + prob[t-1]
            end
            prob /= prob[end]
            n = searchsortedlast(prob, rand())
            Lf[o,n,t] = ΔLfs[o]
        end
    end
end

function household_portfolio_func(t)
    # 金融資産額の推定
    Vs = sum(Mh[:,:,t-1], dims=1)+sum(Eh[:,:,t-1], dims=1)+sum(Fh[:,:,t-1], dims=1)+NLh[:,t]
    # ポートフォリオ配分先の確率の重みの共通部分を計算
    index = (Pf[:,t] - I[:,t]).*(P[:,t]-Pf[:,t])./(sum(Eh[:,:,t-1], dims=2)+sum(Eb[:,:,t-1], dims=2))  # 利益☓配当率
    append!(index, (rL.*L[:,t-1]+sum(Pb[:,:,t], dims=2)).*sum(S[:,:,t], dims=1)./sum(Fh[:,:,t-1], dims=2))
    index /= sum(index)
    # ポートフォリオ配分先の数Int64(x[l])を決めるための準備。
    x = 0.5*Vs./mean(sum(C[:,:,t], dims=1))
    # 家計lのポートフォリオ計算
    for l=1:L
        # ポートフォリオ配分先に選ぶ確率の重み付けを決める
        tmp = Eh[:,l,t-1]
        append!(tmp, Fh[:,l,t-1])
        prob = index + tmp/(sum(Eh[:,l,t-1])+sum(Fh[:,l,t-1]))
        prob /= sum(prob)
        # ポートフォリオ配分先のリストを作る
        lst = sample(prob, Int64(x[l]))
        if length(lst)==0
            continue
        end
        # 収益率から、資産に占める株式の割合の目標を計算し、保有額を決める
        EF_volume = Vs[l]*(λ1 + λ2*(sum(Ph[l,:,t])+sum(S[l,:,t]))/(sum(Eh[:,l,t-1])+sum(F[:,l,t-1])))
        # EF_volume/length(lst)を単位としてlstの企業または銀行の株の保有割合を決める
        for on in lst
            if on <= O
                Eh[on,l,t] += EF_volume/length(lst)
            else
                Fh[on-O,l,t] += EF_volume/length(lst)
            end
        end
    end
end

function one_season(TIMERANGE)
    for t=TIMERANGE
        p[:,t] = λp*(1+ν1+ν2*sum(Lf[:,:,t-1], dims=2)./(sum(C[:,:,t-1], dims=2)+G[:,t-1])).*(sum(W[:,:,t-1], dims=1)+Tv[:,t-1]+Tc[:,t-1]+δ*k[:,t-1])./(uT*γ1*k[:,t-1])+(1-λp)*ν3*(p[:,t-1].-mean(p[:,t-1]))
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
        i[:,t] = δ*k[:,t-1]+(u[:,t-1]-uT).*γ2.*k[:,t-1]+γ3*(sum(Mf[:,:,t-1], dims=1)-sum(Lf[:,:,t-1], dims=2))./p[:,t-1]
        I[:,t] = p[:,t].*i[:,t]
        k[:,t] = (1-δ).*k[:,t-1]+i[:,t]
        K[:,t] = p[:,t].*k[:,t]
        w_and_W_func(t)
        P[:,t] = sum(C[:,:,t], dims=2)+G[:,t]+I[:,t]-sum(W[:,:,t], dims=1)-Tc[:,t]-Tv[:,t]-rL*sum(Lf[:,:,t-1], dims=2)
        for o=1:O
            Ph[:,o,t] = max(0.0, θ1*(P[o,t]-I[o,t])+θ2*(sum(Mf[:,o,t-1])-sum(Lf[o,:,t-1]))).*eh[o,:,t]./e[o,t]
            Pb[:,o,t] = max(0.0, θ1*(P[o,t]-I[o,t])+θ2*(sum(Mf[:,o,t-1])-sum(Lf[o,:,t-1]))).*eb[o,:,t]./e[o,t]
        end
        Pf[:,t] = P[:,t] - sum(Ph[:,:,t], dims=1) - sum(Pb[:,:,t], dims=1)
        for n=1:N
            S[:,n,t] = (θ3*(rL*L[n,t-1]+sum(Pb[n,:,t]))+θ4*sum(Eb[:,n,t-1])).*fb[:,n,t-1]./f[n,t-1]
        end
        NLh[:,t] = -sum(C[:,:,t], dims=1) + sum(W[:,:,t], dims=2)-Ti[:,t]-Ta[:,t]-rL*sum(Lh[:,:,t-1], dims=2)+sum(Ph[:,:,t], dims=2)+sum(S[:,:,t], dims=2)
        NLf[:,t] = -I[:,t] + Pf[:,t]
        NLb[:,t] = rL*L[:,t-1] + sum(Pb[:,:,t], dims=1) - sum(S[:,:,t], dims=1)
        NLg[t] = -sum(G[:,t])+sum(Ti[:,t])+sum(Ta[:,t])+sum(Tv[:,t])+sum(Tc[:,t])
        Lh_func(t)
        ΔLh[:,:,t] = Lh[:,:,t] - Lh[:,:,t-1]
        ΔLf_and_Lf_func(t)
        L[:,t] = sum(Lh[:,:,t], dims=1) + sum(Lf[:,:,t], dims=1)
        ΔL[:,t] = L[:,t] - L[:t-1]
        Δe[:,t] = 1/p[:,t-1]*(1-λ3-λ4*((P[:,t]-Pf[:,t])./(sum(Eh[:,:,t-1], dims=2)+sum(Eb[:,:,t-1], dims=2))-rL)).*(I[:,t]+sum(W[:,:,t], dims=1)+Tv[:,t]+Tc[:,t]+rL*sum(Lf[:,:,t-1], dims=2)-ϕ*sum(Mf[:,:,t-1], dims=1))
        E[:,t] = E[:,t-1] + pe[:,t-1].*Δe[:,t]
        e[:,t] = e[:t-1] + Δe[:,t]
        household_portfolio_func(t)
    end
end