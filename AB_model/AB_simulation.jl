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
    global G
    Gsum = G0*(1+β1)^(t-1)
    G_calc_item .*= 1.0 .+ β2.*(-1.0 .+ β3*randn(N))
    G_calc_item = max.(1.0, G_calc_item)
    G_potential = max.(0.0, G_calc_item .- 2.0)
    G[:,t] = G_potential.*Gsum/sum(G_potential)
end

function c_func(t)
    global c
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
    global w, W, EMP
    for l=1:L
        if EMP[l,t-1] > 0 # 前期就業していた人
            if rand() < ζ3 # 今期失業していた場合
                continue
            else # 今期も同じ企業で就業する場合
                EMP[l,t] = EMP[l,t-1]
                w[l,t] = (1-ζ1)*w[l,t-1]+ζ1*w[l,t-1]*(1+ζ2*abs(randn()))
                W[l,EMP[l,t],t] = w[l,t]
            end
        else # 前期失業していた人
            # 求人数を作る
            offers = [Int64(max(0, (u[o,t-1]*k[o,t-1]-A[o,t]*sum(w[:,o,t-1]>0))/A[o,t])) for o=1:O]
            # 応募確率を作る
            prob = deepcopy(offers)
            for o=2:O
                prob[o,t] += prob[o,t-1]
            end
            prob /= prob[end]
            # 応募先を割り振る
            appli = [[] for _ = 1:O]
            for l=1:L
                if EMP[l,t-1] == 0
                    o = searchsortedlast(prob, rand())
                    push!(appli[o], l)
                end
            end
            # マッチング
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
    # 従業員が0にならないように対策
    s = Set(EMP[:,t])
    UE = findall(x -> x == 0, EMP[:,t])
    for o=1:O
        if !(o in s)
            l = sample(UE)
            setdiff!(UE, [l])
            EMP[l,t] = o
            w[l,t] = (1-ζ1)*w[l,t-1]+ζ1*w[l,t-1]*(1+ζ2*abs(randn()))
            W[l,EMP[l,t],t] = w[l,t]
        end
    end
    # 失業者の要求賃金を下げる
    for l=1:L
        if EMP[l,t]==0
            w[l,t] = (1-ζ1)*w[l,t-1]+ζ1*w[l,t-1]*(1-ζ3*abs(randn()))
        end
    end
end

function Lh_func(t)
    global Lh
    Lhs = ϵ1*NLh[:,t] + ϵ2*sum(C[:,:,t], dims=1)
    for l=1:L
        for n=1:N
            if Mh[n,l,t-1] > 0.0
                Lh[l,n,t] = Lhs[l]
                break
            end
        end
    end
end

function ΔLf_and_Lf_func(t)
    global ΔLf, Lf
    ΔLfs = max.(-sum(Lf[:,:,t-1], dims=2), (λ3+λ4*((P[:,t]-Pf[:,t])./(sum(Eh[:,:,t-1], dims=2)+sum(Eb[:,:,t-1], dims=2))-rL)).*(I[:,t]+sum(W[:,:,t], dims=1)+Tv[:,t]+Tc[:,t]+rL*sum(Lf[:,:,t-1], dims=2)-ϕ*sum(Mf[:,:,t-1], dims=1)))
    for o=1:O
        for n=1:N
            if Mf[n,o,t-1] > 0.0
                if Lf[o,n,t-1] > 0.0
                    Lf[o,n,t] = Lf[o,n,t-1] + ΔLfs[o]
                    ΔLf[o,n,t] = ΔLfs[o]
                else
                    for n2=1:N
                        if Lf[o,n2,t-1] > 0.0
                            Lf[o,n2,t] = Lf[o,n,t-1] + ΔLfs[o]
                            ΔLf[o,n,t-1] = -Lf[o,n,t-1]
                            ΔLf[o,n2,t] = Lf[o,n,t-1] + ΔLfs[o]
                            break
                        end
                    end
                end
                break
            end
        end
    end
end

function household_portfolio_func(t)
    global Eh, Fh
    # 金融資産額の推定
    Vs = sum(Mh[:,:,t-1], dims=1)+sum(Eh[:,:,t-1], dims=1)+sum(Fh[:,:,t-1], dims=1)+NLh[:,t]
    # ポートフォリオ配分先の確率の重みの共通部分を計算
    index = ψ1*(Pf[:,t] - I[:,t])+ψ2*(sum(Mf[:,:,t-1], dims=1)-sum(Lf[:,:,t-1], dims=2))+ψ3*(P[:,t]-Pf[:,t])./(sum(Eh[:,:,t-1], dims=2)+sum(Eb[:,:,t-1], dims=2))
    append!(index, ψ1*(rL.*L[:,t-1]+sum(Pb[:,:,t], dims=2))+ψ3*sum(S[:,:,t], dims=1)./sum(Fh[:,:,t-1], dims=2))
    index = max.(zeros(O), index)
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
        lst = sample(1:(O+N), Weights(prob), Int64(x[l]))
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

function Eb_func(t)
    global Eb
    # 企業の株の保有額の見積もり
    Vs = sum(Eb[:,:,t-1], dims=1)+NLh[:,t]
    # ポートフォリオ配分先の確率の重みの共通部分を計算
    index = ψ1*(Pf[:,t] - I[:,t])+ψ2*(sum(Mf[:,:,t-1], dims=1)-sum(Lf[:,:,t-1], dims=2))+ψ3*(P[:,t]-Pf[:,t])./(sum(Eh[:,:,t-1], dims=2)+sum(Eb[:,:,t-1], dims=2))
    index = exp.(index./(sum(Eh[:,:,t])+sum(Eb[:,:,t])))
    index /= sum(index)
    # 銀行nのポートフォリオ計算
    for n=1:N
        # ポートフォリオ配分割合を決める
        prob = index*abs.(1+0.1*randn(N)) + Eb[:,n,t-1]/(sum(Eb[:,n,t-1]))
        prob /= sum(prob)
        # 収益率から、保有する株式の総額目標を計算
        E_volume = Vs[n]*(1 + λ2*sum(Pb[n,:,t])/sum(Eb[:,n,t-1]))
        # ポートフォリオ決定
        Eb[:,n,t] = E_volume*prob
    end
end

function ΔMh_and_Mh_func(t)
    global ΔMh, Mh
    change = rand(L) < ξ1
    prob = max.(zeros(N), NWb[:,t-1])
    prob /= sum(prob)
    for l=1:L
        ΔMh_sum = NLh[l,t]-sum(pe[:,t-1].*Δeh[:,l,t])+sum(ΔLh[l,:,t])
        last_n = 0
        for n=1:N
            if Mh[n,l,t-1] > 0.0
                last_n = n
                break
            end
        end
        if change
            n = sample(1:N, Weights(prob))
            Mh[n,l,t] = Mh[last_n,l,t-1] + ΔMh_sum
            ΔMh[n,l,t] = Mh[last_n,l,t-1] + ΔMh_sum
            ΔMh[last_n,l,t] = -Mh[last_n,l,t-1]
        else
            Mh[n,l,t] = Mh[n,l,t-1] + ΔMh_sum
            ΔMh[n,l,t] = ΔMh_sum
        end
    end
end

function ΔMf_and_Mf_func(t)
    global ΔMf, Mf
    change = rand(O) < ξ2
    prob = max.(zeros(N), NWb[:,t-1])
    prob /= sum(prob)
    for o=1:O
        ΔMf_sum = NLf[o,t]+sum(ΔLf[o,:,t])+sum(pe[:,t-1].*Δe[:,t])
        last_n = 0
        for n=1:N
            if Mf[n,o,t-1] > 0.0
                last_n = n
                break
            end
        end
        if change
            n = sample(1:N, Weights(prob))
            Mf[n,o,t] = Mf[last_n,o,t-1] + ΔMf_sum
            ΔMf[n,o,t] = Mf[last_n,o,t-1] + ΔMf_sum
            ΔMf[last_n,o,t] = -Mf[last_n,o,t-1]
        else
            Mf[n,o,t] = Mf[n,o,t-1] + ΔMf_sum
            ΔMf[n,o,t] = ΔMf_sum
        end
    end
end

function insolvency_disposition(t)
    global k, K, dis_k, dis_K
    global Mf, M, ΔMf, ΔM, dis_Mf, dis_ΔMf
    dis_lst, srv_lst = [], []
    for o=1:O
        if sum(Mf[:,o,t])<ϕ2*sum(Lf[o,:,t])
            push!(dis_lst, o)
        else
            push!(srv_lst, o)
        end
    end
    for o in dis_lst
        for n=1:N
            M[n,t] -= Mf[n,o,t]
            ΔM[n,t] -= ΔMf[n,o,t]
        end
    end
    
    k, K = k[srv_lst, :], K[srv_lst, :]
    Mf, ΔMf = Mf[srv_lst, :], ΔMf[srv_lst, :]
    dis_k, dis_K = cat(dis_k,k[dis_lst,:],dims=1), cat(dis_K,K[dis_lst,:],dims=1)
    dis_Mf, dis_ΔMf = cat(dis_Mf,Mf[:,dis_lst,:],dims=2), cat(dis_ΔMf,ΔMf[:,dis_lst,:],dims=2)

    return new_
end

function one_season(TIMERANGE)
    global p, pe, pf
    global Ti, Ta, Tv, Tc, w, W
    global g, G, c, C, i, I
    global k, K
    global P, Ph, Pf, Pb, S
    global A, u, DE, DF
    global NLh, NLf, NLb, NLg
    global Lh, Lf, L, ΔLh, ΔLf, ΔL
    global Mh, Mf, M, ΔMh, ΔMf, ΔM
    global eh, eb, e, Δeh, Δeb, Δe, Eh, Eb, E
    global fh, f, Δfh, Δf, Fh, F
    global H, ΔH
    global NWh, NWf, NWb, NWg, NW
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
        Eb_func(t)
        pe[:,t] = (sum(Eh[:,:,t], dims=2)+sum(Eb[:,:,t], dims=2))./e[:,t]
        for l=1:L
            eh[:,l,t] = Eh[:,l,t]./pe[:,t]
        end
        for n=1:N
            eb[:,n,t] = Eb[:,n,t]./pe[:,t]
        end
        Δeh[:,:,t] = eh[:,:,t] - eh[:,:,t-1]
        Δeb[:,:,t] = eb[:,:,t] - eb[:,:,t-1]
        ΔMh_and_Mh_func(t)
        ΔMf_and_Mf_func(t)
        M[:,t] = sum(Mf[:,:,t], dims=2)+sum(Mf[:,:,t], dims=2)
        ΔM[:,t] = M[:,t] - M[:,t-1]
        ΔH[:,t] = NLb[:,t]-[sum(pe[:,t].*Δeb[:,n,t]) for n=1:N]+ΔM[:,t]-ΔL[:,t]
        H[:,t] = H[:,t-1] + ΔH[:,t]
        NWh[:,t] = sum(Mh[:,:,t], dims=1)-sum(Lh[:,:,t], dims=2)+sum(Eh[:,:,t], dims=1)
        NWf[:,t] = K[:,t]+sum(Mf[:,:,t], dims=1)+sum(Lf[:,:,t], dims=2)-E[:,t]
        NWb[:,t] = -M[:,t]+L[:,t]+sum(Eb[:,:,t], dims=1)+H[:,t]
        NWg[t] = sum(H[:,t])
        DE[t] = sum(E[:,t])+sum(Eh[:,:,t])+sum(Eb[:,:,t])
        DF[t] = sum(Fh[:,:,t])-sum(F[:,t])
        # 倒産処理

        # 起業
    end
end