using StatsPlots
using StatsBase
using Random

J, O, N, TIME,  = 20, 5, 2, 11
OBuffer = 50
Oin = 1

C, c, G, g, I, i = zeros(OBuffer, J, TIME), zeros(OBuffer, J, TIME), zeros(OBuffer, TIME), zeros(OBuffer, TIME), zeros(OBuffer, TIME), zeros(OBuffer, TIME)
w, W, Ti, Ta, Tv, Tc = zeros(J,TIME), zeros(J, OBuffer, TIME), zeros(J, TIME), zeros(J, TIME), zeros(OBuffer, TIME), zeros(OBuffer, TIME)
Ph, P, Pf, Pb = zeros(J, OBuffer, TIME), zeros(OBuffer, TIME), zeros(OBuffer, TIME), zeros(N, OBuffer, TIME)
S = zeros(J, N, TIME)
NLh, NLf, NLb, NLg = zeros(J, TIME), zeros(OBuffer, TIME), zeros(N, TIME), zeros(TIME)
p, pe, pf = zeros(OBuffer, TIME), zeros(OBuffer, TIME), zeros(N, TIME)
eh, e, eb, Δeh, Δe, Δeb = zeros(OBuffer, J, TIME), zeros(OBuffer, TIME), zeros(OBuffer, N, TIME), zeros(OBuffer, J, TIME), zeros(OBuffer, TIME), zeros(OBuffer, N, TIME)
Eh, E, Eb = zeros(OBuffer, J, TIME), zeros(OBuffer, TIME), zeros(OBuffer, N, TIME)
Fh, F, f, fh, Δfh = zeros(N, J, TIME), zeros(N, TIME), zeros(N, TIME), zeros(N, J, TIME), zeros(N, J, TIME)
Mh, Mf, M, ΔMh, ΔMf, ΔM = zeros(N, J, TIME), zeros(N, OBuffer, TIME), zeros(N, TIME), zeros(N, J, TIME), zeros(N, OBuffer, TIME), zeros(N, TIME)
Lh, Lf, L, ΔLh, ΔLf, ΔL = zeros(J, N, TIME), zeros(OBuffer, N, TIME), zeros(N, TIME), zeros(J, N, TIME), zeros(OBuffer, N, TIME), zeros(N, TIME)
H, ΔH = zeros(N, TIME), zeros(N, TIME)
K, k = zeros(OBuffer, TIME), zeros(OBuffer, TIME)
DE, DF = zeros(TIME), zeros(TIME)
NWh, NWf, NWb, NWg, NWg, NW = zeros(J, TIME), zeros(OBuffer, TIME), zeros(N, TIME), zeros(TIME), zeros(TIME), zeros(TIME)
u, v, A = zeros(OBuffer, TIME), zeros(OBuffer, TIME), zeros(OBuffer, TIME)
ΔZb, ΔZh = zeros(N, TIME), zeros(J, TIME)

EMP = zeros(Int64, J,TIME)
G_calc_item, G_potential = ones(O), ones(O)
os = [o for o=1:O]
last_os = [deepcopy(os)]

rL = 0.02
G0 = 1.0*O
uT = 0.8
c_base = 0.1

α1, α2, α3, α4 = 0.8, 0.02, 2.0, 0.01
β1, β2, β3, β4 = 0.02, 0.01, 0.2, 1.5 # 0.02, 0.02, 0.2, 2.0
γ1, γ2 = 1.0, 0.2
δ = 0.05
ϵ1, ϵ2 = 1.0, 1.0
λ1, λ2, λ3, λ4, λ5, λ6, λ7, λ8 = 0.3, 1.0, 0.5, 5.0, 0.5, 1.0, 0.1, 0.3
τ1, τ2, τ3, τ4 = 0.3, 0.02, 0.1, 0.2
μ1, μ2 = 0.0, 0.1
ν1, ν2 = 0.3, 0.5
θ1, θ2 = 0.2, 0.1
ϕ, ϕ2 = 1.0, 1.0
ψ1, ψ2, ψ3, ψ4 = 0.1, 0.01, 1.0, 1.0
ζ1, ζ2, ζ3 = 0.1, 0.02, 0.04 # 0.05,0.02,0.04
ξ1, ξ2 = 0.05, 0.05

# 関数定義

function G_func(t)
    global G, G_calc_item, G_potential
    Gsum = G0*(1+β1)^(t-1)
    G_calc_item .*= 1.0 .- β2 .+ β3*randn(O)
    G_calc_item = max.(1.0, G_calc_item)
    G_potential = max.(0.0, G_calc_item .- β4)
    if sum(G_potential)==0.0
        G_potential[sample(1:O)] = 1.0
    end
    G[os,t] = G_potential.*Gsum/sum(G_potential)
end

function c_func(t)
    global c
    # 消費先を変える場合どこに変えるか、の確率を決める
    r = zeros(O)
    for (q, o) in enumerate(os)
        if q==1
            r[1] = sum(c[1,:,t-1])+α4*(sum(c[os,:,t-1])+sum(g[os,t-1]))/J
        else
            r[q] = sum(c[o,:,t-1])+α4*(sum(c[os,:,t-1])+sum(g[os,t-1]))/J + r[q-1]
        end
    end
    r /= r[end]
    # 消費額を計算
    Cs = α1*(dropdims(sum(W[:,os,t-1], dims=2);dims=2)-Ta[:,t]-Ti[:,t]-rL*dropdims(sum(Lh[:,:,t-1], dims=2);dims=2)+dropdims(sum(Ph[:,os,t-1], dims=2);dims=2)+dropdims(sum(S[:,:,t-1], dims=2);dims=2)) 
        +α2*(dropdims(sum(Eh[os,:,t-1], dims=1);dims=1)+dropdims(sum(Mh[:,:,t-1], dims=1);dims=1)-dropdims(sum(Lh[:,:,t-1], dims=2);dims=2))
    Cs = max.(Cs, c_base*mean(p[os,t]))
    # 消費先の企業を選ぶ
    for j=1:J
        # 前期の消費先の特定
        last_o = findall(x -> x > 0.0, c[:,j,t-1])[1]
        # 消費先変更するかどうかの決定
        trans = rand() < α3*(p[last_o,t]-mean(p[os,t]))/p[last_o,t]
        if !(last_o in os)
            trans = true
        end
        # 消費先の選択
        if trans # 変更する場合
            new_o = os[1+searchsortedlast(r, rand())]
            c[new_o,j,t] = Cs[j]/p[new_o,t]
        else # 前回と同じ消費先を選ぶ場合
            c[last_o,j,t] = Cs[j]/p[last_o,t]
        end
    end
end

function w_and_W_func(t, flotation_info)    # wとWの関数の分離を検討
    global w, W, EMP, v
    flotations_js = [info[2] for info in flotation_info]
    flotations_os = [info[3] for info in flotation_info]
    # 求人数を作る
    offers = [Int64(round(max(0, (u[o,t-1]*k[o,t-1]-A[o,t-1]*sum(W[:,o,t-1].>0))/A[o,t-1]))) for o in os]
    # 応募確率を作る
    prob = zeros(O)
    for q in 1:O
        if q==1
            prob[1] = offers[1]
        else
            prob[q] = prob[q-1] + offers[q]
        end
    end
    if prob[end] > 0.0
        prob /= prob[end]
    end
    # 応募先リストの箱を用意
    appli = [[] for _ = 1:O]
    for j=1:J
        if EMP[j,t-1] > 0 # 前期就業していた人
            if !(EMP[j,t-1] in os)  # 勤め先が倒産した場合
                continue
            elseif rand() < ζ1 # 今期失業する場合
                W[j,EMP[j,t-1],t] = v[EMP[j,t-1],t-1]*w[j,t-1]
            else # 今期も同じ企業で就業する場合
                EMP[j,t] = EMP[j,t-1]
                w[j,t] = w[j,t-1]*(1+ζ2*abs(randn()))
                W[j,EMP[j,t-1],t] = v[EMP[j,t-1],t-1]*w[j,t-1]
            end
        elseif j in flotations_js # 起業メンバー
            for (x,j1) in enumerate(flotations_js)
                if j1==j
                    EMP[j,t-1] = flotations_os[x]
                    EMP[j,t] = EMP[j,t-1]
                    w[j,t] = w[j,t-1]*(1+ζ2*abs(randn()))
                    break
                end
            end
        else # 前期失業していた人
            if prob[end] == 0.0 # 募集がなければパス
                continue
            end
            # 応募先を割り振る
            q = 1+searchsortedlast(prob, rand())
            push!(appli[q], j)
        end
    end
    # マッチング
    for q=1:O
        if offers[q] > 0 & length(appli[q]) > 0
            much = collect(Set(sample(appli[q], offers[q])))
            for x=1:length(much)
                j = much[x]
                EMP[j,t] = os[q]
                w[j,t] = w[j,t-1]*(1+ζ2*abs(randn()))
            end
        end
    end
    # 従業員が0にならないように対策
    s = Set(EMP[:,t])
    UE = findall(x -> x == 0, EMP[:,t])
    if length(UE)==0
        return
    end
    for o in os
        if !(o in s)
            j = sample(UE)
            setdiff!(UE, [j])
            EMP[j,t] = o
            w[j,t] = w[j,t-1]*(1+ζ2*abs(randn()))
            if (EMP[j,t-1] > 0) & (EMP[j,t-1] in os)
                W[j,EMP[j,t-1],t] = v[EMP[j,t-1]]*w[j,t-1]
            end
        end
    end
    # 失業者の要求賃金率を下げる
    for j=1:J
        if EMP[j,t]==0
            w[j,t] = w[j,t-1]*(1-ζ3*abs(randn()))
        end
    end
end

function Lh_func(t)
    global Lh
    tmp = (dropdims(sum(C[:,:,t], dims=1);dims=1)+Ta[:,t]+Ti[:,t]+rL*dropdims(sum(Lh[:,:,t-1],dims=2);dims=2)) - dropdims(sum(Mh[:,:,t-1],dims=1);dims=1)
    Lhs = max.(0.0, ϵ1*NLh[:,t] + ϵ2*dropdims(sum(C[os,:,t], dims=1);dims=1), dropdims(sum(Lh[:,:,t-1],dims=2);dims=2)+1.01*tmp)
    for j=1:J
        for n=1:N
            if Mh[n,j,t-1] > 0.0
                Lh[j,n,t] = Lhs[j]
                break
            end
        end
    end
end

function ΔLf_and_Lf_func(t, FRL)
    global ΔLf, Lf
    ΔLfs = max.(-dropdims(sum(Lf[os,:,t-1], dims=2);dims=2), FRL)
    for (q,o) in enumerate(os)
        nLf, nMf = findfirst(x -> x > 0.0, Lf[o,:,t-1]), findfirst(x -> x > 0.0, Mf[:,o,t-1])
        if (nLf == nMf) | (nLf == nothing)
            Lf[o,nMf,t] = Lf[o,nMf,t-1] + ΔLfs[q]
            ΔLf[o,nMf,t] = ΔLfs[q]
        else
            Lf[o,nMf,t] = Lf[o,nLf,t-1] + ΔLfs[q]
            ΔLf[o,nMf,t] = Lf[o,nLf,t-1] + ΔLfs[q]
            ΔLf[o,nLf,t] = -Lf[o,nLf,t-1]
        end
    end
end

function household_portfolio_func(t, o_j_value_n_info)
    global Eh, Fh
    # 起業する企業の株式を除く、金融資産額の推定
    Vs = max.(0.0, 
            -dropdims(sum(Lh[:,:,t-1],dims=2);dims=2)
            +dropdims(sum(Mh[:,:,t-1], dims=1);dims=1)
            +dropdims(sum(Eh[os,:,t-1], dims=1);dims=1)
            +dropdims(sum(Fh[:,:,t-1], dims=1);dims=1)
            +NLh[:,t])
    for (o,j,value,n) in o_j_value_n_info
        Vs[j] -= value
    end
    # ポートフォリオ配分先の確率の重みの共通部分を計算
    index = ψ1*(Pf[os,t] - I[os,t])
            +ψ2*(dropdims(sum(Mf[:,os,t-1], dims=1);dims=1)-dropdims(sum(Lf[os,:,t-1], dims=2);dims=2))
            +ψ3*(P[os,t]-Pf[os,t])./(dropdims(sum(Eh[os,:,t-1], dims=2);dims=2)+dropdims(sum(Eb[os,:,t-1], dims=2);dims=2))
    append!(index, 
            ψ1*(rL.*L[:,t-1]+dropdims(sum(Pb[:,os,t], dims=2);dims=2))
            +ψ3*dropdims(sum(S[:,:,t], dims=1);dims=1)./dropdims(sum(Fh[:,:,t-1], dims=2);dims=2))
    index = max.(zeros(O+N), index)
    index /= sum(index)
    # ポートフォリオ配分先の数Int64(x[j])を決めるための準備。
    x = max.(0.0, 1.0*(Vs.-1.0*mean(sum(C[os,:,t], dims=1)))./mean(sum(C[os,:,t], dims=1)))
    # 家計jのポートフォリオ計算
    for j=1:J
        lastEhFhSum = sum(Eh[os,j,t-1])+sum(Fh[:,j,t-1])
        # ポートフォリオ配分先に選ぶ確率の重み付けを決める
        tmp = Eh[os,j,t-1]
        append!(tmp, Fh[:,j,t-1])
        prob = index
        if lastEhFhSum !== 0.0
            prob = index + tmp/(lastEhFhSum)
        end
        prob /= sum(prob)
        # ポートフォリオ配分先のリストを作る
        lst = sample(1:(O+N), Weights(prob), Int64(round(x[j])))
        if length(lst)==0
            continue
        end
        # 収益率から、資産に占める株式の割合の目標を計算し、保有額を決める
        EF_volume = Vs[j]*λ1
        if lastEhFhSum!==0.0
            EF_volume = Vs[j]*(λ1 + λ2*(sum(Ph[j,os,t])+sum(S[j,:,t]))/(lastEhFhSum))
        end
        # EF_volume/length(lst)を単位としてlstの企業または銀行の株の保有割合を決める
        for on in lst
            if on <= O
                Eh[os[on],j,t] += EF_volume/length(lst)
            else
                Fh[on-O,j,t] += EF_volume/length(lst)
            end
        end
        # 起業する企業の株式保有高を計算
        for (o,j,value, n) in o_j_value_n_info
            Eh[o,j,t] = value
        end
    end
end

function Eb_func(t)
    global Eb
    # ポートフォリオ配分先の確率の重みの共通部分を計算
    index = ψ1*(Pf[os,t] - I[os,t]) 
            +ψ2*(dropdims(sum(Mf[:,os,t-1], dims=1);dims=1)-dropdims(sum(Lf[os,:,t-1], dims=2);dims=2)) 
            +ψ3*(P[os,t]-Pf[os,t])./(dropdims(sum(Eh[os,:,t-1], dims=2);dims=2)+dropdims(sum(Eb[os,:,t-1], dims=2);dims=2)) 
            +ψ4*(dropdims(sum(Eh[os,:,t-1], dims=2);dims=2)+dropdims(sum(Eb[os,:,t-1], dims=2);dims=2))/(sum(Eh[os,:,t-1])+sum(Eb[os,:,t-1]))
    index = exp.(O*index./(sum(Eh[os,:,t-1])+sum(Eb[os,:,t-1])))
    index /= sum(index)
    # 銀行nのポートフォリオ計算
    for n=1:N
        # ポートフォリオ配分割合を決める
        prob = index.*abs.(1.0.+0.1*randn(O)) + Eb[os,n,t-1]/(sum(Eb[os,n,t-1]))
        prob /= sum(prob)
        # 収益率から、保有する株式の総額目標を決定
        EF_volume = (sum(Eb[os,n,t-1])+NLb[n,t])*(1.0-rL.+λ6*(sum(Pb[n,os,t])/sum(Eb[os,n,t-1]).-rL))
        # ポートフォリオ決定
        Eb[os,n,t] = EF_volume*prob
    end
end

function ΔMh_and_Mh_func(t, with_flo_os, o_j_value_n_info)
    global ΔMh, Mh
    change = rand(J) .< ξ1
    prob = max.(zeros(N), NWb[:,t-1])
    if prob[end]==0.0
        change .= false
    else
        prob /= sum(prob)
    end
    for j=1:J
        ΔMh_sum = NLh[j,t]-sum(pe[with_flo_os,t-1].*Δeh[with_flo_os,j,t])-sum(pf[:,t-1].*Δfh[:,j,t])+sum(ΔLh[j,:,t])-ΔZh[j,t]
        last_n = 0
        for n=1:N
            if Mh[n,j,t-1] > 0.0
                last_n = n
                break
            end
        end
        if change[j]
            new_n = sample(1:N, Weights(prob))
            Mh[new_n,j,t] = Mh[last_n,j,t-1] + ΔMh_sum
            ΔMh[new_n,j,t] = Mh[last_n,j,t-1] + ΔMh_sum
            ΔMh[last_n,j,t] -= Mh[last_n,j,t-1]
        else
            Mh[last_n,j,t] = Mh[last_n,j,t-1] + ΔMh_sum
            ΔMh[last_n,j,t] = ΔMh_sum
        end
    end
end

function ΔMf_and_Mf_func(t, o_j_value_n_info, os_adds)
    global ΔMf, Mf
    change = rand(O) .< ξ2
    prob = max.(zeros(N), NWb[:,t-1])
    if prob[end]==0.0
        change .= false
    else
        prob /= sum(prob)
    end
    for (x, o) in enumerate(os)
        ΔMf_sum = NLf[o,t]+sum(ΔLf[o,:,t])+sum(pe[o,t-1].*Δe[o,t])
        last_n = 0
        if o in os
            for n=1:N
                if Mf[n,o,t-1] > 0.0
                    last_n = n
                    break
                end
            end
        else    # 起業した企業
            for (new_o,j,value,n) in o_j_value_n_info
                if new_o==o
                    last_n = n
                    break
                end
            end
        end
        if change[x]
            new_n = sample(1:N, Weights(prob))
            Mf[new_n,o,t] = Mf[last_n,o,t-1] + ΔMf_sum
            ΔMf[new_n,o,t] = Mf[last_n,o,t-1] + ΔMf_sum
            ΔMf[last_n,o,t] -= Mf[last_n,o,t-1]
        else
            Mf[last_n,o,t] = Mf[last_n,o,t-1] + ΔMf_sum
            ΔMf[last_n,o,t] = ΔMf_sum
        end
    end
    for o in os_adds
        Mf[:,o,t] = ΔMf[:,o,t]
    end
end

function insolvency_disposition(t)
    global os, EMP, G_calc_item, G_potential, O
    global ΔMf, ΔLf, Δe, Δeh, Δeb
    global ΔZb, ΔZh
    q = 1
    while length(os) >= q
        o = os[q]
        if (sum(Mf[:,o,t])<ϕ2*sum(Lf[o,:,t])) & (P[o,t]-I[o,t]<0.0)
            setdiff!(os, o)
            deleteat!(G_calc_item, q)
            deleteat!(G_potential, q)
            EMP[EMP[:,t].==o, t] .= 0
            if t<TIME
                ΔMf[:,o,t+1] = -Mf[:,o,t]
                ΔLf[o,:,t+1] = -Lf[o,:,t]
                Δe[o,t+1] = -e[o,t]
                Δeh[o,:,t+1] = -eh[o,:,t]
                Δeb[o,:,t+1] = -eb[o,:,t]
                ΔZb[:,t+1] += ΔMf[:,o,t+1]
                ΔZb[:,t+1] -= ΔLf[o,:,t+1]
                ΔZb[:,t+1] .-= pe[o,t]*Δeb[o,:,t+1]
                ΔZh[:,t+1] .-= pe[o,t]*Δeh[o,:,t+1]
            end
            O -= 1
        else
            q += 1
        end
    end
end

function flotation(t) # 起業
    global G_calc_item, G_potential
    global ΔMh, ΔMf
    global Δe, Δeh
    global k, pe, p, A
    floations_info, o_j_value_n_info, os_adds = [], [], []
    next_o = 1
    for q=1:size(k)[1]
        if sum(k[q,:]) == 0.0
            next_o = q
            break
        end
    end
    prob = max.(0.0, NWh[:, t-1])
    prob /= sum(prob)
    floation_count = min(Oin, sum(EMP[:,t].==0.0))
    capitalists_js = sample(1:J, Weights(prob), floation_count, replace=false)
    workers_js = sample(findall(x -> x == 0, EMP[:,t]), floation_count, replace=false)
    for (x, o) in enumerate(next_o:next_o+floation_count-1)
        work_j = workers_js[x]
        cap_j = capitalists_js[x]

        # 企業情報リストに追加
        push!(floations_info, [cap_j, work_j, o])
        
        # 必要なストックの初期状態の定義
        cap_Mf = 0.0
        for n=1:N
            if Mh[n,cap_j,t-1] > 0.0
                # 企業の預金増加に伴う変化
                cap_Mf = 0.5*Mh[n,cap_j,t-1]
                #ΔMh[n,cap_j,t] -= cap_Mf
                ΔMf[n,o,t] = cap_Mf   # 資本家と企業が同じ銀行に口座を持つこととする

                # 投資した家計の情報を追加
                push!(o_j_value_n_info, (o, cap_j, cap_Mf, n))
                break
            end
        end
        # 資本家が株を持つことに伴う変化
        pe[o,t-1:t] .= 1.0
        Δe[o,t] = cap_Mf
        Δeh[o,cap_j,t] = cap_Mf
        # 企業が価格のつかない資本を持つことに伴う変化
        k[o,t] = 1.0
        A[o,t] = sample(A[os,t-1])

        # 希望価格
        p[o,t] = ν2*(1+ν1)*(w[work_j,t]+δ*k[o,t])/(uT*γ1*k[o,t])+(1-ν2)*(λ5*mean(p[o,t])+(1-λ5)*p[o,t])

        # 生存企業のリストにIDを追加
        push!(os_adds, o)
    end
    # 政府支出受注額決定のための配列を更新
    append!(G_calc_item, [1.0 for _=1:floation_count])
    append!(G_potential, [0.0 for _=1:floation_count])

    return floations_info, o_j_value_n_info, os_adds
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
    global last_os, O

    flotation_info = []
    for t=TIMERANGE
        println("##################### t=",t," ######################")

        p[os,t] = ν2*(1+ν1).*(dropdims(sum(W[:,os,t-1], dims=1);dims=1)+Tv[os,t-1]+Tc[os,t-1]+δ*k[os,t-1])./(uT*γ1*k[os,t-1])+(1-ν2)*(λ5*mean(p[os,t-1]).+(1-λ5)*p[os,t-1])
        w_and_W_func(t, flotation_info)
        Ti[:,t] = τ1*(dropdims(sum(W[:,os,t-1], dims=2);dims=2)+dropdims(sum(Ph[:,os,t-1], dims=2);dims=2)+dropdims(sum(S[:,:,t-1], dims=2);dims=2))
        Ta[:,t] = τ2*(dropdims(sum(Eh[os,:,t-1], dims=1);dims=1)+dropdims(sum(Fh[:,:,t-1], dims=1);dims=1)+dropdims(sum(Mh[:,:,t-1], dims=1);dims=1)-dropdims(sum(Lh[:,:,t-1], dims=2);dims=2))
        Tv[os,t] = τ3*(dropdims(sum(C[os,:,t-1], dims=2);dims=2)+I[os,t-1]+G[os,t-1])
        Tc[os,t] = max.(0.0, τ4*(dropdims(sum(C[os,:,t-1], dims=2);dims=2)+G[os,t-1]+I[os,t-1]-dropdims(sum(W[:,os,t-1], dims=1);dims=1)-Tv[os,t-1]))
        G_func(t)
        g[os,t] = G[os,t]./p[os,t]
        c_func(t)
        for o in os
            C[o,:,t] = p[o,t]*c[o,:,t]
        end
        A[os,t] = A[os,t-1].*(1+μ1.+μ2*i[os,t-1]./k[os,t-1])
        i[os,t] = max.(0.0, δ*k[os,t-1]+(u[os,t-1].-uT).*γ1.*k[os,t-1]+γ2*(dropdims(sum(Mf[:,os,t-1], dims=1);dims=1)-dropdims(sum(Lf[os,:,t-1], dims=2);dims=2))./p[os,t])
        I[os,t] = p[os,t].*i[os,t]
        u[os,t] = (dropdims(sum(c[os,:,t], dims=2);dims=2)+i[os,t]+g[os,t])./(γ1*k[os,t-1])
        k[os,t] = (1-δ).*k[os,t-1]+i[os,t]
        K[os,t] = p[os,t].*k[os,t]
        P[os,t] = dropdims(sum(C[os,:,t], dims=2);dims=2)+G[os,t]+I[os,t]-dropdims(sum(W[:,os,t], dims=1);dims=1)-Tv[os,t]-Tc[os,t]-rL*dropdims(sum(Lf[os,:,t-1], dims=2);dims=2)
        for o in os
            Ph[:,o,t] = max(0.0, θ1*(P[o,t]-I[o,t])+θ2*(sum(Mf[:,o,t-1])-sum(Lf[o,:,t-1]))).*eh[o,:,t-1]./e[o,t-1]
            Pb[:,o,t] = max(0.0, θ1*(P[o,t]-I[o,t])+θ2*(sum(Mf[:,o,t-1])-sum(Lf[o,:,t-1]))).*eb[o,:,t-1]./e[o,t-1]
        end
        Pf[os,t] = P[os,t] - dropdims(sum(Ph[:,os,t], dims=1);dims=1) - dropdims(sum(Pb[:,os,t], dims=1);dims=1)
        for n=1:N
            S[:,n,t] = max.(0.0, θ1*(rL*L[n,t-1]+sum(Pb[n,os,t]))+θ2*sum(Eb[os,n,t-1])).*fh[n,:,t-1]/f[n,t-1]
        end
        NLh[:,t] = -dropdims(sum(C[os,:,t], dims=1);dims=1) + dropdims(sum(W[:,os,t], dims=2);dims=2)-Ti[:,t]-Ta[:,t]-rL*dropdims(sum(Lh[:,:,t-1], dims=2);dims=2)+dropdims(sum(Ph[:,os,t], dims=2);dims=2)+dropdims(sum(S[:,:,t], dims=2);dims=2)
        NLf[os,t] = -I[os,t] + Pf[os,t]
        NLb[:,t] = rL*(dropdims(sum(Lh[:,:,t-1],dims=1);dims=1).+dropdims(sum(Lf[os,:,t-1],dims=1);dims=1)) + dropdims(sum(Pb[:,os,t], dims=2);dims=2) - dropdims(sum(S[:,:,t], dims=1);dims=1)
        NLg[t] = -sum(G[os,t])+sum(Ti[:,t])+sum(Ta[:,t])+sum(Tv[os,t])+sum(Tc[os,t])
        flotation_info, o_j_value_n_info, os_adds = flotation(t)   # 起業
        with_flo_os = cat(os, os_adds, dims=1)
        Lh_func(t)
        ΔLh[:,:,t] = Lh[:,:,t] - Lh[:,:,t-1]
        FR = I[os,t]+dropdims(sum(W[:,os,t], dims=1);dims=1)+Tv[os,t]+Tc[os,t]+rL*dropdims(sum(Lf[os,:,t-1], dims=2);dims=2)-ϕ*dropdims(sum(Mf[:,os,t-1], dims=1);dims=1)
        Δe[os,t] = min.(max.(-λ7*e[os,t-1], 1/p[os,t-1]*(1-λ3.-λ4*((P[os,t]-Pf[os,t])./(dropdims(sum(Eh[os,:,t-1], dims=2);dims=2)+dropdims(sum(Eb[os,:,t-1], dims=2);dims=2)).-rL)).*FR), λ8*e[os,t-1])
        ΔLf_and_Lf_func(t, FR.-pe[os,t-1].*Δe[os,t])
        L[:,t] = dropdims(sum(Lh[:,:,t], dims=1);dims=1) + dropdims(sum(Lf[:,:,t], dims=1);dims=1)
        ΔL[:,t] = dropdims(sum(ΔLh[:,:,t],dims=1);dims=1) + dropdims(sum(ΔLf[:,:,t],dims=1);dims=1)
        E[with_flo_os,t] = E[with_flo_os,t-1] + pe[with_flo_os,t-1].*Δe[with_flo_os,t]
        e[:,t] = e[:,t-1] + Δe[:,t]
        household_portfolio_func(t, o_j_value_n_info)
        f[:,t] = f[:,t-1]
        F[:,t] = F[:,t-1]
		pf[:,t] = dropdims(sum(Fh[:,:,t], dims=2);dims=2)./f[:,t]
        for j=1:J
            fh[:,j,t] = Fh[:,j,t]./pf[:,t]
        end
        Δfh[:,:,t] = fh[:,:,t] .- fh[:,:,t-1]
        Eb_func(t)
        pe[with_flo_os,t] = (dropdims(sum(Eh[with_flo_os,:,t], dims=2);dims=2)+dropdims(sum(Eb[with_flo_os,:,t], dims=2);dims=2))./e[with_flo_os,t]
        for j=1:J
            eh[with_flo_os,j,t] = Eh[with_flo_os,j,t]./pe[with_flo_os,t]
        end
        for n=1:N
            eb[os,n,t] = Eb[os,n,t]./pe[os,t]
        end
        Δeh[with_flo_os,:,t] = eh[with_flo_os,:,t] - eh[with_flo_os,:,t-1]
        Δeb[os,:,t] = eb[os,:,t] - eb[os,:,t-1]
        ΔMh_and_Mh_func(t, with_flo_os, o_j_value_n_info)
        ΔMf_and_Mf_func(t, o_j_value_n_info, os_adds)
        M[:,t] = dropdims(sum(Mh[:,:,t], dims=2);dims=2)+dropdims(sum(Mf[:,:,t], dims=2);dims=2)
        ΔM[:,t] = dropdims(sum(ΔMh[:,:,t],dims=2);dims=2) + dropdims(sum(ΔMf[:,:,t],dims=2);dims=2)
        ΔH[:,t] = NLb[:,t]-[sum(pe[with_flo_os,t-1].*Δeb[with_flo_os,n,t]) for n=1:N]+ΔM[:,t]-ΔL[:,t]-ΔZb[:,t]
        H[:,t] = H[:,t-1] + ΔH[:,t]
        NWh[:,t] = dropdims(sum(Mh[:,:,t], dims=1);dims=1)-dropdims(sum(Lh[:,:,t], dims=2);dims=2)+dropdims(sum(Eh[with_flo_os,:,t], dims=1);dims=1)+dropdims(sum(Fh[:,:,t], dims=1);dims=1)
        NWf[with_flo_os,t] = K[with_flo_os,t]+dropdims(sum(Mf[:,with_flo_os,t], dims=1);dims=1)-dropdims(sum(Lf[with_flo_os,:,t], dims=2);dims=2)-E[with_flo_os,t]
        NWb[:,t] = -M[:,t]+L[:,t]+dropdims(sum(Eb[os,:,t], dims=1);dims=1)-F[:,t]+H[:,t]
        NWg[t] = -sum(H[:,t])
        NW[t] = sum(NWh[:,t])+sum(NWf[with_flo_os,t])+sum(NWb[:,t])+NWg[t]
        DE[t] = -sum(E[with_flo_os,t])+sum(Eh[with_flo_os,:,t])+sum(Eb[os,:,t])
        DF[t] = sum(Fh[:,:,t])-sum(F[:,t])
        
        for o in os
            v[o,t] = (u[o,t]*k[o,t])./(A[o,t]*sum(EMP[:,t].==o))
        end
        
        println("transaction consistency : ",sum(NLh[:,t])+sum(NLf[:,t])+sum(NLb[:,t])+sum(NLg[t]))
        println("h flow consistency : ",sum(NLh[:,t])-sum([sum(pe[with_flo_os,t-1].*Δeh[with_flo_os,j,t]) for j=1:J])-sum([sum(pf[:,t-1].*Δfh[:,j,t]) for j=1:J])+sum(ΔLh[:,:,t])-sum(ΔMh[:,:,t])-sum(ΔZh[:,t]))
        println("f flow consistency : ",sum(NLf[:,t])+sum(pe[with_flo_os,t-1].*Δe[with_flo_os,t])+sum(ΔLf[with_flo_os,:,t])-sum(ΔMf[:,with_flo_os,t]))
        println("b flow consistency : ",sum(NLb[:,t])-sum([sum(pe[with_flo_os,t-1].*Δeb[with_flo_os,n,t]) for n=1:N])-sum(ΔL[:,t])+sum(ΔM[:,t])-sum(ΔH[:,t])-sum(ΔZb[:,t]))
        println("g flow consistency : ",NLg[t]+sum(ΔH[:,t]))
        println("stock consistency : ", sum(K[:,t])+DE[t]+DF[t]-NW[t])

        append!(os, os_adds)
        O += length(os_adds)
        push!(last_os, deepcopy(os))
        insolvency_disposition(t)   # 倒産処理
        println(os, last_os[end])
    end
end

# 初期値設定
function initialise()
    global EMP
    global Eh, Eb, E, pe, eh, eb, e
    global Fh, F, f, fh
    global Lh, Lf, L
    global K, k, p, H
    global NWh, NWf, NWb, NWg, NW

    for _=1:10
        G_func(1)
    end

    lstj, lsto, lstn = zeros(J), zeros(O), zeros(N)
    lstj[1], lsto[1], lstn[1] = 1.0, 1.0, 1.0
    #* 消費ネットワーク
    for j=1:J
        c[os,j,1] = shuffle(lsto)
    end
    #* 雇用ネットワーク
    EMP[:,1] = rand(0:O,J)
    w[:,1] .= 1.0
    for j=1:J
        if EMP[j,1]!==0
            W[j,EMP[j,1],1] = w[j,1]
        end
    end
    #* 企業株式保有ネットワーク
    for o in os
        Eh[o,:,1] = shuffle(lstj)
    end
    Eb[os,:,1] = rand(O,N)
    DIV = sum(Eh[os,:,1], dims=2)+sum(Eb[os,:,1], dims=2)
    for o in os
        Eh[o,:,1] /= DIV[o]
        Eb[o,:,1] /= DIV[o]
    end
    E[os,1] .= 1.0
    pe[os,1] .= 1.0
    eh[os,:,1], eb[os,:,1], e[os,1] = Eh[os,:,1], Eb[os,:,1], E[os,1]
    #* 銀行株式保有ネットワーク
    Fh[:,:,1] = rand(N,J)
    for j=1:J
        Fh[:,j,1] /= sum(Fh[:,j,1])
    end
    F[:,1] = sum(Fh[:,:,1],dims=2)
    fh[:,:,1] .= 1.0/J
    f[:,1] .= 1.0
    pf[:,1] = F[:,1]./f[:,1]
    #* 預金貸借ネットワーク
    for j=1:J
        Mh[:,j,1] = shuffle(10.0.*lstn)
    end
    for o in os
        Mf[:,o,1] = shuffle(lstn)
    end
    M[:,1] = sum(Mh[:,:,1],dims=2)+sum(Mf[:,:,1],dims=2)
    #* 初期資本と初期価格
    k[os,1] .= J/O
    for o in os
        p[o,1] = (1+ν1)*(sum(W[:,o,1])+δ*k[o,1])/(uT*γ1*k[o,1])
    end
    K[os,1] = p[os,1].*k[os,1]
    A[os,1] .= 1.0
    #* 現金ネットワーク
    H[:,1] .= dropdims(sum(Mh[:,:,1],dims=2);dims=2) + dropdims(sum(Mf[:,:,1],dims=2);dims=2)
                -dropdims(sum(Lh[:,:,1],dims=1);dims=1) - dropdims(sum(Lf[:,:,1],dims=1);dims=1)
    #* ストックの整合性を保証
    NWh[:,1] = dropdims(sum(Mh[:,:,1],dims=1); dims=1)-dropdims(sum(Lh[:,:,1],dims=2);dims=2)+dropdims(sum(Eh[:,:,1],dims=1);dims=1)+dropdims(sum(Fh[:,:,1],dims=1);dims=1)
    NWf[os,1] = K[os,1]+dropdims(sum(Mf[:,os,1],dims=1);dims=1)-dropdims(sum(Lf[os,:,1],dims=2),dims=2)-E[os,1]
    NWb[:,1] = -M[:,1]+L[:,1]+dropdims(sum(Eb[os,:,1],dims=1);dims=1)-F[:,1]+H[:,1]
    NWg[1] = -sum(H[:,1])
    NW[1] = sum(K[:,1])
end

initialise()


# シミュレーション
one_season(2:TIME)

# プロット
plot(dropdims(sum(NLh[:,:],dims=1);dims=1), label="NLh")
plot!(dropdims(sum(NLf[:,:],dims=1);dims=1), label="NLf")
plot!(dropdims(sum(NLb[:,:],dims=1);dims=1), label="NLb")
plot!(NLg, label="NLg")
savefig("AB_model/figs/NL_sum.png")

plot(dropdims(sum(NWh[:,:],dims=1);dims=1), label="NWh")
plot!(dropdims(sum(NWf[:,:],dims=1);dims=1), label="NWf")
plot!(dropdims(sum(NWb[:,:],dims=1);dims=1), label="NWb")
plot!(NWg, label="NWg")
plot!(NW, label="NW")
savefig("AB_model/figs/NW_sum.png")

println("###########",NLh[1,TIME]-sum(pe[last_os[TIME],TIME-1].*Δeh[last_os[TIME],1,TIME])-sum(pf[:,TIME-1].*Δfh[:,1,TIME])+sum(ΔLh[1,:,TIME])-sum(ΔMh[:,1,TIME])-ΔZh[1,TIME])
println("###########",NLf[last_os[TIME][1],TIME]+pe[last_os[TIME][1],TIME-1]*Δe[last_os[TIME][1],TIME]+sum(ΔLf[last_os[TIME][1],:,TIME])-sum(ΔMf[:,last_os[TIME][1],TIME]))
println("###########",NLb[1,TIME]-sum(pe[last_os[TIME],TIME-1].*Δeb[last_os[TIME],1,TIME])-ΔL[1,TIME]+ΔM[1,TIME]-ΔH[1,TIME]-ΔZb[1,TIME])
println("###########",NLg[TIME]+sum(ΔH[:,TIME]))
println("%%%%%%%%%%%",sum(NLh[:,TIME])-sum([sum(pe[last_os[TIME],TIME-1].*Δeh[last_os[TIME],j,TIME]) for j=1:J])+sum(ΔLh[:,:,TIME])-sum(ΔMh[:,:,TIME])-sum(ΔZh[:,TIME]))
println("%%%%%%%%%%%",sum(NLf[last_os[TIME],TIME])+sum(pe[last_os[TIME],TIME-1].*Δe[last_os[TIME],TIME])+sum(ΔLf[last_os[TIME],:,TIME])-sum(ΔMf[:,last_os[TIME],TIME]))
println("%%%%%%%%%%%",sum(NLb[:,TIME])-sum([sum(pe[last_os[TIME],TIME-1].*Δeb[last_os[TIME],n,TIME]) for n=1:N])-sum(ΔL[:,TIME])+sum(ΔM[:,TIME])-sum(ΔH[:,TIME])-sum(ΔZb[:,TIME]))
println("%%%%%%%%%%%",NLg[TIME]+sum(ΔH[:,TIME]))

plot(2:TIME, dropdims(sum(NLh[:,2:TIME],dims=1);dims=1).+dropdims(sum(NLf[:,2:TIME],dims=1);dims=1).+dropdims(sum(NLb[:,2:TIME],dims=1);dims=1).+NLg[2:TIME])
savefig("AB_model/figs/transaction_consistency.png")

plot(2:TIME,[sum(NLh[:,t])-sum([sum(pe[last_os[t],t-1].*Δeh[last_os[t],j,t]) for j=1:J])-sum([sum(pf[:,t-1].*Δfh[:,j,t]) for j=1:J])-sum(ΔMh[:,:,t])+sum(ΔLh[:,:,t])-sum(ΔZh[:,t]) for t=2:TIME],label="h")
plot!(2:TIME,[sum(NLf[:,t])+sum(pe[last_os[t],t-1].*Δe[last_os[t],t])+sum(ΔLf[last_os[t],:,t])-sum(ΔMf[:,last_os[t],t]) for t=2:TIME],label="f")
plot!(2:TIME,[sum(NLb[:,t])-sum([sum(pe[last_os[t],t-1].*Δeb[last_os[t],n,t]) for n=1:N])-sum(ΔL[:,t])+sum(ΔM[:,t])-sum(ΔH[:,t])-sum(ΔZb[:,t]) for t=2:TIME],label="b")
plot!(2:TIME,[NLg[t]+sum(ΔH[:,t]) for t=2:TIME],label="g")
savefig("AB_model/figs/flow_consistency.png")

plot(2:TIME, [dropdims(sum(p[last_os[t],t],dims=1);dims=1)./length(last_os[t]) for t=2:TIME], label="p_average")
savefig("AB_model/figs/p_average.png")
plot(2:TIME, [dropdims(sum(pe[last_os[t],t],dims=1);dims=1)./length(last_os[t]) for t=2:TIME], label="pe_average")
savefig("AB_model/figs/pe_average.png")
plot(2:TIME, dropdims(sum(pf[:,2:TIME],dims=1);dims=1)./N, label="pf_average")
savefig("AB_model/figs/pf_average.png")

plot(dropdims(sum(K[:,:],dims=1);dims=1).+DE.+DF.-NW, label="all")
plot!([sum(Mh[:,:,t])-sum(Lh[:,:,t])+sum(Eh[:,:,t])+sum(Fh[:,:,t])-sum(NWh[:,t]) for t=1:TIME], label="h")
plot!([sum(K[:,t])+sum(Mf[:,:,t])-sum(Lf[:,:,t])-sum(E[:,t])-sum(NWf[:,t]) for t=1:TIME], label="f")
plot!([-sum(M[:,t])+sum(L[:,t])+sum(Eb[:,:,t])-sum(F[:,t])+sum(H[:,t])-sum(NWb[:,t]) for t=1:TIME], label="b")
plot!([-sum(H[:,t])-NWg[t] for t=1:TIME], label="g")
savefig("AB_model/figs/stock_consistency.png")

plot(2:TIME, [sum(L[:,t])-(sum(L[:,t-1])+sum(ΔL[:,t])) for t=2:TIME], label="L")
plot!(2:TIME, [sum(M[:,t])-(sum(M[:,t-1])+sum(ΔM[:,t])) for t=2:TIME], label="M")
savefig("AB_model/figs/TMP.png")

plot(K[1,:],label="1")
for o=2:os[end]
    plot!(K[o,:], label=string(o))
end
savefig("AB_model/figs/Ks.png")

plot(p[1,:],label="1")
for o=2:os[end]
    plot!(p[o,:], label=string(o))
end
savefig("AB_model/figs/ps.png")
