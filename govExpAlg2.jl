using StatsPlots

N, T, G0 = 10^5, 1000, 1.0
β1, β2, β3 = 0.02, 0.02, 10.0
Gsum = [G0*(1+β1)^(t-1) for t=1:T]
x, z, Gm = ones(N), ones(N), zeros(N)
for t =1:T
    x .*= 1.0 .+ β2.*(-1.0 .+ β3*randn(N))
    x = max.(1.0, x)
    z = max.(0.0, x .- 2.0)
    Gm = z.*Gsum[t]/sum(z)
end

scatter(
    1:length(Gm[Gm.>0.0]), 
    sort(Gm[Gm.>0.0], rev=true), 
    yscale=:log10, 
    xscale=:log10, 
    label=nothing, 
    xlabel="Rank", 
    ylabel="Value"
)
histogram(Gm[Gm.>0.0], bins=50, yscale=:log10, label=nothing, xlabel="Value", ylabel="Intensity")
println(length(Gm[Gm.>0.0])/length(Gm))