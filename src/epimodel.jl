
f(u, p, t) = begin
  ϕ₁ = p[1]*u[1]*u[2]/sum(u)
  ϕ₂ = p[2]*u[2]
  return (-ϕ₁, ϕ₁-ϕ₂, ϕ₂)
end


function euler(f, u₀, p, tspan, Δt)
  tsteps = tspan[1]:Δt:tspan[2]
  U = zeros(eltype(u₀), length(u₀), length(tsteps))
  U[:, 1] = u₀
  for (i,t) in enumerate(tsteps[2:end])
    u = U[:,i]
    uₙ = u .+ Δt .* f(u, p, t)
    U[:,i+1] = uₙ
  end
  return tsteps, U
end


t₀ = 0.0
t₁ = 1.0
Δt = 0.01
u₀ = [100.0, 1, 0]
p = [25.0, 3.5]
tsteps, soln = euler(f, u₀, p, (t₀,t₁), Δt)


using Plots
using PlotThemes
theme(:default, w=5)


plot(tsteps, soln')


fsirs(u, p, t) = begin
  ϕ₁ = p[1]*u[1]*u[2]/sum(u)
  ϕ₂ = p[2]*u[2]
  ϕ₃ = p[3]*u[3]
  return (ϕ₃-ϕ₁, ϕ₁-ϕ₂, ϕ₂-ϕ₃)
end

psirs = [25.0, 3.5, 2.0]
tsteps, soln = euler(fsirs, u₀, psirs, (t₀,t₁), Δt)


plot(tsteps, soln')


tsteps, soln = euler(fsirs, u₀, psirs, (t₀,2t₁), Δt)
plot(tsteps, soln')


fsirs(u₀, psirs, 0)
tsteps, soln = euler(fsirs, u₀, psirs, (t₀,2t₁), Δt)
uₑ = soln[:, end]
plot(tsteps, soln')


fsirs(uₑ, psirs, 2t₁)
tsteps, soln = euler(fsirs, u₀, psirs, (t₀,3t₁), Δt)
uₑ = soln[:, end]
plot(tsteps, soln')


fsirs(uₑ, psirs, 3t₁)
derivs = zeros(Float64, length(tsteps), length(u₀))
map(1:length(tsteps)) do i
  derivs[i,:] .= fsirs(soln[:, i], psirs, tsteps[i])
end

derivs
plot(tsteps, derivs)


using NLsolve

# We have to wrap our vector field function to satisfy the API expected by NLsolve
fsirs!(du, u) = begin
  du .= fsirs(u, psirs, 0)
  return du
end


result = nlsolve(fsirs!, u₀/sum(u₀), autodiff = :forward, method=:newton)


result.zero


tsteps, soln = euler(fsirs, u₀/sum(u₀), psirs, (t₀,t₁/3), Δt)
uₑ = soln[:, end]
result = nlsolve(fsirs!, uₑ/sum(uₑ), autodiff = :forward, method=:newton)
result.zero


tsteps, soln = euler(fsirs, u₀/sum(u₀), psirs, (t₀,3t₁), Δt)
uₑ = soln[:, end]
plot(tsteps, soln')


norm(x) = sqrt(x'x)
uₑ/norm(uₑ) - result.zero/norm(result.zero)


newtonstep(f, J, xₖ) = xₖ + (J(xₖ)) \ -collect(f(xₖ))

using ForwardDiff
using LinearAlgebra
fsirs₁(u) = collect(fsirs(u, psirs, 0))
Jsirs(x) = ForwardDiff.jacobian(fsirs₁, x)


newtonstep(fsirs₁, Jsirs, u₀/sum(u₀))


function newtoniteration(f, J, x₀, nsteps)
  x = x₀
  for i in 1:nsteps
    if norm(f(x)) < 1e-8 # we did it!
      return x
    end
    x .= newtonstep(f,J, x) # keep going!
  end
  return x
end
newtoniteration(fsirs₁, Jsirs, soln[:, 250], 10)


try 
  newtoniteration(fsirs₁, Jsirs, soln[:, 249], 10)
catch SingularException
  println("Bad News Bears")
end


Jₑ = Jsirs(soln[:,end])


Λ  = eigvals(Jₑ)


abs.(Λ)


plot(tsteps[1:end÷2], real.(exp.(-Λ'.*(tsteps[1:end÷2]))))

