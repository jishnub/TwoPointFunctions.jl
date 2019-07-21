module TwoPointFunctions

using PointsOnASphere
using ForwardDiff

export cosχ

cosχ(x::Vector{<:Real}) = cos(x[1])cos(x[3]) + sin(x[1])sin(x[3])cos(x[2]-x[4])

∂cosχ(x) = ForwardDiff.gradient(cosχ,x)
∂²cosχ(x) = ForwardDiff.hessian(cosχ,x)

for fn in ["cosχ","∂cosχ","∂²cosχ"]
	@eval $(Symbol(fn))((θ₁,ϕ₁)::Tuple{<:Real,<:Real},(θ₂,ϕ₂)::Tuple{<:Real,<:Real}) = $(Symbol(fn))([θ₁,ϕ₁,θ₂,ϕ₂])
	@eval $(Symbol(fn))(θ₁::Real,ϕ₁::Real,θ₂::Real,ϕ₂::Real) = $(Symbol(fn))((θ₁,ϕ₁),(θ₂,ϕ₂))
	@eval $(Symbol(fn))((θ₁,ϕ₁)::Tuple{<:Real,<:Real},θ₂::Real,ϕ₂::Real) = $(Symbol(fn))((θ₁,ϕ₁),(θ₂,ϕ₂))
	@eval $(Symbol(fn))(θ₁::Real,ϕ₁::Real,(θ₂,ϕ₂)::Tuple{<:Real,<:Real}) = $(Symbol(fn))((θ₁,ϕ₁),(θ₂,ϕ₂))
	@eval $(Symbol(fn))((θ₁,ϕ₁)::Tuple{<:Real,<:Real},n::SphericalPoint) = $(Symbol(fn))((θ₁,ϕ₁),(n.θ,n.ϕ))
	@eval $(Symbol(fn))(θ₁::Real,ϕ₁::Real,n::SphericalPoint) = $(Symbol(fn))((θ₁,ϕ₁),(n.θ,n.ϕ))
	@eval $(Symbol(fn))(n::SphericalPoint,(θ₂,ϕ₂)::Tuple{<:Real,<:Real}) = $(Symbol(fn))((n.θ,n.ϕ),(θ₂,ϕ₂))
	@eval $(Symbol(fn))(n::SphericalPoint,θ₂::Real,ϕ₂::Real) = $(Symbol(fn))((n.θ,n.ϕ),(θ₂,ϕ₂))
	@eval $(Symbol(fn))(n₁::SphericalPoint,n₂::SphericalPoint) = $(Symbol(fn))((n₁.θ,n₁.ϕ),(n₂.θ,n₂.ϕ))
end

# use subscripted variables
subscripts = Dict(1=>'\u2081',2=>'\u2082') # \_1 and \_2
powers = Dict(2=>'\u00B2',3=>'\u00B3')
for n in 4:9
	powers[n] = '\u2070' + n
end

# First derivatives
for pt in 1:2, (coord_ind,coord) in enumerate(["θ","ϕ"])
	fname = Symbol("∂$(coord)$(subscripts[pt])cosχ")
	flat_ind = (pt-1)*2 + coord_ind 
	@eval ($fname)(x) = ∂cosχ(x)[$flat_ind]
	@eval ($fname)(x...) = ∂cosχ(x...)[$flat_ind]
	@eval export $fname
end

# Second derivatives
for pt1 in 1:2, (coord1_ind,coord1) in enumerate(["θ","ϕ"])
	for pt2 in 1:2, (coord2_ind,coord2) in enumerate(["θ","ϕ"])
		fname = Symbol("∂$(coord1)$(subscripts[pt1])∂$(coord2)$(subscripts[pt2])cosχ")
		coord1_flat_ind = (pt1-1)*2 + coord1_ind
		coord2_flat_ind = (pt2-1)*2 + coord2_ind
		@eval ($fname)(x) = ∂²cosχ(x)[$coord1_flat_ind,$coord2_flat_ind]
		@eval ($fname)(x...) = ∂²cosχ(x...)[$coord1_flat_ind,$coord2_flat_ind]
		@eval export $fname
	end
end

list_of_function_names = []

# Second θ-derivative wrt one point
for pt=1:2
	@eval $(Symbol("∂²θ$(subscripts[pt])cosχ"))((θ₁,ϕ₁),(θ₂,ϕ₂)) = $(Symbol("∂θ$(subscripts[pt])∂θ$(subscripts[pt])cosχ"))((θ₁,ϕ₁),(θ₂,ϕ₂))
	push!(list_of_function_names,"∂²θ$(subscripts[pt])cosχ")
end

# The following functions are not defined if the gradients are evaluated at the poles

∇ϕ₁cosχ((θ₁,ϕ₁)::Tuple{<:Real,<:Real},(θ₂,ϕ₂)::Tuple{<:Real,<:Real}) = -sin(θ₂)sin(ϕ₁-ϕ₂)
∇ϕ₂cosχ((θ₁,ϕ₁)::Tuple{<:Real,<:Real},(θ₂,ϕ₂)::Tuple{<:Real,<:Real}) = sin(θ₁)sin(ϕ₁-ϕ₂)
push!(list_of_function_names,"∇ϕ₁cosχ")
push!(list_of_function_names,"∇ϕ₂cosχ")

for pt in 1:2
	
	fname = "∂²ϕ$(subscripts[pt])cosχ"
	push!(list_of_function_names,fname)

	@eval $(Symbol(fname))(x) = $(Symbol("∂ϕ$(subscripts[pt])∂ϕ$(subscripts[pt])cosχ"))(x)
	@eval $(Symbol(fname))(x...) = $(Symbol("∂ϕ$(subscripts[pt])∂ϕ$(subscripts[pt])cosχ"))(x...)
end

∇ϕ₁∂ϕ₁cosχ((θ₁,ϕ₁)::Tuple{<:Real,<:Real},(θ₂,ϕ₂)::Tuple{<:Real,<:Real}) = -sin(θ₂)cos(ϕ₁-ϕ₂)
∇ϕ₁∂ϕ₂cosχ((θ₁,ϕ₁)::Tuple{<:Real,<:Real},(θ₂,ϕ₂)::Tuple{<:Real,<:Real}) = sin(θ₂)cos(ϕ₁-ϕ₂)
∂ϕ₁∇ϕ₂cosχ((θ₁,ϕ₁)::Tuple{<:Real,<:Real},(θ₂,ϕ₂)::Tuple{<:Real,<:Real}) = sin(θ₁)cos(ϕ₁-ϕ₂)
∇ϕ₂∂ϕ₂cosχ((θ₁,ϕ₁)::Tuple{<:Real,<:Real},(θ₂,ϕ₂)::Tuple{<:Real,<:Real}) = -sin(θ₂)cos(ϕ₁-ϕ₂)
∂ϕ₁∇ϕ₁cosχ((θ₁,ϕ₁)::Tuple{<:Real,<:Real},(θ₂,ϕ₂)::Tuple{<:Real,<:Real}) = ∇ϕ₁∂ϕ₁cosχ((θ₁,ϕ₁),(θ₂,ϕ₂))
∂ϕ₂∇ϕ₂cosχ((θ₁,ϕ₁)::Tuple{<:Real,<:Real},(θ₂,ϕ₂)::Tuple{<:Real,<:Real}) = ∇ϕ₂∂ϕ₂cosχ((θ₁,ϕ₁),(θ₂,ϕ₂))
push!(list_of_function_names,"∇ϕ₁∂ϕ₂cosχ")
push!(list_of_function_names,"∂ϕ₁∇ϕ₁cosχ")
push!(list_of_function_names,"∂ϕ₁∇ϕ₂cosχ")
push!(list_of_function_names,"∇ϕ₁∂ϕ₁cosχ")
push!(list_of_function_names,"∇ϕ₂∂ϕ₂cosχ")
push!(list_of_function_names,"∂ϕ₂∇ϕ₂cosχ")

# Higher-order derivatives wrt. ϕ
for pt in 1:2
	for order in 3:2:9
		fn = "∂$(powers[order])ϕ$(subscripts[pt])cosχ"
		f = @eval $(Symbol(fn))(x) = $(Symbol("∂ϕ$(subscripts[pt])cosχ"))(x) * real(im^($order-1))
		push!(list_of_function_names,f)
		@eval $(Symbol(fn))(x...) = $(Symbol("∂ϕ$(subscripts[pt])cosχ"))(x...) * real(im^($order-1))
	end
	for order in 4:2:8
		fn = "∂$(powers[order])ϕ$(subscripts[pt])cosχ"
		f = @eval $(Symbol(fn))(x) = $(Symbol("∂²ϕ$(subscripts[pt])cosχ"))(x) * real(im^($order-2))
		@eval $(Symbol(fn))(x...) = $(Symbol("∂²ϕ$(subscripts[pt])cosχ"))(x...) * real(im^($order-2))
		push!(list_of_function_names,f)
	end
end

# Add methods
for fn in list_of_function_names
	@eval $(Symbol(fn))(θ₁::Real,ϕ₁::Real,θ₂::Real,ϕ₂::Real) = $(Symbol(fn))((θ₁,ϕ₁),(θ₂,ϕ₂))
	@eval $(Symbol(fn))((θ₁,ϕ₁)::Tuple{<:Real,<:Real},θ₂::Real,ϕ₂::Real) = $(Symbol(fn))((θ₁,ϕ₁),(θ₂,ϕ₂))
	@eval $(Symbol(fn))(θ₁::Real,ϕ₁::Real,(θ₂,ϕ₂)::Tuple{<:Real,<:Real}) = $(Symbol(fn))((θ₁,ϕ₁),(θ₂,ϕ₂))
	@eval $(Symbol(fn))((θ₁,ϕ₁)::Tuple{<:Real,<:Real},n::SphericalPoint) = $(Symbol(fn))((θ₁,ϕ₁),(n.θ,n.ϕ))
	@eval $(Symbol(fn))(θ₁::Real,ϕ₁::Real,n::SphericalPoint) = $(Symbol(fn))((θ₁,ϕ₁),(n.θ,n.ϕ))
	@eval $(Symbol(fn))(n::SphericalPoint,(θ₂,ϕ₂)::Tuple{<:Real,<:Real}) = $(Symbol(fn))((n.θ,n.ϕ),(θ₂,ϕ₂))
	@eval $(Symbol(fn))(n::SphericalPoint,θ₂::Real,ϕ₂::Real) = $(Symbol(fn))((n.θ,n.ϕ),(θ₂,ϕ₂))
	@eval $(Symbol(fn))(n₁::SphericalPoint,n₂::SphericalPoint) = $(Symbol(fn))((n₁.θ,n₁.ϕ),(n₂.θ,n₂.ϕ))
end

for fn in list_of_function_names
	fname = Symbol(fn)
	@eval export $fname
end

end