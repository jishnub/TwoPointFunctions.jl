module TwoPointFunctions

using ForwardDiff
using StaticArrays

export cosχ
export ∇ϕ₁cosχ
export ∇ϕ₂cosχ
export ∂²θ₁cosχ
export ∂²θ₂cosχ
export ∂²ϕ₁cosχ
export ∂²ϕ₂cosχ
export ∇₁cosχ
export ∇₂cosχ

_cosχ2(x) = sin(x[1])sin(x[3])cos(x[2] - x[4])
_cosχ2((θ₁, ϕ₁), (θ₂, ϕ₂)) = sin(θ₁)sin(θ₂)cos(ϕ₁ - ϕ₂)
cosχ(x::AbstractVector{<:Real}) = cos(x[1])cos(x[3]) + _cosχ2(x)

∂cosχ(x) = ForwardDiff.gradient(cosχ, x)
∂²cosχ(x) = ForwardDiff.hessian(cosχ, x)

"""
	cosχ((θ₁,ϕ₁), (θ₂,ϕ₂))

Cosine of the angle between the unit vectors denoted by the spherical coordinates 
``(\\theta_1,\\phi_1)`` and ``(\\theta_2,\\phi_2)``. This is equivalent to the cosine of the 
angular distance along a great circle between these two points on a sphere.

Numerically this is given by
```math
\\cos\\chi((\\theta_1,\\phi_1), (\\theta_2,\\phi_2)) = \\cos(\\theta_1)\\cos(\\theta_2) + 
\\sin(\\theta_1)\\sin(\\theta_2)\\cos(\\phi_1 - \\phi_2)
```

The unicode character at the end of the function's name is `\\chi`, and not `x`.

# Examples
```jldoctest
julia> cosχ((0,0), (0,0))
1.0

julia> cosχ((0,0), (pi,0))
-1.0
```
"""
function cosχ end

"""
	∂θ₁cosχ((θ₁,ϕ₁), (θ₂,ϕ₂))

Partial derivative of [`cosχ`](@ref) with respect to `θ₁`.

# Examples
```jldoctest
julia> ∂θ₁cosχ((0,0), (0,0))
0.0

julia> ∂θ₁cosχ((0,0), (pi/2,0))
1.0
```

See also: [`∂θ₂cosχ`](@ref)
"""
function ∂θ₁cosχ end

"""
	∂θ₂cosχ((θ₁,ϕ₁), (θ₂,ϕ₂))

Partial derivative of [`cosχ`](@ref) with respect to `θ₂`.

# Examples
```jldoctest
julia> ∂θ₂cosχ((0,0), (0,0))
0.0

julia> ∂θ₂cosχ((0,0), (pi/2,pi/2))
-1.0
```

See also: [`∂θ₁cosχ`](@ref) 
"""
function ∂θ₂cosχ end

"""
	∂ϕ₁cosχ((θ₁,ϕ₁), (θ₂,ϕ₂))

Partial derivative of [`cosχ`](@ref) with respect to `ϕ₁`.

# Examples
```jldoctest
julia> ∂ϕ₁cosχ((pi/2, 0), (pi/2, pi/3)) ≈ sin(pi/3)
true

julia> ∂ϕ₁cosχ((pi/2, 0), (pi/2, pi/2))
1.0
```

See also: [`∂ϕ₂cosχ`](@ref)
"""
function ∂ϕ₁cosχ end

"""
	∂ϕ₂cosχ((θ₁,ϕ₁), (θ₂,ϕ₂))

Partial derivative of [`cosχ`](@ref) with respect to `ϕ₂`.

# Examples
```jldoctest
julia> ∂ϕ₂cosχ((pi/2, 0), (pi/2, pi/3)) ≈ -sin(pi/3)
true

julia> ∂ϕ₂cosχ((pi/2, 0), (pi/2, pi/2))
-1.0
```

See also: [`∂ϕ₁cosχ`](@ref)
"""
function ∂ϕ₂cosχ end

for fn in ["cosχ","∂cosχ","∂²cosχ"]
	@eval $(Symbol(fn))((θ₁,ϕ₁)::Tuple{<:Real,<:Real}, (θ₂,ϕ₂)::Tuple{<:Real,<:Real}) = $(Symbol(fn))(@SVector [θ₁,ϕ₁,θ₂,ϕ₂])
end

# use subscripted variables
subscripts = Dict(1=>'\u2081',2=>'\u2082') # \_1 and \_2
powers = Dict(2=>'\u00B2',3=>'\u00B3')
for n in 4:9
	powers[n] = '\u2070' + n
end

# First derivatives
for pt in 1:2, (coord_ind,coord) in enumerate(("θ","ϕ"))
	fname = Symbol("∂$(coord)$(subscripts[pt])cosχ")
	flat_ind = (pt-1)*2 + coord_ind 
	@eval ($fname)(n1, n2) = ∂cosχ(n1, n2)[$flat_ind]
	@eval export $fname
end

# Second derivatives
for pt1 in 1:2, (coord1_ind,coord1) in enumerate(("θ","ϕ"))
	for pt2 in 1:2, (coord2_ind,coord2) in enumerate(("θ","ϕ"))
		fname = Symbol("∂$(coord1)$(subscripts[pt1])∂$(coord2)$(subscripts[pt2])cosχ")
		coord1_flat_ind = (pt1-1)*2 + coord1_ind
		coord2_flat_ind = (pt2-1)*2 + coord2_ind
		@eval ($fname)(n1, n2) = ∂²cosχ(n1, n2)[$coord1_flat_ind,$coord2_flat_ind]
	end
end

"""
	∂²θ₁cosχ((θ₁,ϕ₁), (θ₂,ϕ₂))

Second derivative of [`cosχ`](@ref) with respect to `θ₁`.

# Examples
```jldoctest
julia> ∂²θ₁cosχ((0,0), (0,0))
-1.0

julia> ∂²θ₁cosχ((0,0), (pi,0))
1.0
```
"""
∂²θ₁cosχ(n1, n2) = -cosχ(n1, n2)

"""
	∂²θ₂cosχ((θ₁,ϕ₁), (θ₂,ϕ₂))

Second derivative of [`cosχ`](@ref) with respect to `θ₂`.

# Examples
```jldoctest
julia> ∂²θ₂cosχ((0,0), (0,0))
-1.0

julia> ∂²θ₂cosχ((0,0), (pi,0))
1.0
```
"""
∂²θ₂cosχ(n1, n2) = -cosχ(n1, n2)

"""
	TwoPointFunctions.∂θ₁∂θ₂cosχ((θ₁,ϕ₁), (θ₂,ϕ₂))

Mixed derivative of [`cosχ`](@ref) with respect to `θ₁` and `θ₂`. 

# Examples
```jldoctest
julia> TwoPointFunctions.∂θ₁∂θ₂cosχ((0,0), (0,0))
1.0

julia> TwoPointFunctions.∂θ₁∂θ₂cosχ((0,0), (pi,0))
-1.0
```
"""
function ∂θ₁∂θ₂cosχ end

# The following functions are not defined if the gradients are evaluated at the poles

"""
	∇ϕ₁cosχ((θ₁,ϕ₁), (θ₂,ϕ₂))

The ``{\\hat{\\phi}}_1`` component of ``\\nabla_1 \\cos \\chi`` in the spherical polar basis.
Numerically this is equivalent to 

```math
\\nabla \\phi_1 \\cos\\chi ((\\theta_1, \\phi_1), (\\theta_2, \\phi_2)) = 
\\frac{1}{\\sin(\\theta_1)} \\partial \\phi_1 \\cos\\chi ((\\theta_1, \\phi_1), (\\theta_2, \\phi_2))
```

except `∇ϕ₁cosχ` is defined to avoid coordinate singularities at the poles.

# Examples
```jldoctest
julia> ∇ϕ₁cosχ((0, 0), (0, 0))
-0.0

julia> ∇ϕ₁cosχ((pi/2, 0), (pi/2, pi/2))
1.0
```

See also: [`∂ϕ₁cosχ`](@ref), [`∇ϕ₂cosχ`](@ref)
"""
∇ϕ₁cosχ((θ₁,ϕ₁)::Tuple{<:Real,<:Real}, (θ₂,ϕ₂)::Tuple{<:Real,<:Real}) = -sin(θ₂)sin(ϕ₁-ϕ₂)

"""
	∇ϕ₂cosχ((θ₁,ϕ₁), (θ₂,ϕ₂))

The ``{\\hat{\\phi}}_2`` component of ``\\nabla_2 \\cos \\chi`` in the spherical polar basis.
Numerically this is equivalent to 

```math
\\nabla \\phi_2 \\cos\\chi ((\\theta_1, \\phi_1), (\\theta_2, \\phi_2)) = 
\\frac{1}{\\sin(\\theta_2)} \\partial \\phi_2 \\cos\\chi ((\\theta_1, \\phi_1), (\\theta_2, \\phi_2))
```

except `∇ϕ₂cosχ` is defined to avoid coordinate singularities at the poles.

# Examples
```jldoctest
julia> ∇ϕ₂cosχ((0, 0), (0, 0))
0.0

julia> ∇ϕ₂cosχ((pi/2, 0), (pi/2, pi/2))
-1.0
```

See also: [`∂ϕ₂cosχ`](@ref), [`∇ϕ₁cosχ`](@ref)
"""
∇ϕ₂cosχ((θ₁,ϕ₁)::Tuple{<:Real,<:Real}, (θ₂,ϕ₂)::Tuple{<:Real,<:Real}) = sin(θ₁)sin(ϕ₁-ϕ₂)

"""
	∂²ϕ₁cosχ((θ₁,ϕ₁), (θ₂,ϕ₂))

Second derivative of [`cosχ`](@ref) with respect to `ϕ₁`.

# Examples
```jldoctest
julia> ∂²ϕ₁cosχ((0,0), (0,0))
-0.0

julia> ∂²ϕ₁cosχ((pi/2,0), (pi/2, 0))
-1.0
```
"""
function ∂²ϕ₁cosχ end

"""
	∂²ϕ₂cosχ((θ₁,ϕ₁), (θ₂,ϕ₂))

Second derivative of [`cosχ`](@ref) with respect to `ϕ₂`.

# Examples
```jldoctest
julia> ∂²ϕ₂cosχ((0,0), (0,0))
-0.0

julia> ∂²ϕ₂cosχ((pi/2,0), (pi/2, 0))
-1.0
```
"""
function ∂²ϕ₂cosχ end

for pt in 1:2
	∂²ϕ = "∂²ϕ$(subscripts[pt])cosχ"
	∂ϕ∂ϕ = "∂ϕ$(subscripts[pt])∂ϕ$(subscripts[pt])cosχ"
	@eval $(Symbol(∂²ϕ))(n1, n2) = -_cosχ2(n1, n2)
end

∇ϕ₁∂ϕ₁cosχ((θ₁,ϕ₁)::Tuple{<:Real,<:Real}, (θ₂,ϕ₂)::Tuple{<:Real,<:Real}) = -sin(θ₂)cos(ϕ₁-ϕ₂)
∇ϕ₁∂ϕ₂cosχ((θ₁,ϕ₁)::Tuple{<:Real,<:Real}, (θ₂,ϕ₂)::Tuple{<:Real,<:Real}) = sin(θ₂)cos(ϕ₁-ϕ₂)
∂ϕ₁∇ϕ₂cosχ((θ₁,ϕ₁)::Tuple{<:Real,<:Real}, (θ₂,ϕ₂)::Tuple{<:Real,<:Real}) = sin(θ₁)cos(ϕ₁-ϕ₂)
∇ϕ₂∂ϕ₂cosχ((θ₁,ϕ₁)::Tuple{<:Real,<:Real}, (θ₂,ϕ₂)::Tuple{<:Real,<:Real}) = -sin(θ₂)cos(ϕ₁-ϕ₂)
∂ϕ₁∇ϕ₁cosχ((θ₁,ϕ₁)::Tuple{<:Real,<:Real}, (θ₂,ϕ₂)::Tuple{<:Real,<:Real}) = ∇ϕ₁∂ϕ₁cosχ((θ₁,ϕ₁),(θ₂,ϕ₂))
∂ϕ₂∇ϕ₂cosχ((θ₁,ϕ₁)::Tuple{<:Real,<:Real}, (θ₂,ϕ₂)::Tuple{<:Real,<:Real}) = ∇ϕ₂∂ϕ₂cosχ((θ₁,ϕ₁),(θ₂,ϕ₂))

# Higher-order derivatives wrt. ϕ
for pt in 1:2
	for order in 3:2:9
		fn = "∂$(powers[order])ϕ$(subscripts[pt])cosχ"
		@eval $(Symbol(fn))(n1, n2) = $(Symbol("∂ϕ$(subscripts[pt])cosχ"))(n1, n2) * real(im^($order-1))
	end
	for order in 4:2:8
		fn = "∂$(powers[order])ϕ$(subscripts[pt])cosχ"
		@eval $(Symbol(fn))(n1, n2) = $(Symbol("∂²ϕ$(subscripts[pt])cosχ"))(n1, n2) * real(im^($order-2))
	end
end

"""
	∇₁cosχ((θ₁,ϕ₁), (θ₂,ϕ₂))

Return 
``(\\hat{\\theta}_1 \\cdot \\nabla_1 \\cos \\chi,\\,
\\hat{\\phi}_1 \\cdot \\nabla_1 \\cos \\chi)``, where ``{\\nabla}_1`` represents 
the gradient on the unit sphere with respect to ``(\\theta_1, \\phi_1)``.

# Examples
```jldoctest
julia> ∇₁cosχ((0,0), (pi/2, pi/3)) == (cos(pi/3), sin(pi/3))
true
```
"""
function ∇₁cosχ end

"""
	∇₂cosχ((θ₁,ϕ₁), (θ₂,ϕ₂))

Return 
``(\\hat{\\theta}_2 \\cdot \\nabla_2 \\cos \\chi,\\, 
\\hat{\\phi}_2 \\cdot \\nabla_2 \\cos \\chi)``, where ``{\\nabla}_2`` represents 
the gradient on the unit sphere with respect to ``(\\theta_2, \\phi_2)``.

# Examples
```jldoctest
julia> ∇₂cosχ((0,0), (pi/2, pi/3))
(-1.0, -0.0)
```
"""
function ∇₂cosχ end

for pt in 1:2
	fname = "∇$(subscripts[pt])cosχ"
	∇θ = "∂θ$(subscripts[pt])cosχ"
	∇ϕ = "∇ϕ$(subscripts[pt])cosχ"
	@eval $(Symbol(fname))(n1, n2) = ($(Symbol(∇θ))(n1, n2), $(Symbol(∇ϕ))(n1, n2))
end

end
