using Test
using TwoPointFunctions

import TwoPointFunctions: ∂θ₁∂θ₂cosχ, ∂θ₂∂θ₁cosχ, ∂ϕ₁∂ϕ₂cosχ, ∂ϕ₂∂ϕ₁cosχ, 
∂ϕ₂∂ϕ₂cosχ, ∂ϕ₁∂ϕ₁cosχ, ∂θ₁∂ϕ₁cosχ, ∂θ₂∂ϕ₂cosχ, ∂ϕ₁∂θ₁cosχ, ∂ϕ₂∂θ₂cosχ

isapproxzero(x) = isapprox(x, zero(x), atol=1e-14, rtol=1e-8)

@testset "cosχ" begin
	# cos(θ₁)cos(θ₂) + sin(θ₁)sin(θ₂)cos(ϕ₁ - ϕ₂)
    @test cosχ((0, 0), (0, 0)) == 1
    @test cosχ((0, 0), (0, pi/3)) == 1
    @test cosχ((0, 0), (pi, 0)) == -1
    @test cosχ((0, 0), (pi, pi/3)) == -1
    
    @test cosχ((0, 0), (0,0)) == 1
    @test cosχ((0, pi/3), (0,0)) == 1
    @test cosχ((pi, 0), (0,0)) == -1
    @test cosχ((pi, pi/3), (0,0)) == -1
    
    @test isapproxzero(cosχ((0, 0), (pi/2, 0)))
    @test isapproxzero(cosχ((0, 0), (pi/2, pi/3)))
    
    @test isapproxzero(cosχ((pi/2, 0), (0, 0)))
    @test isapproxzero(cosχ((pi/2, pi/3), (0, 0)))

    @test cosχ((0, 0), (pi/3, 0)) ≈ cos(pi/3)
    @test cosχ((0, 0), (pi/3, pi/3)) ≈ cos(pi/3)
end

@testset "first derivatives" begin
    @testset "∂θ₁cosχ" begin
		# -sin(θ₁)cos(θ₂) + cos(θ₁)sin(θ₂)cos(ϕ₁ - ϕ₂)
	    @test isapproxzero(∂θ₁cosχ((0, 0), (0, 0)))
	    @test isapproxzero(∂θ₁cosχ((0, 0), (0, pi/3)))
	    @test isapproxzero(∂θ₁cosχ((0, 0), (pi, 0)))
	    @test isapproxzero(∂θ₁cosχ((0, 0), (pi, pi/3)))
	    @test ∂θ₁cosχ((0, 0), (pi/2, 0)) ≈ 1
	    @test ∂θ₁cosχ((0, 0), (pi/2, pi/3)) ≈ 0.5
	    @test ∂θ₁cosχ((pi/2, 0), (0, 0)) ≈ -1
	    @test ∂θ₁cosχ((pi/2, pi/3), (0, 0)) ≈ -1
	end

	@testset "∂θ₂cosχ" begin
		# -cos(θ₁)sin(θ₂) + sin(θ₁)cos(θ₂)cos(ϕ₁ - ϕ₂)
	    @test isapproxzero(∂θ₂cosχ((0, 0), (0, 0)))
	    @test isapproxzero(∂θ₂cosχ((0, 0), (0, pi/3)))
	    @test isapproxzero(∂θ₂cosχ((0, 0), (pi, 0)))
	    @test isapproxzero(∂θ₂cosχ((0, 0), (pi, pi/3)))
	    @test ∂θ₂cosχ((0, 0), (pi/2, 0)) ≈ -1
	    @test ∂θ₂cosχ((0, 0), (pi/2, pi/3)) ≈ -1
	    @test ∂θ₂cosχ((pi/2, 0), (0, 0)) ≈ 1
	    @test ∂θ₂cosχ((pi/2, pi/3), (0, 0)) ≈ 0.5
	end

	@testset "∂ϕ₁cosχ" begin
		# -sin(θ₁)sin(θ₂)sin(ϕ₁ - ϕ₂)
	    @test isapproxzero(∂ϕ₁cosχ((0, 0), (0, 0)))
	    @test isapproxzero(∂ϕ₁cosχ((0, 0), (0, pi/3)))
	    @test isapproxzero(∂ϕ₁cosχ((0, 0), (pi, 0)))
	    @test isapproxzero(∂ϕ₁cosχ((0, 0), (pi, pi/3)))
	    
	    @test isapproxzero(∂ϕ₁cosχ((0, 0), (0, 0)))
	    @test isapproxzero(∂ϕ₁cosχ((0, pi/3), (0, 0)))
	    @test isapproxzero(∂ϕ₁cosχ((pi, 0), (0, 0)))
	    @test isapproxzero(∂ϕ₁cosχ((pi, pi/3), (0, 0)))
	    
	    @test isapproxzero(∂ϕ₁cosχ((pi/2, 0), (pi/2, 0)))
	    @test ∂ϕ₁cosχ((pi/2, 0), (pi/2, pi/6)) ≈ 0.5
	    @test ∂ϕ₁cosχ((pi/2, pi/6), (pi/2, 0)) ≈ -0.5
	    @test ∂ϕ₁cosχ((pi/2, pi/2), (pi/2, 0)) ≈ -1
	end

	@testset "∇ϕ₁cosχ" begin
	    # -sin(θ₂)sin(ϕ₁ - ϕ₂)
	    @test isapprox(∇ϕ₁cosχ((0, 0), (0, 0)), 0, atol=1e-14)
	    @test isapprox(∇ϕ₁cosχ((0, 0), (0, pi/3)), 0, atol=1e-14)
	    @test isapprox(∇ϕ₁cosχ((0, 0), (pi, 0)), 0, atol=1e-14)
	    @test isapprox(∇ϕ₁cosχ((0, 0), (pi, pi/3)), 0, atol=1e-14)
	    
	    @test isapprox(∇ϕ₁cosχ((0, 0), (0, 0)), 0, atol=1e-14)
	    @test isapprox(∇ϕ₁cosχ((0, pi/3), (0, 0)), 0, atol=1e-14)
	    @test isapprox(∇ϕ₁cosχ((pi, 0), (0, 0)), 0, atol=1e-14)
	    @test isapprox(∇ϕ₁cosχ((pi, pi/3), (0, 0)), 0, atol=1e-14)
	    
	    @test isapprox(∇ϕ₁cosχ((pi/2, 0), (pi/2, 0)), 0, atol=1e-14)
	    @test ∇ϕ₁cosχ((pi/2, 0), (pi/2, pi/6)) ≈ 0.5
	    @test ∇ϕ₁cosχ((pi/2, pi/6), (pi/2, 0)) ≈ -0.5
	    @test ∇ϕ₁cosχ((pi/2, pi/2), (pi/2, 0)) ≈ -1
	end

	@testset "∂ϕ₂cosχ" begin
		# sin(θ₁)sin(θ₂)sin(ϕ₁ - ϕ₂)
	    @test isapprox(∂ϕ₂cosχ((0, 0), (0, 0)), 0, atol=1e-14)
	    @test isapprox(∂ϕ₂cosχ((0, 0), (0, pi/3)), 0, atol=1e-14)
	    @test isapprox(∂ϕ₂cosχ((0, 0), (pi, 0)), 0, atol=1e-14)
	    @test isapprox(∂ϕ₂cosχ((0, 0), (pi, pi/3)), 0, atol=1e-14)
	    
	    @test isapprox(∂ϕ₂cosχ((0, 0), (0, 0)), 0, atol=1e-14)
	    @test isapprox(∂ϕ₂cosχ((0, pi/3), (0, 0)), 0, atol=1e-14)
	    @test isapprox(∂ϕ₂cosχ((pi, 0), (0, 0)), 0, atol=1e-14)
	    @test isapprox(∂ϕ₂cosχ((pi, pi/3), (0, 0)), 0, atol=1e-14)
	    
	    @test isapprox(∂ϕ₂cosχ((pi/2, 0), (pi/2, 0)), 0, atol=1e-14)
	    @test ∂ϕ₂cosχ((pi/2, 0), (pi/2, pi/6)) ≈ -0.5
	    @test ∂ϕ₂cosχ((pi/2, pi/6), (pi/2, 0)) ≈ 0.5
	    @test ∂ϕ₂cosχ((pi/2, pi/2), (pi/2, 0)) ≈ 1
	end

	@testset "∇ϕ₂cosχ" begin
	    # sin(θ₁)sin(ϕ₁ - ϕ₂)
	    @test isapprox(∇ϕ₂cosχ((0, 0), (0, 0)), 0, atol=1e-14)
	    @test isapprox(∇ϕ₂cosχ((0, 0), (0, pi/3)), 0, atol=1e-14)
	    @test isapprox(∇ϕ₂cosχ((0, 0), (pi, 0)), 0, atol=1e-14)
	    @test isapprox(∇ϕ₂cosχ((0, 0), (pi, pi/3)), 0, atol=1e-14)
	    
	    @test isapprox(∇ϕ₂cosχ((0, 0), (0, 0)), 0, atol=1e-14)
	    @test isapprox(∇ϕ₂cosχ((0, pi/3), (0, 0)), 0, atol=1e-14)
	    @test isapprox(∇ϕ₂cosχ((pi, 0), (0, 0)), 0, atol=1e-14)
	    @test isapprox(∇ϕ₂cosχ((pi, pi/3), (0, 0)), 0, atol=1e-14)
	    
	    @test isapprox(∇ϕ₂cosχ((pi/2, 0), (pi/2, 0)), 0, atol=1e-14)
	    @test ∇ϕ₂cosχ((pi/2, 0), (pi/2, pi/6)) ≈ -0.5
	    @test ∇ϕ₂cosχ((pi/2, pi/6), (pi/2, 0)) ≈ 0.5
	    @test ∇ϕ₂cosχ((pi/2, pi/2), (pi/2, 0)) ≈ 1
	end
end

@testset "second derivatives" begin
    @testset "∂²θ₁cosχ" begin
        # -cos(θ₁)cos(θ₂) -sin(θ₁)sin(θ₂)cos(ϕ₁ - ϕ₂) == -cosχ
        @test isapproxzero(∂²θ₁cosχ((pi/2, 0), (0, 0)))
        @test ∂²θ₁cosχ((pi/2, 0), (pi/2, 0)) == -1
        @test ∂²θ₁cosχ((pi/2, 0), (pi/2, pi/3)) ≈ -0.5

        @test ∂²θ₁cosχ((0, 0), (0, 0)) ≈ -cosχ((0, 0), (0, 0))
        @test ∂²θ₁cosχ((0, 0), (pi, 0)) ≈ -cosχ((0, 0), (pi, 0))
        @test ∂²θ₁cosχ((0, 0), (pi/2, 0)) ≈ -cosχ((0, 0), (pi/2, 0))
        @test ∂²θ₁cosχ((0, pi/3), (pi/2, 0)) ≈ -cosχ((0, pi/3), (pi/2, 0))
    end
    
    @testset "∂²θ₂cosχ" begin
        # -cos(θ₁)cos(θ₂) - sin(θ₁)sin(θ₂)cos(ϕ₁ - ϕ₂) = -cosχ
        @test isapproxzero(∂²θ₂cosχ((pi/2, 0), (0, 0)))
        @test ∂²θ₂cosχ((pi/2, 0), (pi/2, 0)) == -1
        @test ∂²θ₂cosχ((pi/2, 0), (pi/2, pi/3)) ≈ -0.5

        @test ∂²θ₂cosχ((0, 0), (0, 0)) ≈ -cosχ((0, 0), (0, 0))
        @test ∂²θ₂cosχ((0, 0), (pi, 0)) ≈ -cosχ((0, 0), (pi, 0))
        @test ∂²θ₂cosχ((0, 0), (pi/2, 0)) ≈ -cosχ((0, 0), (pi/2, 0))
        @test ∂²θ₂cosχ((0, pi/3), (pi/2, 0)) ≈ -cosχ((0, pi/3), (pi/2, 0))
    end
    
    @testset "∂θ₁∂θ₂cosχ" begin
    	# sin(θ₁)sin(θ₂) + cos(θ₁)cos(θ₂)cos(ϕ₁ - ϕ₂)

    	@test ∂θ₁∂θ₂cosχ((0, 0), (0, 0)) ≈ 1
    	@test ∂θ₁∂θ₂cosχ((0, pi/3), (0, 0)) ≈ cos(pi/3)
    	@test ∂θ₁∂θ₂cosχ((0, pi/3), (pi, 0)) ≈ -cos(pi/3)

    	@test isapproxzero(∂θ₁∂θ₂cosχ((pi/2, pi/3), (0, 0)))
    	@test isapproxzero(∂θ₁∂θ₂cosχ((pi/2, pi/3), (0, pi/3)))
    	@test isapproxzero(∂θ₁∂θ₂cosχ((pi/2, pi/3), (pi, 0)))
    	@test isapproxzero(∂θ₁∂θ₂cosχ((pi/2, 0), (pi, 0)))
    	@test isapproxzero(∂θ₁∂θ₂cosχ((pi/2, 0), (pi, pi/3)))

    	@test ∂θ₁∂θ₂cosχ((pi/2, 0), (pi/2, 0)) ≈ 1
    	@test ∂θ₁∂θ₂cosχ((pi/2, 0), (pi/2, pi/3)) ≈ 1
        
        @testset "∂θ₂∂θ₁cosχ" begin 
        	for θ₁ in LinRange(0,π,10), ϕ₁ in LinRange(0,2π,20), 
            	θ₂ in LinRange(0,π,10), ϕ₂ in LinRange(0,2π,20)

            	@test isapprox(∂θ₂∂θ₁cosχ((θ₁, ϕ₁), (θ₂, ϕ₂)) , ∂θ₁∂θ₂cosχ((θ₁, ϕ₁), (θ₂, ϕ₂)), atol=1e-14, rtol=1e-8)
        	end
        end
    end

    @testset "∂²ϕ₁cosχ" begin
        # -sin(θ₁)sin(θ₂)cos(ϕ₁ - ϕ₂)

        @test isapproxzero(∂²ϕ₁cosχ((0, 0), (0, 0)))
        @test isapproxzero(∂²ϕ₁cosχ((0, 0), (pi,0)))
        @test isapproxzero(∂²ϕ₁cosχ((pi, π/4), (pi,0)))
        @test isapproxzero(∂²ϕ₁cosχ((0, pi/3), (0, 0)))
        @test isapproxzero(∂²ϕ₁cosχ((0, pi/3), (0, pi/2)))
        @test isapproxzero(∂²ϕ₁cosχ((pi/4, pi/2), (pi/4, 0)))
        
        @test ∂²ϕ₁cosχ((pi/4, 0), (pi/4, 0)) ≈ -sin(π/4)^2
        @test ∂²ϕ₁cosχ((pi/4, pi/6), (pi/4, pi/6)) ≈ -sin(π/4)^2

        @test ∂²ϕ₁cosχ((pi/4, pi/3), (pi/4, pi/6)) ≈ -sin(π/4)^2*cos(π/6)

        @testset "∂ϕ₁∂ϕ₁cosχ" begin 
        	for θ₁ in LinRange(0,π,10), ϕ₁ in LinRange(0,2π,20), 
            	θ₂ in LinRange(0,π,10), ϕ₂ in LinRange(0,2π,20)

            	@test isapprox(∂²ϕ₁cosχ((θ₁, ϕ₁), (θ₂, ϕ₂)) , ∂ϕ₁∂ϕ₁cosχ((θ₁, ϕ₁), (θ₂, ϕ₂)), atol=1e-14, rtol=1e-8)
        	end
        end
    end
    
    @testset "∂²ϕ₂cosχ" begin
        # -sin(θ₁)sin(θ₂)cos(ϕ₁ - ϕ₂)

        @test isapproxzero(∂²ϕ₂cosχ((0, 0), (0, 0)))
        @test isapproxzero(∂²ϕ₂cosχ((0, 0), (pi,0)))
        @test isapproxzero(∂²ϕ₂cosχ((pi, π/4), (pi,0)))
        @test isapproxzero(∂²ϕ₂cosχ((0, pi/3), (0, 0)))
        @test isapproxzero(∂²ϕ₂cosχ((0, pi/3), (0, pi/2)))
        @test isapproxzero(∂²ϕ₂cosχ((pi/4, pi/2), (pi/4, 0)))
        
        @test ∂²ϕ₂cosχ((pi/4, 0), (pi/4, 0)) ≈ -sin(π/4)^2
        @test ∂²ϕ₂cosχ((pi/4, pi/6), (pi/4, pi/6)) ≈ -sin(π/4)^2

        @test ∂²ϕ₂cosχ((pi/4, pi/3), (pi/4, pi/6)) ≈ -sin(π/4)^2*cos(π/6)

        @testset "∂ϕ₂∂ϕ₂cosχ" begin 
        	for θ₁ in LinRange(0,π,10), ϕ₁ in LinRange(0,2π,20), 
            	θ₂ in LinRange(0,π,10), ϕ₂ in LinRange(0,2π,20)

            	@test isapprox(∂²ϕ₂cosχ((θ₁, ϕ₁), (θ₂, ϕ₂)) , ∂ϕ₂∂ϕ₂cosχ((θ₁, ϕ₁), (θ₂, ϕ₂)), atol=1e-14, rtol=1e-8)
        	end
        end
    end

    @testset "∂ϕ₁∂ϕ₂cosχ" begin
        # sin(θ₁)sin(θ₂)cos(ϕ₁ - ϕ₂)

        @test isapproxzero(∂ϕ₁∂ϕ₂cosχ((0, 0), (0, 0)))
        @test isapproxzero(∂ϕ₁∂ϕ₂cosχ((0, 0), (pi,0)))
        @test isapproxzero(∂ϕ₁∂ϕ₂cosχ((pi, π/4), (pi,0)))
        @test isapproxzero(∂ϕ₁∂ϕ₂cosχ((0, pi/3), (0, 0)))
        @test isapproxzero(∂ϕ₁∂ϕ₂cosχ((0, pi/3), (0, pi/2)))
        @test isapproxzero(∂ϕ₁∂ϕ₂cosχ((pi/4, pi/2), (pi/4, 0)))
        
        @test ∂ϕ₁∂ϕ₂cosχ((pi/4, 0), (pi/4, 0)) ≈ sin(π/4)^2
        @test ∂ϕ₁∂ϕ₂cosχ((pi/4, pi/6), (pi/4, pi/6)) ≈ sin(π/4)^2

        @test ∂ϕ₁∂ϕ₂cosχ((pi/4, pi/3), (pi/4, pi/6)) ≈ sin(π/4)^2*cos(π/6)

        @testset "∂ϕ₂∂ϕ₁cosχ" begin 
        	for θ₁ in LinRange(0,π,10), ϕ₁ in LinRange(0,2π,20), 
            	θ₂ in LinRange(0,π,10), ϕ₂ in LinRange(0,2π,20)

            	@test isapprox(∂ϕ₂∂ϕ₁cosχ((θ₁, ϕ₁), (θ₂, ϕ₂)) , ∂ϕ₁∂ϕ₂cosχ((θ₁, ϕ₁), (θ₂, ϕ₂)), atol=1e-14, rtol=1e-8)
        	end
        end
    end

    @testset "∂θ₁∂ϕ₁cosχ" begin

        @testset "∂ϕ₁∂θ₁cosχ" begin 
        	for θ₁ in LinRange(0,π,10), ϕ₁ in LinRange(0,2π,20), 
            	θ₂ in LinRange(0,π,10), ϕ₂ in LinRange(0,2π,20)

            	@test isapprox(∂θ₁∂ϕ₁cosχ((θ₁, ϕ₁), (θ₂, ϕ₂)) , ∂ϕ₁∂θ₁cosχ((θ₁, ϕ₁), (θ₂, ϕ₂)), atol=1e-14, rtol=1e-8)
        	end
        end
    end

    @testset "∂θ₂∂ϕ₂cosχ" begin

        @testset "∂ϕ₂∂θ₂cosχ" begin 
        	for θ₁ in LinRange(0,π,10), ϕ₁ in LinRange(0,2π,20), 
            	θ₂ in LinRange(0,π,10), ϕ₂ in LinRange(0,2π,20)

            	@test isapprox(∂θ₂∂ϕ₂cosχ((θ₁, ϕ₁), (θ₂, ϕ₂)) , ∂ϕ₂∂θ₂cosχ((θ₁, ϕ₁), (θ₂, ϕ₂)), atol=1e-14, rtol=1e-8)
        	end
        end
    end
end