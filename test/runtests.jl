#! /usr/bin/env julia

using Arcus
using Base.Test

import Base.isapprox
isapprox(a::Float64, b::Float64) = isapprox(a, b, atol = 2eps(typeof(a))) || isapprox(1/a, 1/b, atol = 2eps())

@testset begin

@testset "creation" begin

  @test Arc(0.0) != nothing
  @test Arc(pi/2) != nothing
  @test Arcd(44) != nothing
  @test Arc(1.0, 1.0) != nothing
  @test_throws MethodError Arc(1.0, 2.0, 3.0)

end

@testset "dedicated values" begin
  @test Arc(0.0) == zero(Arc)
  @test Float64(Arcd(90)) == pi / 2
  @test Float64(Arcd(180)) == 1pi
  @test Float64(Arcd(270)) == -pi / 2
  @test Arcd(360) == Arcd(0)
end

# @testset "arithmetic" begin
  @test -Arc(pi*3/8) == Arc(-pi*3/8)
  @test -Arc(pi*3/4) == Arc(-pi*3/4)
  @test -Arc(pi*5/4) == Arc(-pi*5/4)
  @test -Arc(pi*7/4) == Arc(-pi*7/4)
  @test -Arc(pi/3) == Arc(5pi/3)
  @test Arc(1) + Arc(1) == 2Arc(1)
  @test Arc(1) + Arc(2) ≈ Arc(3)
  @test Arc(3) - Arc(2) ≈ Arc(1)
  @test -Arc(1) + Arc(3) ≈ Arc(2)
  @test Arcd(190) / 170 ≈ Arcd(-1.0)
  @test Arc(0.3) * 30 ≈ Arc(mod(9.0, 2pi))
  @test 20Arcd(260) == Arcd(260) * 20
# end

@testset "representation as sine 0" begin
  @test reinterpret(Float64, Arc(0)) == 0
end

@testset "representation as sine $(w)" for w in linspace(-pi/4, pi/4, 11)
  @test reinterpret(Float64, Arc(w)) ≈ sin(w)
end

const logsp = logspace(-16, 306, 5) * (nextfloat(0.0)*1e16)

@testset "representation as sine $(w)" for w in logsp
  @test reinterpret(Float64, Arc(w)) == sin(w)
end

const ep = eps(2.0)
@testset "trigonometric functions $w" for w in 0:45:360
  @test sin(Arcd(w)) ≈ sind(w)
  @test cos(Arcd(w)) ≈ cosd(w)
  @test tan(Arcd(w)) ≈ tand(w)
  @test csc(Arcd(w)) ≈ cscd(w)
  @test sec(Arcd(w)) ≈ secd(w)
  @test cot(Arcd(w)) ≈ cotd(w)
end

@testset "trigonometric functions $w" for w in 2/7:2/3:2pi
  @test sin(Arc(w)) ≈ sin(w)
  @test cos(Arc(w)) ≈ cos(w)
  @test tan(Arc(w)) ≈ tan(w)
  @test csc(Arc(w)) ≈ csc(w)
  @test sec(Arc(w)) ≈ sec(w)
  @test cot(Arc(w)) ≈ cot(w)
end

end


