
using Arcus
using Base.Test

import Base.isapprox
isapprox(a::Arc, b::Arc) = abs(mod(Float64(a) - Float64(b), 2pi)) <= 10eps()

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

@testset "arithmetic" begin
  @test -Arc(pi/3) == Arc(5pi/3)
  @test Arc(1) + Arc(1) == 2Arc(1)
  @test Arc(1) + Arc(2) ≈ Arc(3)
  @test Arc(3) - Arc(2) ≈ Arc(1)
  @test -Arc(1) + Arc(3) ≈ Arc(2)
  @test Arcd(190) / 170 ≈ Arcd(-1.0)
  @test Arc(0.3) * 100 ≈ Arc(mod(30.0, 2pi))
  @test 20Arcd(260) == Arcd(260) * 20
end

@testset "representation as sine $(w)" for w in linspace(-pi/4 + 2eps(), pi/4, 11)
  @test reinterpret(Float64, Arc(w)) ≈ sin(w)
end

@testset "trigonometric functions $w" for w in linspace(0.0, 2pi-10eps(), 18)
  @test sin(Arc(w)) ≈ sin(w)
  @test cos(Arc(w)) ≈ cos(w)
  @test tan(Arc(w)) ≈ tan(w)
  @test csc(Arc(w)) ≈ csc(w)
  @test sec(Arc(w)) ≈ sec(w)
  @test cot(Arc(w)) ≈ cot(w)
end

nothing



