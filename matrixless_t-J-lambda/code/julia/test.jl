module Foo
    export makeA, makeB, model

    function makeA()
        return 1
    end
    function makeB()
        return 7
    end
    function model(x::Vector{ComplexF64}, a::Int64, b::Int64)::Vector{ComplexF64}
        return reverse(x)
    end
end
