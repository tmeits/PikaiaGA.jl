tests = ["main"] 

for t in tests
    fp = joinpath("test", "test_$t.jl")
    println("$fp ...")
    include(fp)
end

