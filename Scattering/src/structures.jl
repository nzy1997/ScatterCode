# SimpleGraph generators:

# Universal computation by quantum walk (Andrew M. Childs)
# arXiv:0806.1972

# Fig. 1 (c)
# U = [-i 1; 1 i]
# 1 -- 5 -- 8 -- 3
#      |    |
#      6    9
#      |    |
# 2 -- 7 --10 -- 4
function andrew_basis_change(weights = ones(Float64, 10)) 
    return ScatterGraph(SimpleGraph(Edge.([(1,5), (2,7), (3, 8), (4, 10), (5,6), (5,8), (6,7), (7,10), (8,9), (9, 10)])),weights ,[1, 2],[3, 4])
end


# Fig. 1 (d)
#      ⋮
#      3    7
#      |    |
#      5 -- 6 -- 8
#      |
# 1 -- 4 -- 2
function andrew_momentum_filter(weights = ones(Float64, 7))
    return ScatterGraph(SimpleGraph(Edge.([(1,4), (2,4), (3, 5), (4, 5), (5, 6), (6,7), (6,8)])),weights, [1, 2], [3, 4])
end


# Fig. 1 (e)
#      7
#      |
#      6 -\
#      |   5 
#      4 -/
#      |
# 1 -- 3 -- 2
function andrew_momentum_separator(weights = ones(Float64, 7))
    return ScatterGraph(SimpleGraph(Edge.([(1,3), (2,3), (3, 4), (4, 5), (4,6), (5,6), (6,7)])),
        weights, [1, 2], [3,4])
end


# Universal computation by multi-particle quantum walk (Andrew M. Childs, David Gosset, Zak Webb)
# arXiv:1205.3782

# Fig. 4
function andrew_hadamard(weights = ones(Float64, 20))
    return ScatterGraph(SimpleGraph(Edge.([(1,7), (5,6),(4,6),(3,11),(6,7),(6,10),(5,9),(5,8),(7,8),(7,10),(8,9),(9,10),(8,11),(9,13),(10,11),(10,13),(8,12),(11,12),(2,12),(12,13)])),weights, [1, 2], [3, 4])
end


# Chinese character "五" (five)
# 1 -- 5 -- 3
#      |    
# 6 -- 7 -- 8
#      |    |
# 2 -- 9 -- 4
function wuzi(weights = ones(Float64, 9))
    return ScatterGraph(SimpleGraph(Edge.([(1,5), (5,7), (3, 5), (4, 9), (7,8), (4,8), (6,7), (7,9), (2,9)])),weights, [1, 2], [3, 4])
end