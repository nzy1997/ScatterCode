# credit: https://github.com/toivoh/ConicHulls.jl/blob/888109ab6fe7d334ffa48e649947f0c293a967ae/src/Common.jl#L37
function isevenperm(p)
    @assert isperm(p)
    n = length(p)
    used = falses(n)
    even = true
    for k = 1:n
        if used[k]; continue; end
        # Each even cycle flips even (an odd number of times)
        used[k] = true
        j = p[k]
        while !used[j]
            used[j] = true
            j = p[j]
            even = !even
        end
    end
    even
end

# upper_triangle_count and its inverse
upper_triangle_count(j::Int, n::Int) = j

function upper_triangle_count(i::Int, j::Int, n::Int)
	@assert 1 <= i < j <= n
	return (2 * n - i) * (i - 1) ÷ 2 + j - i
end

function upper_triangle_count_inv(k::Int, n::Int)
	@assert 1 <= k <= matrix_size(4,n)
	i = Int(ceil((2 * n - 1 - sqrt((2 * n - 1)^2 - 8 * k)) / 2))
	j = k + i - (2 * n - i) * (i - 1) ÷ 2
	return i, j
end