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

