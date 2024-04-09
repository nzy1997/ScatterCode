
function plot_static_line(data;save_name::String)
    N = size(data, 2)
    ts = size(data, 1)
    @show N,ts
    
    xs = [0.1*(i-1) for i in 1:ts]
    fig = Figure()
    # Create an axis with title and labels
    ax = Axis(fig[1, 1], title = "Momentum Distribution", xlabel = "time", ylabel = L"\langle n_k\rangle")
    for i in 1:N
        lines!(ax, xs, data[:,i], label = "k=$i")
    end
    # Add a lengend at the bottom right with label size 15
    axislegend(ax; position = :rb, labelsize = 15)
    save(save_name, fig)
end

function plot_line_model(psis; framerate::Int64=100,save_name::String)
    N = size(psis, 2)
    ts = size(psis, 1)
    @show N,ts
    time = Observable(1)
    
    xs = range(1, N, length=N)
    ys_1 = @lift(psis[$time, :])
    
    
    fig=lines(xs, ys_1, color=:blue, linewidth=4,
        axis=(title=@lift("t = $(round($time, digits = 1))"),))
    
    timestamps = range(1, ts, step=1)
    
    record(fig, save_name, timestamps;
        framerate=framerate) do t
        time[] = t
    end
end

function plot_2line_model(psis1,psis2; framerate::Int64=100,save_name::String)
    N = size(psis1, 2)
    ts = size(psis1, 1)
    @show N,ts
    @assert N == size(psis2, 2)
    @assert ts == size(psis2, 1)


    time = Observable(1)
    
    xs = range(1, N, length=N)
    ys_1 = @lift(psis1[$time, :])
    ys_2 = @lift(psis2[$time, :])


    fig=lines(xs, ys_1, color=:blue, linewidth=4,
        axis=(title=@lift("t = $(round($time, digits = 1))"),))
    lines!(xs, ys_2, color=:red, linewidth=4)
    timestamps = range(1, ts, step=1)
    
    record(fig, save_name, timestamps;
        framerate=framerate) do t
        time[] = t
    end
end