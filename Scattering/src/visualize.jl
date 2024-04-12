function animate_wave(gt::GraphWithTails, waves::AbstractVector;
        name="chains", pathname="", vertex_size=10.0, layout=:stress,
        step=ceil(Int, length(waves)/100), framerate=10)
    locs = assign_locations(gt, layout)
    edgs = [(e.src, e.dst) for e in edges(gt.graph)]
    config = LuxorGraphPlot.graphsizeconfig(locs)
    defaultunit = 50
    width, height = round(Int, config.Dx * defaultunit), round(Int, config.Dy * defaultunit)
    nframes = ceil(Int, length(waves) / step)

    function frame(scene, framenumber::Int)
        Luxor.origin(0, 0)
        vertex_sizes=vertex_size .* abs2.(waves[min(length(waves), 1+step*(framenumber-1))])
        vertex_stroke_colors=[i in gt.center ? "red" : "black" for i=1:nv(gt.graph)]
        LuxorGraphPlot._show_graph(locs, edgs; texts=fill("", nv(gt.graph)),
            vertex_sizes, vertex_stroke_colors)
    end
    mov = Movie(width, height, name)
    animate(mov, [Scene(mov, frame, 1:nframes)]; creategif=true, pathname, framerate)
end

function assign_locations(gt::GraphWithTails, layout::Symbol=:manual)
    if layout == :manual
        locs = fill((0.0, 0.0), nv(gt.graph))
        cg = center_graph(gt)
        locs0 = LuxorGraphPlot.autolocs(cg, nothing, layout, 1.0, trues(nv(cg)))
        cloc = reduce((x, y)->x.+y, locs0) ./ length(locs0)
        locs[gt.center] .= locs0
        for tail in gt.tails
            loc1 = locs[tail[1]]
            dx, dy = loc1 .- cloc
            angle = atan(dy, dx)
            for (j, v) in enumerate(tail[2:end])
                locs[v] = loc1 .+ (cos(angle), sin(angle)) .* j .* 0.03
            end
        end
    else
        return LuxorGraphPlot.autolocs(gt.graph, nothing, layout, 0.03, trues(nv(gt.graph)))
    end
end

function scatter(x, y; color="black")
    graph = SimpleGraph(length(x))
    xmin, xmax = minimum(x), maximum(x)
    ymin, ymax = minimum(y), maximum(y)
    show_graph(graph; locs=zip(x, -y), vertex_stroke_color=color, vertex_size=0.003, texts=fill("", length(x))) do transform
        setcolor("black")
        arrow(Point(transform((min(0.0, xmin), 0.0))), Point(transform((xmax+0.2*(xmax-xmin), 0.0))))
        arrow(Point(transform((0.0, max(0.0, -ymin)))), Point(transform((0.0, -ymax-0.2*(ymax-ymin)))))
    end
end