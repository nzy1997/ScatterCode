using YaoXYModel

n=10
k=3
time_stop=500
@info "parameters: n=$n, k=$k, time_stop=$time_stop."

for Ω in [0.2, 0.5]
    @info "Running Yao simulation for Ω=$Ω."
    disYao =sim_model(;n,k,time_stop,Ω)

    @info "Plotting results for Ω=$Ω. Results are saved in examples/."
    plot_static_line(disYao[:,1:n];save_name="examples/momentum_dis_omega$Ω.png")
    plot_line_model(disYao[:,1:n];save_name="examples/momentum_dis_omega$Ω.mp4")
    plot_line_model(disYao[:,n+1:2*n];save_name="examples/expectz_dis_omega$Ω.mp4")

    @info "Running two level approximate simulation for Ω=$Ω."
    dis2 =twolevel_sim_model(;n,k,time_stop,Ω)

    @info "Plotting results for Ω=$Ω. Results are saved in examples/."
    plot_static_line(dis2;save_name="examples/momentum_dis_two_omega$Ω.png")
    plot_line_model(dis2;save_name="examples/momentum_dis_two_omega$Ω.mp4")

    @info "Plotting Yao vs two level approximate results for Ω=$Ω. Results are saved in examples/."
    plot_2line_model(disYao[:,1:n],dis2;save_name="examples/Yaovstwo_omega$Ω.mp4")
end