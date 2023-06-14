function [] = diagnostics(rho,u,p,E,grid)


% Run only at select iterations:
if (mod ( grid.iter, grid.Output_interval ) == 0 || grid.iter == grid.NT)

    % Clear the figure
    clf()

    % Save/ load the exact result
    if(grid.NT == grid.iter && grid.Nx == 3000)
        x_exact = grid.x;
        rho_exact = rho;
        u_exact = u;
        p_exact = p;
        E_exact = E;
        save('exact_soln','x_exact','rho_exact','u_exact','p_exact','E_exact')
    end


    %Plot the diagnostic output comparison to fig1 GA Sod
    subplot(2,3,1)
    plot(grid.x,rho,'red')
    title("Density")
    hold on



    subplot(2,3,2)
    plot(grid.x,u,'red')
    title("Velocity")
    hold on


    subplot(2,3,3)
    plot(grid.x,p,'red')
    title("Pressure")
    hold on


    subplot(2,3,4)
    plot(grid.x,E - (1/2)*(u.*u),'red')
    title("Energy")
    hold on



    if (grid.NT == grid.iter)
        load('exact_soln','x_exact','rho_exact','u_exact','p_exact','E_exact')
        %Plot the diagnostic output comparison to fig1 GA Sod
        subplot(2,3,1)
        plot(x_exact,rho_exact,":black")
        gridsize = sprintf("Grid = %d",grid.Nx);
        legend(gridsize,"Exact")

        subplot(2,3,2)
        plot(x_exact,u_exact,":black")
        legend(gridsize,"Exact",'Location','southwest')
        hold on

        subplot(2,3,3)
        plot(x_exact,p_exact,":black")
        legend(gridsize,"Exact")
        hold on

        subplot(2,3,4)
        plot(x_exact,E_exact - (1/2)*(u_exact.*u_exact),":black")
        legend(gridsize,"Exact",'Location','northwest')
        hold on
    end



    pause(0.01)

end



end