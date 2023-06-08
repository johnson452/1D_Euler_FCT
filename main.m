%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: Grant Johnson
%Date: 5/30/2022
%Yee algorithm *1D*
%Non-relativistic Euler
%With FCT: https://pdf.sciencedirectassets.com/272570/1-s2.0-S0021999100X02859/1-s2.0-0021999179900512/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEO7%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIQDd8BdayoozwKUoM1CmwTdviKp6qo9Ia3h%2B1FiHFoEZjAIgMWTQLrWmaN8vWMbBPzvFbJrr29hxGP71T0I18rP9xYcqswUINhAFGgwwNTkwMDM1NDY4NjUiDGztTMnwtLBRAPIrdiqQBWNxqti6BKGTVdk3KpaQDWlrywAxn4Y1t%2FPjz9G%2BjjQ9qFpCyQkSryKD6coHMmTHjhAn795ND6AJdiO%2FeFt7RcWk77qP8LtqUtw8YNl03K78cXpaxa7n4emZhjNwBo7jgk7%2BUx6WcGgtz21rpDQ2uR%2FCSmunBzvxMfJwT2VKOVtcEbQ3RErWaEJEYxMrt%2BKsQU2fdDO2dPev%2BPt9khy%2FB0pTCUZ%2BpK1ZVdDUO6nc41M11LWKypeAKjIKmdyEeADL17xHeA8cWDG7f%2BSfIwDXv9491WKcEZs5NkdoLEsTv5dRVicgsrTNWSmESmsgxumg5JwlfWofYs5BtJaKDkGXitiOldeKfKFL1Tm%2Bmr1mTB3H%2BehapZu2LsUZm62zFnZUDifQDrkhWStBK9jCyPzCXFcei46gozWlf%2Bbf7Q%2FS7EZ7CkN3mLjhPGA0AyUey7ZWX6vT0tsq2W3Cozv1q1SRtb4MuGvh8ONo1C1cWeuKuG88uVSuxFMVXqKXsee0qKi5fqmNl2VukJnbeetv6NhUSdRfX5qG5AvD2%2BYAtu2i6l9xxE0hxsdQWbyyw25eqm9dDadx5CtFREoWhJO9EJ5tsUgL039s6mQ6ISzA5x6r%2BUVIRl38u8YG4Vb%2FxctTM6jTPXy2b3H3n4ytK5JT%2F4bKSrZEgQkx569kFgyaxYsEfEqewDIWyG7sgLgU%2BtBY1e3hU8dgNEwmYEAQb0P3wtQiowtrz1KEcTnWO6VPs0IzX5Nc8pRDTOXL2YMDh5fes4XEpilagy2XaG%2FSN%2BFFo2HCamTBYNDCgauySdTuokJgzZ0UzfGf8QpBlwBnNRjVWqh4ePsd3AhQ73pO%2FVFt5NZbyKpM%2BCG3UZXV9SL77XBdYvANMNSe%2BaMGOrEB%2B9j9Yp39Yy7E3R8K9cV2ECj04VPkTs9gD4a%2FeOD7HthXoR4JgdqjE5qrIj71M%2FZ37pKV5OaxGlX6j790LUU0CfaXRvq6akN8b77yWYMk6Bk%2FF6U9PH2tZfYx%2BbevHsGfrkuR%2BHgHdpMj0Bw1IqeodlfF4DuQDPxCJCgJ5KFFbQ%2BCXl8dkmNAvn7DO%2Fdai%2F%2B7vbaqDXm%2Bfgauhp4MaI39%2B1BTgKwmDsozXBPrKGkChubk&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20230605T221714Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYX6UT6RNP%2F20230605%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=58e694cc967a550e40e4f6ce95030561576e66daaa1b80a263f620f6f0ed4e74&hash=0fa4626b4a6ef5040dee43c261ea3ac495898f88afc2b248f4f5b1725b9dc324&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=0021999179900512&tid=spdf-6248abf5-0a3f-43c2-950a-232e37cf6e4c&sid=da0f7889150d13473f2aef1-529c91f6051cgxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=131559025c0107570755&rr=7d2bc8d61e3630b8&cc=us

%Notes:
%-1D
% Fluid Scheme: https://ammar-hakim.org/sj/hancock-muscl.html
% Good example: https://en.wikipedia.org/wiki/MUSCL_scheme
% https://www.cambridge.org/core/services/aop-cambridge-core/content/view/8F5CD408E7073099BFDE1E409C1E79AB/S0022112077001463a.pdf/a-numerical-study-of-a-converging-cylindrical-shock.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

%Setup Initial shock-tube problem 
[rho,u,p,E,grid] = make_grid();

%Make the diagnostic Figure
figure('units','normalized','outerposition',[0 0 0.5 0.5])

%%% Time loop %%%
while(grid.time < grid.t_max)
    
    %Call i/o and diagnostics
    diagnostics(rho,u,p,E,grid);
    
    %Update the gridtime
    grid.time = grid.time + grid.dt;
    
    %Update the iterator
    grid.iter = grid.iter + 1;
 
    %Updater - updates all quantities simultaneosly
    % n -> n + 1 all quantities
    [rho,u,p,E,grid] = push_all(rho,u,p,E,grid);

     %BC - All outflow (copy)
    [rho,u,p,E,grid] = BC(rho,u,p,E,grid);

end
%%% End Time Loop %%%
%%% End main %%%