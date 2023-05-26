function [x_in,y_in,x_out,y_out] = cylinder(imax,r_in,r_out)

    [x_in,y_in,x_out,y_out] = deal(zeros(imax,1));
   
    for i=1:imax
        x_in(i) = r_in*cos((i-1)*(-2*pi/(imax-1)));
        y_in(i) = r_in*sin((i-1)*(-2*pi/(imax-1)));
    end
    for i=1:imax
        x_out(i) = r_out*cos((i-1)*(-2*pi/(imax-1)));
        y_out(i) = r_out*sin((i-1)*(-2*pi/(imax-1)));
    end

end