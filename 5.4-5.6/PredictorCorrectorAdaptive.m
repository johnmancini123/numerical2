function [t, w, y] = PredictorCorrectorAdaptive(a, b, yo, f, realF, tol, hmin, hmax)
    %4th order predictor corrector Adams method
    %with variable step size
    h = hmax;  
    y = zeros(1, 100); %keeping track of y exact values 
    t = zeros(1, 100);
    w = zeros(1, 100);
    t(1) = a;
    w(1) = yo;
    y(1) = yo;
    
    for j = 1:3
        k1 = h*f(t(j), w(j));
        k2 = h*f(t(j) + h/2, w(j) + k1/2);
        k3 = h*f(t(j) + h/2, w(j) + k2/2);
        k4 = h*f(t(j) + h, w(j) + k3);
        w(j+1) = w(j) + (k1 + 2*k2 + 2*k3 + k4)/6;
        t(j+1) = t(j) + h;
        y(j+1) = abs((w(j) - realF(t(j))));
    end
  
    
    i = 5;
    
    FLAG = 1;
    LAST = 0;
    NFLAG = 1;
    
    T = t(4) + h;
    
    while(FLAG)  
        %this is how to be computationally efficient doing PC
        %here, ft0 means f(t-0*h, y_t-0*h)
        ft0 = f(t(i-1), w(i-1));
        ft1 = f(t(i-2), w(i-2));
        ft2 = f(t(i-3), w(i-3));
        %calculating wPredictor
        wPredictor = w(i-1) + (h/24)*(55*ft0 - 59*ft1 +37*ft2 - 9*f(t(i-4), w(i-4)));
        %since wCorrector is implicit, t must be iterated
        
        %calculating wCorrector
        wCorrector = w(i-1) + (h/24)*(9*f(T, wPredictor) + 19*ft0 - 5*ft1 + ft2);
        disp(wPredictor);
        fprintf("-----\n");
        disp(wCorrector);
        %determines if we keep our values or change the step size
        delta = 19*abs(wCorrector - wPredictor)/(270*h);
        
        if delta <= tol %value has been accepted
           w(i) = wCorrector;
           t(i) = T;
           y(i) = abs(w(i) -realF(t(i)));
           %%fprintf("Output coming\n");
           if (NFLAG) %if we have calculated previous 4 vals, output them
               
               for j = i-3:i
                   %fprintf("%.7f %.7f %.7f\n", t(j), w(j), y(j));
                   
               end
           else %else just output the calculated value
               %fprintf("%.7f %.7f %.7f\n", t(j), w(j), y(j));
           end
           
           if LAST == 1
               FLAG = 0;
           else
              NFLAG = 0;
              i = i+1;
          
               if delta <= .1*tol || (t(i-1)+h) > b
                   q = (tol/(2*delta))^.25;
                   if q > 4
                       h = 4*h;
                   else
                       h = q*h;
                   end
                   if h > hmax
                       h = hmax;
                   end
                   
                   if t(i-1) + 4*h > b
                       h = (b - (t(i-1)))/4;
                       LAST = 1;
                   end
                   
                   %[tempT, tempW] = rungeKutta4Helper(t(i-1), h, w(i-1), f);
                   %t = horzcat(t, tempT);
                   %w = horzcat(w, tempW);
                   %implementing rk in real time
                   for j = (i-1):(i+2)
                        k1 = h*f(t(j), w(j));
                        k2 = h*f(t(j) + h/2, w(j) + k1/2);
                        k3 = h*f(t(j) + h/2, w(j) + k2/2);
                        k4 = h*f(t(j) + h, w(j) + k3);
                        w(j+1) = w(j) + (k1 + 2*k2 + 2*k3 + k4)/6;
                        t(j+1) = t(j) + h;
                        y(j+1) = abs((w(j) - realF(t(j))));
                   end
                   
                   NFLAG = 1;
                   i = i+3;
               end
           end
               
            
        else %if result was not accepted
            q = (tol/(2*delta))^.25;
            if q < .1
                h = .1*h;
            else 
                h = q*h;
            end
            
            if h < hmin
                FLAG = 0;
                fprintf("hmin exceeded\n");
            else
                if (NFLAG)
                    i = i-3;
                end
                    %[tempT, tempW] = rungeKutta4Helper(t(i-1), h, w(i-1), f);
                    %t = horzcat(t, tempT);
                    %w = horzcat(w, tempW);
                    %implementing runge kutta 4 in real time
                    for j = (i-1):(i+2)
                        k1 = h*f(t(j), w(j));
                        k2 = h*f(t(j) + h/2, w(j) + k1/2);
                        k3 = h*f(t(j) + h/2, w(j) + k2/2);
                        k4 = h*f(t(j) + h, w(j) + k3);;
                        w(j+1) = w(j) + (k1 + 2*k2 + 2*k3 + k4)/6;
                        t(j+1) = t(j) + h;
                        y(j+1) = abs((w(j) - realF(t(j))));
                    end
                    i = i+3;
                    NFLAG = 1;
             
            end
                    
            
        end
        T = t(i-1) + h;
    end
    
end


