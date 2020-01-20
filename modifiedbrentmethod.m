function [root,info] = modifiedbrentmethod(func, Int, params)

%values to be used
function_tol = params.func_tol;
root_tol = params.root_tol;
maxit = params.maxit;
a = Int.a;
b = Int.b;
fa = func(a);
fb = func(b); 
assert(fa * fb <= 0, "Error, root not bracketed");

%make sure f(a) < f(b) if not switch them
if abs(fa) < abs(fb)
   [b, a] = deal(a, b);
end

c = a;
total_steps = 0;
interpolation_steps = 0;
info.flag = 1;
biInt.a = a;
biInt.b = b;
root = b; 

while (abs(b - a) > root_tol) && (abs(root) > function_tol) && (total_steps < maxit) 
    
    fa = func(a);
    fb = func(b);
    fc = func(c);
    %perform inverse quadratic interpolation 
    if fa ~= fc && fb ~= fc
        root = (a * fb * fc) / ((fa - fb) * (fa - fc))...
            + (b * fa * fc) / ((fb - fa) * (fb - fc))...
            + (c * fa * fb) / ((fc - fa) * (fc - fb));
        interpolation_steps = interpolation_steps + 1;
        
        if interpolation_steps == 5 && ((abs((Int.b - Int.a) / (b - a)) < 2) &&...
           (abs((biInt.b - biInt.a) / (b - a)) < 2))
        
           root = (a + b) / 2;
           info.flag = 1;
           interpolation_steps = 0;
           %keeping track of interval created by bisection step
           biInt.a = a;
           biInt.b = b;
           if (fa * fr) < 0
              biInt.b = root;
           else 
              biInt.a = root;
           end
    
           if abs(func(biInt.a)) < abs(func(biInt.b))
              [biInt.b, biInt.a] = deal(biInt.a, biInt.b);
           end
        %modification 2. If forcing bisection within
        elseif (abs(fr) / abs(old_root)) > .5
            root = (a + b) / 2;
            info.flag = 1;
            interpolation_steps = 0;
            %keeping track of interval created by bisection step
            biInt.a = a;
            biInt.b = b;
            if (fa * fr) < 0
                biInt.b = root;
            else 
                biInt.a = root;
            end
    
            if abs(func(biInt.a)) < abs(func(biInt.b))
                [biInt.b, biInt.a] = deal(biInt.a, biInt.b);
            end
        end
   
    %if these conditions fail do secant method
    else
        root = b - ((fb * (b - a)) / (fb - fa));
        interpolation_steps = 0;
    end 
     

    c = b; 
    fr = func(root);
    
    if (fa * fr) < 0
        b = root;
    else 
        a = root;
    end
    
    if abs(fa) < abs(fb)
        [b, a] = deal(a, b);
    end
    
    total_steps = total_steps + 1;
    old_root = root;
end
if info.flag == 0
   info.flag = 1;
else 
   info.flag = 0;
end
%disp(total_steps);
%disp(biInt);
end