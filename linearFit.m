function [a, b, da, db] = linearFit(x, y)
    n = length(x);
    if (n ~= length(y))
        error('linearFit : x and y do not have the same size');
    end
    x_mean = mean(x);
    y_mean = mean(y);
    a = 0;
    for i = 1:n
        a = a + (x(i)-x_mean)*(y(i)-y_mean);
    end
    div = 0;
    for i = 1:n
        div = div + (x(i)-x_mean)^2;
    end
    if(div == 0)
        error('linearFit : div = 0');
    end
    a = a/div;
    b = y_mean - a * x_mean;
    
    da = 0;
    for i = 1:n
        da = da + (y(i) - a*x(i) - b)^2;
    end
    if(n <= 2)
        error('linearFit : length(x) is not grater than 2');
    end
    da = da / (n-2);
    da = da / div;
    da = sqrt(da);
    da = 1.96 * da;
    
    db = 0;
    sum_x_squared = 0;
    for i = 1:n
        db = db + (y(i) - b - a * x(i))^2;
        sum_x_squared = sum_x_squared + x(i)^2;
    end
    db = db * sum_x_squared / ((n-2) * (n * sum_x_squared - (n * x_mean)^2));
    db = sqrt(db);
    db = 1.96 * db;
end