function out = Test()

a = 1;
b = 2;

y = 0;

out.x0 = 1;
out.f0 = f(out.x0);
[out.x1,rc] = csolvevb(@f,out.x0,[],1e-10,1000);
out.f1 = f(out.x1);
out.y = y;

function z = f(x)
    y = a+b*x;
    z = 2*y;
end

end
