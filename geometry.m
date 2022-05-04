classdef geometry
properties (Access = private)
    n
    I
    k
    Ia1
    Ia2
    ca1
    ca2
    na1
    na2
end
methods
function self = geometry(n, I, k, Ia1, Ia2, ca1, ca2, na1, na2)
    self.n = n;
    self.I = I;
    self.k = k;
    self.Ia1 = Ia1;
    self.Ia2 = Ia2;
    self.ca1 = ca1;
    self.ca2 = ca2;
    self.na1 = na1;
    self.na2 = na2;
end
end
end
