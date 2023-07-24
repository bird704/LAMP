% Feed rate

function F_in = feed_rate(t, split)
% split(1) = batch(gly) -> fed-batch(gly)
% split(2) = fed-batch(gly) -> fed-batch(meth)
% split(3) = fed-batch(meth) -> constant-fed(meth)
% F_in -> mL/h

if (t >= split(1)) && (t < split(2))
    F_in = 0;

elseif (t >= split(2)) && (t < split(3))
    F_in = 0.0385*exp(0.18*(t-split(2)));

elseif (t >= split(3)) && (t < split(4))
    F_in = 0.0105*exp(0.42*(t-split(3)));

else
    F_in = 0.0105*exp(0.42*(split(4)-split(3)));

end
